"""
auther: rongbupt
espprc 求解使用C++实现 https://github.com/DouYishun/vrp-espprc
分支定界框架使用pybnb
"""
import pandas as pd
from ortools.linear_solver import pywraplp
import time
import random
import copy
from collections import defaultdict
from espprc import calculate
import pybnb
random.seed(1)


class VRPTW(pybnb.Problem):
    def generate_distance(self, file_path):
        """
        读取算例文件
        """
        data = pd.read_csv(file_path)

        def disCalcu(i, j):
            node1 = coordinates[i]
            node2 = coordinates[j]
            result = ((node1[0] - node2[0]) ** 2 + (node1[1] - node2[1]) ** 2) ** 0.5
            return round(result, 2)
        
        coordinates = data['coordinates']  # 点坐标
        coordinates = coordinates.map(lambda x: (float(x.split(',')[0][1:]), float(x.split(',')[1][1:-1]))).to_dict()
        demand = data['demand'].to_dict()  # 需求量
        demand[len(demand)] = demand[0]
        ready_time = data['ready_time'].to_dict()  # 最早开始时间
        ready_time[len(ready_time)] = ready_time[0]
        due_time = data['due_time'].to_dict()  # 最晚开始时间
        due_time[len(due_time)] = due_time[0]
        service_time = data['service_time'].to_dict()  # 服务时长
        service_time[len(service_time)] = service_time[0]

        dis_mat = []  # 距离矩阵
        for i in range(len(data)):
            dis_mat_row = []
            for j in range(len(data)):
                dis_mat_row.append(disCalcu(i, j))
            dis_mat.append(dis_mat_row)

        self.numnodes = len(coordinates)  # 点的总数,包含0点一次

        dis_mat.append(copy.deepcopy(dis_mat[0]))
        for i in dis_mat:
            i.append(copy.deepcopy(i[0]))
        return dis_mat

    def generate_initial_solution(self):
        """
        生成初始路径池
        """
        path_list = []
        for i in range(1, self.numnodes):
            path_list.append([i])
        return path_list


    def __init__(self, file_path, file_path_cpp):
        self.file_path_cpp = file_path_cpp
        self.objective_value = float('inf')  # 记录当前状态的目标函数值
        self.objective_routes = []  # 记录当前解的路径方案
        self._bound = float('-inf')  # 记录当前状态的目标函数值下界
        self.distance = self.generate_distance(file_path)
        self.largenum = 1e10

        self.ini_continuous_routes_pool = []  # 初始路径池
        self.ini_ban_edges = []  # 去除的边
        self.ini_set_edges = []  # 保留的边
        self.edge_handled = []  # 已经处理过的边

        self.solved_continuous_routes_pool = self.generate_initial_solution()
        self.solved_ban_edges = []
        self.solved_set_edges = []
        self.solved_continuous_coeff = []



    def get_pi(self, master):
        """
        获取对偶变量
        """
        pi = [c.dual_value() for c in master.constraints()]
        pi.insert(0, 0)
        pi.append(0)
        return pi

    def get_route_length(self, route, numnodes, distance):
        """
        计算路径长度
        """
        ksit = [0] * numnodes
        lastcustomer = 0
        route_length = 0
        for i, cust in enumerate(route):
            ksit[cust] += 1
            route_length += distance[lastcustomer][cust]
            lastcustomer = cust
        route_length += distance[lastcustomer][0]
        return route_length, ksit

    def get_reduced_cost_mat(self, pi, distance):
        """
        计算路径reduced_cost
        """
        reduced_cost_mat = []
        for i in range(len(pi)):
            reduced_cost_mat_temp = []
            for j in range(len(pi)):
                if i == j:
                    reduced_cost_mat_temp.append(0)
                else:
                    reduced_cost_mat_temp.append(distance[i][j] - pi[i] / 2 - pi[j] / 2)
            reduced_cost_mat.append(reduced_cost_mat_temp)
        return reduced_cost_mat


    def SUBP(self, pi, distance):
        """
        子问题求解，返回新路径和其reduced_cost
        """
        ban_flat = []
        set_flat = []
        for i in self.ini_ban_edges:
            ban_flat.append(i[0])
            ban_flat.append(i[1])
        for i in self.ini_set_edges:
            set_flat.append(i[0])
            set_flat.append(i[1])
        new_path = calculate(self.file_path_cpp, pi[:-1], ban_flat, set_flat)  # 调用C++库的calculate方法
        if new_path == ():
            return [], 0
        new_path = [i for i in new_path]
        espprc_value = 0
        for i in range(len(new_path) - 1):
            espprc_value += distance[new_path[i]][new_path[i + 1]] - pi[new_path[i]] / 2 - pi[new_path[i + 1]] / 2
        espprc_value += distance[new_path[-1]][0] - pi[new_path[-1]] / 2
        return new_path, espprc_value

    def update_distance(self):
        """
        对于去除的边，将其长度设置为一个很大的数，
        对于保留的边(a,b),把(a,x)和(x,b)设置为一个很大的数
        """
        self.distance_local = copy.deepcopy(self.distance)
        for ban_edge in self.ini_ban_edges:
            self.distance_local[ban_edge[0]][ban_edge[1]] = self.largenum
        for set_edge in self.ini_set_edges:
            self.ini_continuous_routes_pool.append([i for i in set_edge])
            for i in range(len(self.distance)):
                if i != set_edge[0] and set_edge[1] != 0:
                    self.distance_local[i][set_edge[1]] = self.largenum
                if i != set_edge[1] and set_edge[0] != 0:
                    self.distance_local[set_edge[0]][i] = self.largenum

    def sense(self):
        return pybnb.minimize

    def objective(self):
        return self.objective_value

    def bound(self):
        """
        节点求解的列生成过程
        """
        sub_objective_list = []
        self.update_distance()
        distance = self.distance_local  # self.distance_local是根据分支修改后的距离矩阵
        routes_pool = self.ini_continuous_routes_pool
        initialroutecount = len(routes_pool)
        col_coeff = defaultdict(list)  # 路径对应的列系数
        route_length = defaultdict(float)
        # 构造主问题
        for r in range(initialroutecount):
            route_length[r], col_coeff[r] = self.get_route_length(routes_pool[r], self.numnodes, distance)
        master = pywraplp.Solver.CreateSolver('CLP_LINEAR_PROGRAMMING')
        y = {}
        for r in range(initialroutecount):
            y[r] = master.NumVar(lb=0.0, ub=float('inf'), name='y_%s' % r)
        custconstr = {}
        for i in range(1, self.numnodes):
            custconstr[i] = master.Add(
                sum(col_coeff[r][i] * y[r] for r in range(initialroutecount)) - 1 >= 0, 'cust_%s' % i)
        master.Minimize(sum(route_length[r] * y[r] for r in range(initialroutecount)))
        master.Solve()
        iteration_cout = 0
        start_time = time.time()
        while True:
            pi = self.get_pi(master)
            result_temp, sub_objective = self.SUBP(pi, distance)
            new_path = [result_temp[1:-1]]
            sub_objective_list.append(sub_objective)
            # print("new_path:", new_path, "sub_objective:", sub_objective)
            if sub_objective >= 0:
                # print("没有负值新路")
                break
            if len(sub_objective_list) > 20 and sub_objective == sub_objective_list[-20]:
                # print("已经连续20次生成重复新路")
                break

            # 把新生成的路径加入到主问题
            K = len(routes_pool)
            for i, j in enumerate(new_path):
                route_length[K], col_coeff[K] = self.get_route_length(j, self.numnodes, distance)
                routes_pool.append(j)
                y[K] = master.NumVar(lb=0.0, ub=float('inf'), name='y_%s' % K)  # 增加新列对应的变量
                for i in range(1, self.numnodes):  # 增加约束新列
                    custconstr[i].SetCoefficient(y[K], col_coeff[K][i])
                master.Objective().SetCoefficient(y[K], route_length[K])
                K += 1
            master.Solve()
            # print('第{}次列生成的RMP目标值'.format(iteration_cout), "  ", master.Objective().Value())
            iteration_cout += 1
        end_time = time.time()
        # print("当前点列生成时长----", end_time - start_time)
        self._bound = master.Objective().Value()  # RMP连续解，就是当前节点的下界
        # print("连续解:    ", self._bound)
        basic_index = []
        for i in y:
            if (y[i].solution_value() != 0):
                basic_index.append(i)
        self.solved_continuous_routes_pool = [routes_pool[i] for i in basic_index]  # 当前状态的连续解路径，该项需要传递到子节点
        self.solved_continuous_coeff = [y[i].solution_value() for i in basic_index]  # 当前状态的连续解路径对应的系数，用于branch
        # print("该节点求解结束时的连续路径：", self.solved_continuous_routes_pool)

        # 把主问题的整数约束加上，求解整数规划的解
        NumberofVariable = len(routes_pool)
        master = pywraplp.Solver.CreateSolver('CBC_MIXED_INTEGER_PROGRAMMING')
        y = {}
        for r in range(NumberofVariable):
            y[r] = master.IntVar(lb=0.0, ub=float('inf'), name='y_%s' % r)
        custconstr = {}
        for i in range(1, self.numnodes):
            custconstr[i] = master.Add(
                sum(col_coeff[r][i] * y[r] for r in range(NumberofVariable)) - 1 >= 0, 'cust_%s' % i)
        master.Minimize(sum(route_length[r] * y[r] for r in range(NumberofVariable)))
        master.Solve()
        int_solution = master.Objective().Value()  # 当前节点的整数解，即当前节点的目标函数值
        self.objective_value = int_solution
        # print("整数解:    ", int_solution)
        path_index = []
        for i in y:
            if (y[i].solution_value() != 0):
                path_index.append(i)
        # node_routes_int = copy.deepcopy([routes_pool[i] for i in path_index])  # 当前状态的整数解路径，该项不需要传递，只需要求解结束时返回
        # print("该节点求解结束时的整数路径：", node_routes_int)
        return self._bound


    def save_state(self, node):
        """
        solved_continuous_routes_pool: 求解后的路径池
        self.solved_ban_edges：子节点中要去除的边
        self.solved_set_edges：子节点中要保留的边
        self.edge_handled：已经处理过的边，避免重复处理
        """
        node.state = (self.solved_continuous_routes_pool, self.solved_ban_edges, self.solved_set_edges, self.edge_handled)


    def load_state(self, node):
        """
        ini_continuous_routes_pool: 求解前的路径池
        ini_ban_edges： 求解前要去除的边
        ini_set_edges： 求解前要保留的边
        edge_handled： 已经处理过的边
        """
        (self.ini_continuous_routes_pool, self.ini_ban_edges, self.ini_set_edges, self.edge_handled) = node.state

    def get_edge_from_route(self, route):
        """
        返回一条路径中包含的边
        """
        result = []
        for i in range(len(route) - 1):
            result.append((route[i], route[i + 1]))
        return result

    def choose_edge_to_branch(self, routes, coeff):
        """
        返回作为分支的边
        """
        max_temp = 100
        index_chosen = 100
        max_loop = 10
        loop_num = 0
        for i in range(len(routes)):
            if max_temp > (coeff[i] - 0.6) ** 2:
                max_temp = (coeff[i] - 0.6) ** 2
                index_chosen = i
        edge2chosen = self.get_edge_from_route([0] + routes[index_chosen] + [0])
        try:
            index_chosen = random.randint(0, len(edge2chosen) - 1)
        except:
            print("错误，此时edge2chosen为", edge2chosen)
            print("错误，此时index_chosen为", index_chosen)
            print("错误，此时routes为", routes)
        result = edge2chosen[index_chosen]
        while (result in self.edge_handled):
            loop_num += 1
            if (loop_num == max_loop):
                loop_num = 0
                for i in range(len(routes)):
                    if max_temp > (coeff[i] - random.random()) ** 2:
                        max_temp = (coeff[i] - random.random()) ** 2
                        index_chosen = i
                edge2chosen = self.get_edge_from_route([0] + routes[index_chosen] + [0])
            index_chosen = random.randint(0, len(edge2chosen) - 1)
            result = edge2chosen[index_chosen]
        return result

    def branch(self):
        """
        分支过程，生成子节点
        """
        child_list = []
        # 找到要分支的边
        edge_to_branch = self.choose_edge_to_branch(self.solved_continuous_routes_pool, self.solved_continuous_coeff)
        self.edge_handled.append(edge_to_branch)
        # 左孩子去除这条边
        child = pybnb.Node()
        new_ban_edges = self.ini_ban_edges + [edge_to_branch]
        new_set_edges = self.ini_set_edges
        child.state = (self.solved_continuous_routes_pool, new_ban_edges, new_set_edges, self.edge_handled)
        child_list.append(child)
        # 右孩子保留这条边
        child = pybnb.Node()
        new_ban_edges = self.ini_ban_edges
        new_set_edges = self.ini_set_edges + [edge_to_branch]
        child.state = (self.solved_continuous_routes_pool, new_ban_edges, new_set_edges, self.edge_handled)
        child_list.append(child)
        return child_list


if __name__ == '__main__':
    import pybnb.misc
    file_path = 'G:/VRP_BBP-master/Rc101Data/40-0.csv'
    file_path_cpp = "G:/VRP_BBP-master/Rc101Data/40-0.txt"

    problem = VRPTW(file_path, file_path_cpp)
    pybnb.misc.create_command_line_solver(problem)