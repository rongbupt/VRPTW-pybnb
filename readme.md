**注意：**

请在python3.9环境下运行，如果是其他python版本，请重新编译_espprc.pyd

**安装需要的包：**

```python
pip install -r requirements.txt
```

Windows平台请先下载安装msmpi

**运行命令：**

```
mpiexec -n 10 python VRPTW-pybnb.py
```

**切换算例：**

请在VRPTW-pybnb.py文件中修改算例文件地址，如下图：

![](G:\行业文章\VRPTW-pybnb-code\address.jpg)

因为一些历史原因，python读取的文件和C++读取的文件格式不同，但其实内容是一样的。切换算例的时候需要两者同时修改。