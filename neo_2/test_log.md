# neo_2/test1 试运行日志
源码直接复制自朱金林，未作任何修改
## part1
输出文件`2x_rev.txt`~`298x_rev.txt`，按原理采集的是最右侧单元在**0.01s**的时间内接收到的信号
1. 画出中间单元`150x_rev.txt`的波形，发现
    - 通过part1，输出文件均为的为2E5个点的列向量
    - 只有前1E5个点有值
2. ~~将原始数据存入`data1`文件夹~~
3. ~~在`data2`内的数据为将once.f90中迭代过程的`t2=t1`均改为`t1=t2`的迭代结果(之后也保持了此修改)~~
4. ~~`data3`：将hysteretic region去除后的迭代结果(即，没有缺陷的情况)~~
5. `data4`：去除迟滞区域和上方的源
6. `data5`:250kHz;源为左侧5~145mm;Q=10;<=5e5 x2/ >5e4 x2;**效果符合论文原图！**(发现是不计经典非线性的！！)
7. `data6`:只将源的频率改为8.5kHz
8. `data7`:频率250kHz,Q=2
9. `data8`:347加入**经典非线性**

画出`data1`,`data2`,`data3`,`data4`中第150个单元的接收信号，**发现均与理论不太符合！！**
决定重新建立工程运行！做适当修改
## part3
matlab部分，读取文件`2x_rev.txt`~`298x_rev.txt`，进行滤波操作，输出文件`2B.txt`~`298B.txt`
1. 画出中间单元`150B.txt`的波形，发现
    - 输出文件为7000个点的列向量=滤波之后截取了前6000个点+在后面接了1000个0点

# neo_2/test2 试运行日志
在Array_Imaging(已经过备注和修改)的源码基础上，对照朱的源码进行修改：(主要是去除相控阵聚焦的部分)
- `deltad=0.2e-2` ,delta=2mm
- `Q=10`
---

- `data1`：**不正确的图像！**
- `data2`：Q=10，f=250KHz，202\208行均改为50000，无迟滞区域，无调制波，不计衰减，**时域图像不太正确，但功率谱可以看出1、2、3、4次谐波**
- `data3`：Q=2，f=8.5KHz，density=5e-13(常数)，202\208行均改为50000，无迟滞区域，无调制波，不计衰减
- `data4`：Q=2，f=8.5KHz，density=5e-13(常数)，202\208行均改为50000，无迟滞区域，无调制波，考虑衰减 **时域上接近论文原图，但是功率谱不对！**
- `data5`：Q=2，f=8.5KHz，density=5e-13(常数)，202\208行均改为500000，无迟滞区域，无调制波，考虑衰减，声源位置改为(100~200,1)(216\225行做相应修改) **时域上接近论文原图，但是功率谱不对！**
- `data6`：Q=10，f=8.5KHz，density=5e-13(常数)，202\208行均改为500000，无迟滞区域，无调制波，考虑衰减，声源位置改为(100~200,1)(216\225行做相应修改)**效果还不如data5**
- `data7`：184、190、202、208：改为5e4(意义是，输入脉冲持续了0.0025s)，deltd=0.5mm？？ **错误的图像**
- `data8`：将184、190改回5E5

参数设置|p1|part1|是否相同
--|--|--|--
$\Delta t,\Delta d$  | 5E-8,5E-4  |  5E-8,5E-4|$\checkmark$
Q  | 2  | 2  |  $\checkmark$
左侧源频率 f | 8.5KHz  | 250KHz  |  
左侧声源位置(m,n)  |  [100~200,1] | [10~290,1]  |  
正对声源的判断条件  |<=5e5,>5e4   | <=5e5,>5e4  |  $\checkmark$
中间区域是否考虑经典非线性  |  是 |  是 |  $\checkmark$
考虑衰减  | 是  | 是  |  $\checkmark$
once中  | t1=t2  | t2=t1  |  

- `data9`:once中改回为`t2=t1`

# neo_3/test1 试运行日志
直接复制neo_2/test2的源码  
在问题已回复的基础上，进行修改
1. density(po,pc,sign)，在实际程序中是一个常数，并非是分段函数，sign没有使用，A、P_up、P_down也没有使用
2. **t1=t2**
3. ep/eps在once中用作循环结束条件
4. `h=68,68`可以去除，包括68/68a/68b.txt均为调试所用
5. 5e5\5e4为声源持续时长判断语句，可以自己定
6. flag的作用：交替计算应力和速度
7. 边界和角落计算采用特殊方式，简化处理，边界不算衰减
8. 显示(180,180)单元的目的在于调试，可以去除
9. 经自己整理，发现`flag=.not.flag`的地方不太对？？？：k=1,计算vxf,vyf --> k=2,计算Txxf,Tyyf,Txyf，并迭代vx、vy和Txxf,Tyyf,Txyf -->k=3,计算下一个vxf,vyf --> ...

## 运行数据/调整
### 工程：p1--用于试验调整应力/应变的计算逻辑
- p1/data1:(···计算中···) **错误！**
- p1/data2:根据以上第9点，修改应力/应变交替计算的逻辑：
  - 发现：导出的数据均为0，没有经过迭代？？？
  - **data2数据已清空！**
```
k=n
-->flag=true,remark=1：vxf,vyf
-->flag=false
-->remark=2:Txxf,Tyyf,Txyf
-->迭代：vx=vxf,vy=vyf,Txx=Txxf,Tyy=Tyyf,Txy=Txyf
-->存入这个时刻的接送信号：v_recsig_x/v_recsig_y
-->flag=true

k=n+1
···
···
```
- p1/data3:计算逻辑仿佛已经正确，外力作用时间可能导致了时域波形后半部分的混乱
- p1/data4:`k<=1e5，有外力`--混乱部分后移、减少，可见与外力时间的设置有关
- p1/data5:`k<=1e5，有外力`，**得出了理想的时域波形！**，时域信号的混乱可能来自于反射？
基于p1对neo_p1进行修改，以同时运行两个程序
- p1/data6:在data5的基础上，加入衰减，**得出了理想的时域波形！**，但频谱可能因为算法问题仍然得不到与论文中相同的结果！

### 工程：neo_p1
- neo_p1/data1:加入了迟滞区域 n:200~300,m:50~250，经典非线性，无衰减
    - 发现：>0.005s时，时域的图像不符合--设想将时间延长
    - 发现：通过同时显示k和kt，**导出v_recsig_x时只导出了前1E5项！！！**，这与之前发现k的步长为2相一致！！！，原因就是上述第9点，k每隔2个步长才传值一次
- neo_p1/data2:加入了迟滞区域 n:200~300,m:50~250，经典非线性，有衰减
- neo_p1/data3:重跑一次，**重复出了与data2相同的结果**，但**没有解决应力/应变交替的问题**，思考可能是matlab中频谱分析的方法不太正确，因此重写了`matlab\neo_test3`
- neo_p1/data4:
    - 思考：之前的模拟重复的是原论文中2.6(a)/(b)两图，都应该加上条件：无破损区域，导出数据的波形大致形状相似，但是出现混乱的现象，原因应该在于：1.加入了破损区域；2.模拟条件中，左端的源应该是一直有信号的  
    —— 针对此，进行data4的模拟计算（暂不改变应力/应变的交替计算逻辑）
    —— **数据错误！**，data4已删除
### 现已根据p1修改了neo_p1！并且把以上错误的数据：data1~data4都删除！
- data1:按照论文修改参数，有10*10mm，中心在(95,95)的缺陷
    - 经过调试，发现：加入缺陷后，调用`pm_caculate()`函数时，存在数据溢出的问题！若缺陷区域不加入非经典非线性，则可以正常计算——说明：**`pm_caculate()`中有问题！！！**
