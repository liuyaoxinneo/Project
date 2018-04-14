# 非线性相控阵成像程序：neo_1
---
## once：4个一次积分
问题：
- 四个积分对应的含义？
- once_jifen所用的积分的数值计算方法是什么？？
- 4：h没有定义类型
- 22：t1没有定义类型
- 40：是否应该是`t1=t2`？？否则28行计算t2之后的值没有传递下去，而且jifen1的返回值为t2,而t2=t1,t1一直没有发生过变化
- 变量ep的含义？？猜想是不是在测量两次递归之间的误差
```
if (t1==0)then
    ep=abs(t2-t1)
elseif(t1==t2)then
    ep=0
else
    ep=abs(1-t2/t1)
endif
```
-

理解：
- 猜想数值积分方法用的是递归梯形公式  
 *——参考《数值方法》P292*
 \[t1=T(0)\\t2=T(1),T(2)....,T(9)\]
 \[once\_jifen1=\int\nolimits_{y_{down}}^x density() dy \]
 \[once\_jifen2=\int\nolimits_x^{y_{up}} density() dy \]
 \[once\_jifen3=\int\nolimits_{x_{donw}}^y density() dx \]
 \[once\_jifen4=\int\nolimits_y^{x_{up}} density() dx \]


程序中的字母|含义|计算位置
--|--|--
k,h  |   |  
t1,t2  | $t1=T(0)$，t2：存储各次迭代的结果，$t2=T(J)$  |  
hh  | 存储各次迭代区间长度，$h_0=hh=a-b,h_1=\frac{h_0}{2},···,h_9=\frac{h_0}{2^9}$  |  
m,n  | m：将积分区间[a,b]分为2m个区间，即，m为递推的次数；n：$2m=2^n$  |  
x  |   |  
x(2),y(2)  | x(2):对x积分时的上下界，y(2):对y积分时的上下界  |  
xx,yy  | 存储各子区间的x和y坐标  |  
once_result  |   |  
y_down,y_up  |   |  
x_down,x_up  |   |  
**sign**  |   |  
ep  |   |  
eps  |   |  
---
## twice：对once中的四个一次积分 按照定义再各积一次
设$density()=f()$
\[twice\_jifen1(xstart,xend,y\_down,sign)=\int\nolimits_{xstart}^{xend}\int\nolimits_{y-down}^x f() dx dy\]
\[twice\_jifen2=\int\nolimits_{xstart}^{xend}\int\nolimits_x^{y-up} f() dx dy\]
\[twice\_jifen3=\int\nolimits_{ystart}^{yend}\int\nolimits_{x-down}^y f() dx dy\]
\[twice\_jifen4=\int\nolimits_{ystart}^{yend}\int\nolimits_y^{x-up} f() dx dy\]
程序中的字母|含义|计算位置
--|--|--
xstart,xend  |   |  
ystart,yend  |   |  
s1,s2  |   |  

---
## pm_density
1. 声明变量：$A=1.5\times10^{8},P_{up}=3\times10^{8},P_{down}=1\times10^{7}$
2. 定义了函数：density(Po,Pc,sign)

问题：
- 定义的A,P_up,P_down在本模块的函数`density`中并没有使用？？函数`density`的返回值只和sign有关，那么设置其输入变量为三个的意义何在？

---
## pm_reaction_rewrite
定义了子程序`pm_caculate`

程序中的字母|含义|计算位置
--|--|--
Khr_line  |应力-应变图的斜率   |  
Kh1,Kh2  |开时的模量，关时的模量   |  
riseflag(m,n)  |T：该单元的应力增长   |  
---
## pmsimu
问题：
- 69：`do h=68,68`的意义？实际效果为只循环一次
- 87~88：根据参考文献44 14(b) 应该为$K_s=2\mu$
- 170：`k=1,40000`的含义？是否应该为`k=1,200000`
- 193：k肯定小于5e5,判断条件多余
- 197、217：成立条件：Txx时间导数为0
- 199：flag的作用
- 200：`t=k*deltt/2`中为什么要除以2？
- 200、201：Focus=$\sqrt{(m-xxx)^2+yyy^2}\times\Delta d$，其中(xxx,yyy)为(m,n)坐标系下的坐标
- 203：计算过程太复杂，可以简略为$t_n=\frac{Focus}{speed}\times\sqrt{1-(degree)^2}$
- 214：比193行少一项-Txx(m,n)
- 220：判断条件少了一个0，应该为`k>500000`！！！
- 230：最左侧不正对源的部分没有计算vx,Txx，意义何在？
- 234~240：初判断条件外，与上几行完全一致，**可以省略**
- 218、252：成立条件：Txx时间导数为0
- 222：成立条件：U为常数
- 260：n取值范围奇怪，是否是因为上方还有一个源？
- 260：与214相似，少一项-Tyy(m,n),是不是默认为Tyy(m,n)=0？
- 264、273：成立条件：Tyy(m,n)=0
- 265、277：成立条件：$\dot{T}_{yy}(m,n)=0$
- flag：是否是用于切换V和T的计算的？
- 272：似乎缺少一行vyf(m-1,n)的计算
- 289：成立条件：Txy(m-1,n)=0，即Tubxy(n)=0
- 287、288：对速度乘衰减因子？
- 296、284：成立条件：vyf(m-1,n)=0
- 298、355：整个区域内各单元的U均相同？
- 298：经典非线性
- 312：判断条件有错，只有前一个有效
- 320：Txxb,Tyyb,Txyb并没有传值和初始化？？
- 320、321：存在的意义？？
- 327~329：使用了张略ref_44 13(b)，$\dot{T}_{v}=K_v\times\dot{\varepsilon}_{v}$ 离散化为：$\dot{T}_{v}=\frac{Tvf-Tv}{\Delta t}$，$\dot{\varepsilon}_{v}=\frac{\varepsilon_{xx}+\varepsilon_{yy}}{\sqrt{2}}$，$\varepsilon_{xx}=\frac{vxf(m,n)-vxf(m,n-1)}{\Delta d}$，$\varepsilon_{yy}=\frac{vyf(m,n)-vyf(m-1,n)}{\Delta d}$， $\varepsilon_{xy}=\frac{vxf(m+1,n)-vxf(m,n)+vyf(m,n+1)-vyf(m,n)}{\Delta d}$
- 341：调用了pm_caculate()
- 361：用K11代替了标准迭代公式中的K1
- 372：如果是最后一列，则判断为corner？？判断条件有误
- 382~383、388~389、394~395、399~400：设速度差是线性递减的,成立条件：$\Delta V_1+\Delta V_3=2\times\Delta V_2$
- 387：成立条件：Txx(m,n)=0
- 405、406：成立条件：vx(1,0)=0,vy(1,0)=0，即，vlbx(1)=0
- 411、417：成立条件：vx(301,1)=0，vy(1,501)=0
- 445：Trr()存储了(100,500)单元的vx的值，意义何在？
- 470：`if(1)`的含义
- 475：2E4或许可能为2E5？？
- 493、494：需要close文件

需要统计每个参数的含义！整理成表格  
变量统计表

程序中的字母|含义|计算位置
--|--|--
K1,K2,U  |  弹性模量：$K_1$,$K_2$,$\mu$ |  
p  | 试块密度(程序中均取为**铁**) $\rho=7700 kg/m^3$  |  
w  | $\omega=2{\pi}f$  |  
Q  |与衰减有关的品质因子，应用于 $\alpha=e^{{-\pi f\Delta t}/2Q}$   |  
  B|  单个网格的面积：$B={\Delta d}^2$ |  
delt  | $delt=\Delta t/\Delta d$  |  
r1,r2  | 改进的PM模型中：$r_1$，$r_2$  |  
Kh1,Kh2  | 改进的PM模型中：$K_{M1}$，$K_{M2}$  |  
initial=1e20  | 用于给Hem,s_Hem来初始化  |  
pm,xx,yy(2)  |   |  
 deltat , deltad| 时间步长：$\Delta t$ , 网格边长：$\Delta d$   |  
Kv(301,501) Kd(301,501) ks(301,501)  |  $K_V=K_1+K_2$, $K_D=K_1-K_2$, $K_S=2\mu$, **此处程序与公式不符！！**    |  
x(m) , xf(m)  |本次迭代前的变量，迭代后的变量  
vx(300,500) , vy(300,500)| 迭代前某单元的$V_x$、$V_y$|单元的边上
vlbx(300)  |left boundary 左侧边界上的$V_x$，相当于vx(300,0)，边界条件  
vuby(500)  |upper boundary 上侧边界上的$V_y$，相当于vy(0,500)，边界条件
Txx(300,500) , Tyy(300,500)|应力，$\sigma_{xx}$，$\sigma_{yy}$|单元中心
Txy(300,500)  |应力，$\sigma_{xy}$   |  单元右下角
Tubxy(500)  |upper boundary 上侧边界上的 $\sigma_{xy}$ ，相当于Txy(0,500)  |  
Tlbxy(300)  |left boundary 左侧边界上的 $\sigma_{xy}$ ，相当于Txy(300,0)  |  
Txxb(300,500) , Tyyb(300,500) , Txyb(300,500)  |   |  
Tv,Td,Ts  | 用特征系统表示弹性张量时的参数，应用于公式：$T_v=T_{xx}+T{yy}/\sqrt{2}$，$T_d=T_{xx}-T{yy}/\sqrt{2}$，$T_s=T_{xy}$  |  
flag  |   |  
cornerflag  |   |  
riseflag  |   |  
initial_state  |   |  
s_riseflag  |   |  
s_initial_state  |   |  
Hem ,s_Hem  |   |  
A_input,f,Fx,Fy  |$F_x=A_{input}sin(2\pi f t),F_y=0$ ，A_input：输入简谐外力的幅值，f：简谐外力的频率    |  
f  | 左侧发射换能器的频率  |  
fL  | 上方发射换能器的频率(用以调制波的研究)  |  
fatten  | 声源的频率？  |  
K11(301,501),bb  | 经典非线性部分对体模量的贡献，应用于公式：$K_C(\sigma)=K_0(1+\beta\sigma+\delta\sigma^2+···)$，其中：K11即为$K_c$，bb为 $\beta$，K1为$K_0$  |  
h  | 写文件时所用的控制变量(见下表)  |  
m,n  | 网格位置循环控制变量，(m,n)  |  
t,k  | $t=\frac{k\times\Delta t}{2}$(程序中的公式)，t为时间计数，k为采样点数目 |  
xxx,yyy  | 相控阵聚焦点坐标：(190,190)  |  
Focus,speed,degree  | 焦距：，材料中的声速，偏转角度的正弦值  |  
tn  |相控阵聚焦于(xxx,yyy)时，从(m*deltad,0)发射的声束的延迟时间  |  
mm,nn  | 计算四个角落上单元时的控制变量(mm,nn)  |  
kt,sig_row  |存储数据时， kt：时间采样点的控制变量，sig_row：行数控制变量  |  
v_recsig_x(300,200000),v_recsig_y(300,200000) | recived signal,在试块右边接收到的速度信号，存储(2~299,500)单元网格的vx和vy,连续采集0.01s(200000*0.5e-7)，采集2e5个时间采样点的信号 |  

文件代码|文件名称|存放的变量|存放内容
--|--|--|--
h  | h.txt  | Trr(其实也是vxf)  | Trr(kt)=vxf(100,column)
h+10  | h+10a.txt  | vxf  | k=5E4时，vxf存储的所有内容
h+20  | h+20b.txt  | vxf  | k=1E5时，vxf存储的所有内容


理解：
- 前提：忽略外力， $F_x=0$，$F_y=0$？并不是
- 使用了假设：在相邻的四个单元内，$\mu$ 为常数
- 90~121：对二维或者三维数组A，`A=0`表示对整个数组赋值，使数组中的每个元素都为0
- 247~250：通过分别计算四个单元的一个速度分量，完成一个单元四条边上速度分量的计算— —防止单元计算重复
- $P=-\sigma$
- 196：使用了条件 $\dot{\sigma_{xx}}(m,n)=0$
- 327~329：$\dot{T_v}=K_v(m,n)\times1/\sqrt{2}\times(\frac{\partial{V_x}}{\partial{x}}+\frac{\partial{V_y}}{\partial{y}})$

---
