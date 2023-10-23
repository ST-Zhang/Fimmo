# Fimmo
This is a code to simulate collisions and laser-plasma interactions of multiple-species plasmas by solving 1D3V Vlasov-Fokker-Planck-Maxwell equations.
主程序是FP_main_vxxx.m，计算前先保证当前文件夹是这个程序的所在文件夹，再运行这个程序。

通过一个输入文件来设定初始条件，该文件命名为INPUT.m。

关于如何写INPUT.m，会在下面给出，在Eaxmples文件夹中有几个例子，这些例子为Vlasov-Fokker-Planck-Maxwell simulations for plasmas in inertial confinement fusion.pdf中的算例。

在Vlasov-Fokker-Planck-Maxwell simulations for plasmas in inertial confinement fusion.pdf中还给出了程序的算法，以及程序的归一化条件。

主程序的第6、7行用于给出INPUT.m的所在目录，这也是输出文件的所在目录。

输出文件包含fxxxxxx.mat和Fieldxxxxxx.mat两种，前者包含分布函数，后者包含各场。关于输出变量和实际分布函数和电磁场的关系，参见Vlasov-Fokker-Planck-Maxwell simulations for plasmas in inertial confinement fusion.pdf。



一、INPUT.m的写法

1.1 % constants 这部分可以定义变量，注意不要与其他变量命名重复。

1.2 % species 这部分定义等离子体

dimsp 粒子种类，至少为2，第一个必须为电子

Z 粒子电荷量

mass 粒子质量与电子质量之比

number_ratio 粒子数之比，第一个元素必须为1

1.3 % control 这部分定义网格

Nx 位形空间网格数

xmax 位形空间长度

Nv 各粒子的速度空间网格数

vmax 各粒子的速度空间最大值

Nt 时间步数

dt 时间步长

lmax 谐波指标l的最大值，为非负整数，对于各向同性等离子体的模拟可设为0，碰撞越弱，需要越多的lmax来保证精度

mmax 谐波指标m的最大值，0≤m≤l，对于轴对称系统（如无横向场且无外加磁场时）可设为0

1.4 % characteristic 这部分给出系统特征以简化计算

flag.notransfield 没有横向场时设为1（如无激光场），否则为0

flag.isuniform 没有对流和纵向场时设为1（如等离子体状态与x无关），否则为0

1.5 % output 这部分设定输出文件

flag.read_tensor 程序首次运行时会在INPUT.m所在文件夹中生成一个tensor.mat文件，当flag.read_tensor=1时，将会监测INPUT.m所在文件夹中是否有该文件，如果有则直接使用，如果无则进行生成；当flag.read_tensor=0时，将不检查直接生成tensor.mat，如果已有则覆盖之。

flag.loop_pause 进行第一个步长后是否自动暂停

flag.write_pause 输出第一个步长后是否自动暂停

flag.init_pause 初始化（包括计算初始分布函数等）后是否暂停

flag.debug_info 在命令行窗口显示一些运行的信息

file_write 写文件间隔，一个dt为一个间隔

file_start_number -1：从分布函数算起；非负：从编号为file_start_number的文件算起

1.6 % interruption 这部分给出程序监测异常中断的条件

number_tol 当电子总数变化率超过此值时结束运行,略大于零，等于0时无效（一般设为0）

max_multiple 当电子分布函数的最大值超过原先的此值的倍数时结束运行，等于0时无效（建议值15或更大）

1.7 % distributions 定义初始分布

dist_ini 初始分布，给定各粒子分布函数的名称和参数，对于一个名称为“xxxxx”的函数，需要额外写一个Distribution_xxxxx.m文件，关于如何写此文件，会在下面给出。主程序所在文件夹中提供了几个例子，也是Example文件夹中用到的初始分布

density_profile_ini 给初始分布附加一个额外的密度轮廓，对于密度轮廓已经在Distribution_xxxxx.m文件写好的情况，可以设定为'default'；其他情况参见INPUT.m的注释。

dist_neu_ini 对某一粒子的初始密度进行调整使得等离子体无净电荷。设定为'none'：不调整；设定为{'overall',Num}对第Num种粒子的分布函数乘以一个常数以实现整体电中性；设定为{'everywhere',Num}对第Num种粒子的分布函数乘以与x有关的常数以实现处处电中性。

1.8 % collisions 碰撞的开关

flag.collision_on 碰撞总开关，1：有碰撞；0：无碰撞

flag.collision_on_ee 电子-电子碰撞开关

flag.collision_on_ie 离子-电子碰撞开关，注意该碰撞无法保证能量守恒，如无必要建议关掉

flag.collision_on_ei 电子-离子碰撞开关

flag.collision_on_ii 离子-离子碰撞开关

1.9 % conservation 这部分决定是否进行守恒性修正，一般不需要修改

1.10 % laser 设定激光参数

flag.laser_on 激光开关

laser_L 左边界入射激光。行：y偏振，z偏振；列：振幅，频率，相位，持续时间

laser_R 右边界入射激光。设定同上。

1.11 % External fields 设定外加场

Bx0 外加纵向静磁场

flag.extenal_fields 是否有外加电场

Ey0_extenal 当flag.extenal_fields=1是启用，用匿名函数定义外加电场

1.12 % boundaries 定义边界

bdy.f 分布函数的边界，'periodic'周期性，'vaccum'真空，'even'边界导数为0，'invariable'保持初始值不变。注意含有reflect的边界目前没进行调试

bdy.Ft 横向场的边界，'periodic'周期性，'outflow'向外传出，'laser'激光

1.13 % IB heating 不要使用

1.14 % Phenomenological heating 通过将一定比例的分布函数替换为给定温度的Maxwellian分布，从而实现唯象加热（并保持温度）的效果

flag.heating 开关

alpha 温度变化率=alpha*(Te_hot-T(t))/dt

Te_hot_eV 用匿名函数规定加热的温度轮廓

1.15 % probes 场的探针

flag.probe_on 开关，为1时在INPUT.m所在文件夹中输出一个probe.mat的文件，给出给定位置probe_x的横向场在各时刻的值

probe_x 插入探针的位置，探针在最接近该位置的网格

1.16 % algorithm 设置算法

alg.timepush 时间推进算法，可选参数：'rk4'(4阶Runge-Kutta),'Heun3'(3阶Heun),'midpoint'(2阶midpoint),'Heun2'(2阶Heun)

alg.convection 对流项差分格式，可选参数：'central'(中心差分),'TVD'(二阶TVD差分)

alg.xgridrbdy 横坐标最大值包不包含xmax，可选参数：'exclude'(横坐标不包含x=xmax点);'include'(横坐标包含x=xmax点);'default'(边界周期时不包含，否则包含)

alg.IBheating 没有用



二、初始分布Distribution_xxxxx.m的写法

2.1 第1行：f=Distribution_xxxxx(x,v,Isp,VAR1,VAR2,...)中，xxxxx表示函数名，VAR1,VAR2,...为需要用到的参数，Isp为粒子种类。

2.2 第2~3行不用修改

2.3 第4行一般为[V,X]=meshgrid(v,x); 分布与X无关时也可以写为[V,~]=meshgrid(v,x);

2.4 从第5行开始写分布函数，分布函数命名为f_ini，位形空间命名为X，速度空间命名为V。注意使用'.*'，'./'，'.^'等运算符。程序会按照number_ratio变量自动归一化，因此可以不归一化

2.5 最后一行不用修改
