% INPUT
% constants

% species:
dimsp=2;
Z=[-1,1];
mass=[1,1836];
number_ratio=[1,1]; %总粒子数之比
% control:
Nx=4;
xmax=20;
Nv=[96,72];
vmax=[0.22,3.2e-3];
Nt=100000;
dt=0.2;
lmax=[0,0];
mmax=[0,0];
% characteristic
flag.notransfield=1;%没有横向场
flag.isuniform=1;%没有对流和纵向场
% output:
flag.read_tensor=1;%1：从文件中读取张量。0：重新计算张量。
flag.loop_pause=0;
flag.write_pause=0;
flag.init_pause=1;
flag.debug_info=1;
file_write=1000;%写文件间隔，一个dt为一个间隔。
file_start_number=80000;%-1：从分布函数算起；非负：从编号file_start_number的文件算起。
% interruption
number_tol=0;%当电子总数变化率超过此值时结束运行,略大于零，等于0时无效。
max_multiple=15;%当电子分布函数的最大值超过原先的此值的倍数时结束运行，等于0时无效。
% distributions:
% 分布函数：函数名，其他参数（热速度等）。
% 密度轮廓：轮廓类型，参数。
%           注意：电子的密度轮廓类型的函数必须对ne归一化。
% 密度轮廓参数：
%           trapezoid（均匀梯形╱▔╲）：[梯形的从大到小的四个横坐标]。电子的中间密度为1
%           trapezoid-nonuniform（非均匀梯形|╱|）：[两点的横坐标]，[两点的纵坐标]
%           quad（四边形）：[四个点的横坐标]，[中间两点的纵坐标]
%           const|constant（常数）:[全空间均匀的密度]。
%           smooth（平滑）:[左，右边界坐标/xmax]，[边界平滑指数（20~200越小越平滑）]
%           default（使用分布函数本身的轮廓）
% 电中性参数：{调整的方法，调整的粒子种类}
%           none：不进行电中性
%           overall：整个等离子体不带净电荷。
%           everywhere：整个等离子体处处电中性。
dist_ini={...
    {'Bump',0.024,0.08,1/2},...
    {'Bump',3.96e-4,1.32e-3,1/2},...
    };
density_profile_ini={...
    {'const',1},...
    {'const',1},...
    };
dist_neu_ini={'overall',2};
% collisions:
flag.collision_on=1;%1：有碰撞；0：无碰撞。
flag.collision_on_ee=1;
flag.collision_on_ie=0;
flag.collision_on_ei=1;
flag.collision_on_ii=1;
ne=1e23;   % 单位cm^-3，n=1时的实际密度。 
Te=1.286e3;    % eV 影响ei碰撞截断，对结果基本无影响，大致给对即可。
Ivc=[3,3];
% conservation
flag.cons_ee_NE=1; %电子-电子碰撞粒子数和能量守恒
flag.cons_ii_NE=1; %同种离子-离子碰撞粒子数和能量守恒
flag.cons_i1i2_NE=1; %不同种离子-离子碰撞粒子数和总能量守恒
flag.cons_ee_P=0;  %电子-电子碰撞动量守恒
flag.cons_E_N=1;  %电场项粒子数守恒
flag.cons_ii_P=0; %同种离子-离子碰撞粒子数和能量守恒
% laser: 行：y偏振，z偏振；列：振幅，频率，相位，持续时间
flag.laser_on=0;
laser_L=[0.0003,sqrt(5),-pi/2,inf;...
         0.000,sqrt(10),0,0];
laser_R=[0.000,sqrt(10),-pi/2,4*2*pi/sqrt(10);...
         0.000,sqrt(10),0,0];
% External fields
Bx0=0.0;
flag.extenal_fields=0;
%Ey0_extenal=@(x,t)0.027*exp(-(x-xmax/2).^2/(50*lambda_De)).*sin(10*t);
% boundaries:
bdy.f='periodic';%periodic;vaccum;vaccum|reflect;reflect|vaccum;reflect;even(对流项边界导数为零);invariable(边界等同于初始边界)
bdy.Ft='periodic';%periodic;outflow;laser;
% IB heating 不要使用
flag.IB_heating=0;
v_OSC=@(x)0.027;%*exp(-(x-xmax/2).^2/(50*lambda_De));
% Phenomenological heating
flag.heating=0;
alpha=1e-3;        %温度变化率=alpha*(Te_hot-T(t))/dt
Te_hot_eV=@(x)(T_e0+(T_e1-T_e0)*exp(-(x-xmax/2).^2/(xmax/heatlength)^2))*1000;   %eV
% probes:
flag.probe_on=0;
probe_x=[0,xmax/4,xmax/2];
% algorithm:
alg.timepush='midpoint';%rk4(4阶Runge-Kutta),Heun3(3阶Heun),midpoint(2阶midpoint),Heun2(2阶Heun)
alg.convection='TVD';%central(中心差分),TVD(二阶TVD差分)
alg.xgridrbdy='default';%exclude(横坐标不包含x=xmax点);include(横坐标包含x=xmax点);default(边界周期时不包含，否则包含)
alg.IBheating='Langdon';%Langdon(Langdon格式，不满足泊松方程);sc(满足泊松方程) %不要使用