load('Info.mat');
file_space=80000:100:156000;
k0=[10.76,12.08];%[10.2625,10.8909];%计算给定k的增长率
envelope_Index_length=20;
% 读取数据并转换到k空间
Ex_k_t=zeros(length(x),length(file_space));
for Ifile=1:length(file_space)
    load(['Field',num2str(file_space(Ifile),'%.6d'),'.mat']);
    Ex_=Ex-mean(Ex);
    Ex_k_t(:,Ifile)=abs(fft(Ex_))/length(Ex_);
end
%pcolor(2*pi/dx*(0:length(Ex)-1)/length(Ex),(1:length(file_space))*dt,Ex_k_t');
k_space=2*pi/dx*(0:length(Ex)/2-1)/length(Ex);
Ex_k_t=Ex_k_t(1:end/2,:);
time=file_space*dt;

% 计算选取k的指标范围，并在这个范围里取平均
[~,Ik0_left]=min(abs(k0(1)-k_space));
[~,Ik0_right]=min(abs(k0(end)-k_space));
Ex_k0_t=mean(Ex_k_t(Ik0_left:Ik0_right,:),1);
% 取包络，取对数，线性拟合出增长率
Ex_k0_t_amp=envelope(Ex_k0_t,envelope_Index_length,'peak');
[p,S]=polyfit(time,log(abs(Ex_k0_t_amp)),1);
disp(S);
p(1)
figure('Color',[1 1 1]);plot(time,Ex_k0_t);xlabel('t');ylabel(['Ex(k=',num2str(k0),')']);
figure('Color',[1 1 1]);plot(time,log(abs(Ex_k0_t)));xlabel('t');ylabel(['ln |Ex(k=',num2str(k0),')|']);
%hold on;plot(time,Ex_k0_t_amp,'r')