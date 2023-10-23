load('Info.mat');
file_space=20000:500:40000;
x_fft_interval=30;%越大dk越小但越不局域
k0_range=[10,15];%重要！记得改
envelope_Index_length=20;

Ex_x_t=zeros(dimx,length(file_space));
for Ifile=1:length(file_space)
    load(['Field',num2str(file_space(Ifile),'%.6d'),'.mat']);
    Ex_=Ex-mean(Ex);
    Ex_x_t(:,Ifile)=Ex_;
end

x_fft_interval_Index=round(x_fft_interval/dx);
gamma_x_k0=zeros(1,dimx-x_fft_interval_Index);
gamma_x_k0_S=zeros(1,dimx-x_fft_interval_Index);
k0_x=zeros(1,dimx-x_fft_interval_Index);
k_space=2*pi/dx*(0:x_fft_interval_Index/2-1)/x_fft_interval_Index;
[~,Ik0_left]=min(abs(k0_range(1)-k_space));
[~,Ik0_right]=min(abs(k0_range(end)-k_space));
for Ix=1:dimx-x_fft_interval_Index+1
    Ex_k_t=zeros(x_fft_interval_Index,length(file_space));
    for It=1:length(file_space)
        Ex_k_t(:,It)=abs(fft(Ex_x_t(Ix:Ix+x_fft_interval_Index-1,It)))/x_fft_interval_Index;
    end
    Ex_k_t=Ex_k_t(1:floor(end/2),:);
    time=file_space*dt;
    k_gamma=zeros(x_fft_interval_Index,1);
    S_gamma=zeros(x_fft_interval_Index,1);
    for Ik=Ik0_left:Ik0_right
        Ex_k0_t_amp=envelope(Ex_k_t(Ik,:),envelope_Index_length,'peak');
        [p,S]=polyfit(time,log(abs(Ex_k0_t_amp)),1);
        k_gamma(Ik)=p(1);
        S_gamma(Ik)=S.normr;
    end
    
    gamma_x_k0(Ix)=max(k_gamma(Ik0_left:Ik0_right));
    [~,Ik0]=min(abs(k_gamma-gamma_x_k0(Ix)));
    k0_x(Ix)=k_space(Ik0);
    gamma_x_k0_S(Ix)=S_gamma(Ik0);

end
figure('Color',[1 1 1]);
plot(x(1:dimx-x_fft_interval_Index+1)+envelope_Index_length/2,gamma_x_k0);
xlabel('x');ylabel('\gamma');
figure('Color',[1 1 1]);
plot(envelope_Index_length/2:dx:envelope_Index_length/2+(length(gamma_x_k0)-1)*dx,k0_x,'.','MarkerSize',10);
xlabel('x');ylabel('k_{max}');
figure('Color',[1 1 1]);
plot(envelope_Index_length/2:dx:envelope_Index_length/2+(length(gamma_x_k0)-1)*dx,gamma_x_k0_S);
xlabel('x');ylabel('S.normr');