load('Info.mat')
file_space=(round(440/dt/100)*100-round(2*pi*200):1:round(440/dt/100)*100);
omega=laser.freq_L_y;
load('probe.mat');
%Z(2)=48;
T_e_eV=vTe^2*511000;

time=(1:length(probe(1).Ex))*dt;
S1=(probe(1).Ey.*probe(1).Bz)-(probe(end).Ey.*probe(end).Bz);%算线性密度
%S1=(probe(2).Ey.*probe(1).Bz)-(probe(4).Ey.*probe(end).Bz);%算均匀密度
%S_bar=(max(S1(file_space))+min(S1(file_space)))/2;  %传入-传出
S_bar=mean(S1(file_space));
S0=laser.amp_L_y^2/2;       %传入
disp(S_bar/S0)
figure('Color',[1 1 1]);plot(time,S1,time(file_space),S1(file_space));xlabel('t');ylabel('\alpha_{abs}')

if flag.collision_on
    n_e_cr=omega^2;
    n_e_cm3=n_e_cr*ne;
    lnL_ei=(T_e_eV>=10*Z(2)^2).*(24-log(n_e_cm3.^(1/2)./T_e_eV))+...
            (T_e_eV<10*Z(2)^2).*(23-log(n_e_cm3.^(1/2)*Z(2).*T_e_eV.^(-3/2)));
    omega_pe_sec=5.64e4*sqrt(n_e_cm3);
    tau_ei_sec=3.44e5*T_e_eV.^(3/2)./n_e_cm3./lnL_ei/Z(2);
    tau_ei=tau_ei_sec.*omega_pe_sec;
    coeff_28=[0.51,0.44,0.40,0.38,NaN,1];%0.29
    %nu_ei_star=coeff_28(end)*omega^2/tau_ei;
    %disp(1-exp(-32*nu_ei_star*L/15))

    nu_ei_star=omega^2/tau_ei;
    disp(1-exp(-32*nu_ei_star*L/15))
end