addpath 'E:\FokkerTransport\Code\amend04_benchmark\06.HeatConduct'
load('Info.mat')
file_space=0:400:Nt;
kappa_simu_t=zeros(1,length(file_space));
kappa_simu_ion_t=zeros(1,length(file_space));
kappa_theo_t=zeros(1,length(file_space));
kappa_theo_ion_t=zeros(1,length(file_space));
qe=zeros(1,length(file_space));
qi=zeros(1,length(file_space));
%sigma_simu_t=zeros(1,length(file_space));
time=zeros(1,length(file_space));
for Ifile=1:length(file_space)
    load(['f',num2str(file_space(Ifile),'%.6d'),'.mat']);
    load(['Field',num2str(file_space(Ifile),'%.6d'),'.mat']);
    time(Ifile)=file_space(Ifile)*dt;
    kappa_sim;
    %kappa_simu(isinf(kappa_simu))=nan;
    %sigma_simu=-u_{1}.*n_{1}./(Ex+Ex0_);
    kappa_simu_t(Ifile)=kappa_simu(dimx/2+1);%mean(kappa_simu,'omitnan');
    kappa_simu_ion_t(Ifile)=kappa_simu_ion(dimx/2+1);
    kappa_theo_t(Ifile)=kappa_theo_e(dimx/2+1);
    if Z(2)==1
        kappa_theo_ion_t(Ifile)=kappa_theo_i(dimx/2+1);
    end
    qe(Ifile)=q{1}(dimx/2+1);
    qi(Ifile)=q{2}(dimx/2+1);
    %sigma_simu_t(Ifile)=sigma_simu(dimx/2+1);%mean(kappa_simu,'omitnan');
end

for Isp=1:2
    n_{Isp}=trapz(v{Isp},4*pi*f{Isp}(:,:,1).*v{Isp}.^2,2); %#ok<SAGROW> 
    u_{Isp}=trapz(v{Isp},4*pi/3*f{Isp}(:,:,2).*v{Isp}.^3,2)./n_{Isp}; %#ok<SAGROW> 
    T_{Isp}=(4*pi/3)*mass(Isp)*trapz(v{Isp},f{Isp}(:,:,1).*v{Isp}.^4,2)./n_{Isp}; %#ok<SAGROW> 
        %-(16*pi^2/27)*trapz(v{Isp},f{Isp}(:,:,2).*v{Isp}.^3,2).^2./n_{Isp}.^2; %#ok<SAGROW> 
end

n_e_cm3=n_{1}*ne;
n_i_cm3=n_{2}*ne;
omega_pe_sec=5.64e4*sqrt(n_e_cm3);
T_e_eV=T_{1}*5.11e5;
T_i_eV=T_{2}*5.11e5;

coeff1=[3.16 4.9  6.1  6.9];
if Z(2)==1
coeff2=3.906;
end

% %NRL
% lnL_ei=(T_e_eV>=10*Z(2)^2).*(24-log(n_e_cm3.^(1/2)./T_e_eV))+...
%             (T_e_eV<10*Z(2)^2).*(23-log(n_e_cm3.^(1/2)*Z(2).*T_e_eV.^(-3/2)));
% tau_ei_sec=3.44e5*T_e_eV.^(3/2)./(n_{1}*ne)./lnL_ei/Z(2);
% tau_ei=omega_pe_sec.*tau_ei_sec;
% lnL_ii=23-log(Z(2)^2./T_i_eV.*sqrt(2*n_i_cm3*Z(2)^2./T_i_eV));
% tau_i_sec=2.09e7*T_i_eV.^(3/2)./(n_{2}*ne)./lnL_ii*sqrt(mass(2)/1836)/Z(2)^2;
% tau_i=omega_pe_sec.*tau_i_sec;
% %Braginskii
% % lnL=(T_e_eV<50).*(23.4-1.15*log(n_e_cm3)/log(10)+3.45*log(T_e_eV)/log(10))+...
% %     (T_e_eV>50).*(25.3-1.15*log(n_e_cm3)/log(10)+2.3*log(T_e_eV)/log(10));
% % tau_ei_sec=3.5e5./lnL.*T_e_eV.^(3/2)/Z(2)./n_e_cm3;
% % tau_ei=omega_pe_sec.*tau_ei_sec;
% % tau_i_sec=3.0e7./lnL.*sqrt(mass(2)/2/1836).*T_i_eV.^(3/2)/Z(2)^3./n_i_cm3;
% % tau_i=omega_pe_sec.*tau_i_sec;
% 
% kappa_theo_e=coeff1(Z(2))*n_{1}.*T_{1}.*tau_ei;
% if Z(2)==1
% kappa_theo_i=coeff2(Z(2))*n_{2}.*T_{2}.*tau_i/mass(2);
% end

figure('color',[1 1 1]);
plot(time,kappa_simu_t,'b','LineWidth',2.0);xlabel('t');
hold on;plot(time,kappa_theo_t,'r--','LineWidth',2.0);ylabel('\kappa_e');
set(gca,'ticklength',[0.03,0.06],'xminortick','on','yminortick','on','linewidth',2.0,'fontsize',20,'fontweight','bold');
if Z(2)==1
figure('color',[1 1 1]);
plot(time,kappa_simu_ion_t,'b','LineWidth',2.0);xlabel('t');
hold on;plot(time,kappa_theo_ion_t,'r--','LineWidth',2.0);ylabel('\kappa_i');
set(gca,'ticklength',[0.03,0.06],'xminortick','on','yminortick','on','linewidth',2.0,'fontsize',20,'fontweight','bold');
end