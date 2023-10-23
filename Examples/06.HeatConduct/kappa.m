load('Info.mat')
for Isp=1:2
    n_{Isp}=trapz(v{Isp},4*pi*f{Isp}(:,:,1).*v{Isp}.^2,2); %#ok<SAGROW> 
    u_{Isp}=trapz(v{Isp},4*pi/3*f{Isp}(:,:,2).*v{Isp}.^3,2)./n_{Isp}; %#ok<SAGROW> 
    T_{Isp}=(4*pi/3)*mass(Isp)*trapz(v{Isp},f{Isp}(:,:,1).*v{Isp}.^4,2)./n_{Isp}...
        -(16*pi^2/27)*trapz(v{Isp},f{Isp}(:,:,2).*v{Isp}.^3,2).^2./n_{Isp}.^2; %#ok<SAGROW> 
    q{Isp}=1/2*mass(Isp)*(...
        -u_{Isp}.*(16*pi/15).*trapz(v{Isp},f{Isp}(:,:,I20(Isp)).*v{Isp}.^4,2)...
        -u_{Isp}.*(8*pi/3).*trapz(v{Isp},f{Isp}(:,:,I00(Isp)).*v{Isp}.^4,2)...
        +u_{Isp}.^2*(4*pi).*trapz(v{Isp},f{Isp}(:,:,I10(Isp)).*v{Isp}.^3,2)...
        +(4*pi/3).*trapz(v{Isp},f{Isp}(:,:,I10(Isp)).*v{Isp}.^5,2)...
        -u_{Isp}.*(2*pi).*trapz(v{Isp},f{Isp}(:,:,I00(Isp)).*v{Isp}.^4,2)...
        ); %#ok<SAGROW> 
end

dT_dx=zeros(dimx,1);
for Ix=1:dimx
    if Ix==1
        dT_dx(1)=(T_{1}(2)-T_{1}(1))/2/dx;
    elseif Ix==dimx
        dT_dx(dimx)=(T_{1}(dimx)-T_{1}(dimx-1))/2/dx;
    else
        dT_dx(Ix-1)=(T_{1}(Ix+1)-T_{1}(Ix-1))/2/dx;
    end
end
coeff1=[0.71  0.9  1.0  1.1  1.5];
qT=q{1}-coeff1(Z(2))*n_{1}.*T_{1}.*u_{1};
% kappa_simu=-q{1}./dT_dx;
kappa_simu=-qT./dT_dx;

n_e_cm3=n_{1}*ne;
n_i_cm3=n_{2}*ne;
omega_pe_sec=5.64e4*sqrt(n_e_cm3);
T_e_eV=T_{1}*5.11e5;
T_i_eV=T_{2}*5.11e5;
coeff2=[3.16,4.9,6.1,6.9,12.5];

lnL_ee=23.5-log(n_e_cm3.^(1/2).*T_e_eV.^(-5/4))-(1e-5+(log(T_e_eV)-2).^2/16).^(1/2);
tau_e_sec=3.44e5*T_e_eV.^(3/2)./(n_{1}*ne)./lnL_ee/Z(2);
tau_e=omega_pe_sec.*tau_e_sec;
lnL_ii=23-log(Z(2)^2./T_i_eV.*sqrt(2*n_i_cm3*Z(2)^2./T_i_eV));
tau_i_sec=2.09e7*T_i_eV.^(3/2)./(n_{2}*ne)./lnL_ii*sqrt(mass(2)/1836)/Z(2)^2;
tau_i=omega_pe_sec.*tau_i_sec;

lnL_ei=(T_e_eV>=10*Z(2)^2).*(24-log(n_e_cm3.^(1/2)./T_e_eV))+...
            (T_e_eV<10*Z(2)^2).*(23-log(n_e_cm3.^(1/2)*Z(2).*T_e_eV.^(-3/2)));
tau_ei_sec=3.44e5*T_e_eV.^(3/2)./(n_{1}*ne)./lnL_ei/Z(2);
tau_ei=omega_pe_sec.*tau_ei_sec;


% lnL_e_Brag=25.3-1.15*log(1e23)/log(10)+2.3*log(T_e_eV)/log(10);
% tau_e_Brag=3.5e5./lnL_e_Brag.*T_e_eV.^(3/2)./n_e_cm3/Z(2) .* omega_pe_sec;

%kappa_theo=coeff(Z(2))*n_{1}.*T_{1}.*tau_e;
% kappa_theo=coeff(Z(2))*n_{1}.*T_{1}.*tau_e_Brag;
kappa_theo=coeff2(Z(2))*n_{1}.*T_{1}.*tau_ei;


figure('Color',[1,1,1]);plot(x,kappa_simu,x,kappa_theo)