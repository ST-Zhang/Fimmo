if size(f{1},3)==1
    f{1}(:,:,2)=0;
    I00(1)=2;
    I10(1)=2;
    I20(1)=2;
end
for Isp=1:2
    n_{Isp}=trapz(v{Isp},4*pi*f{Isp}(:,:,1).*v{Isp}.^2,2); %#ok<SAGROW> 
    u_{Isp}=trapz(v{Isp},4*pi/3*f{Isp}(:,:,2).*v{Isp}.^3,2)./n_{Isp}; %#ok<SAGROW> 
    T_{Isp}=(4*pi/3)*mass(Isp)*trapz(v{Isp},f{Isp}(:,:,1).*v{Isp}.^4,2)./n_{Isp}; %#ok<SAGROW> 
        %-(16*pi^2/27)*trapz(v{Isp},f{Isp}(:,:,2).*v{Isp}.^3,2).^2./n_{Isp}.^2; %#ok<SAGROW> 
    q{Isp}=1/2*mass(Isp)*(...
        -u_{Isp}.*(16*pi/15).*trapz(v{Isp},f{Isp}(:,:,I20(Isp)).*v{Isp}.^4,2)...
        -u_{Isp}.*(8*pi/3).*trapz(v{Isp},f{Isp}(:,:,I00(Isp)).*v{Isp}.^4,2)...
        +u_{Isp}.^2*(4*pi).*trapz(v{Isp},f{Isp}(:,:,I10(Isp)).*v{Isp}.^3,2)...
        +(4*pi/3).*trapz(v{Isp},f{Isp}(:,:,I10(Isp)).*v{Isp}.^5,2)...
        -u_{Isp}.*(2*pi).*trapz(v{Isp},f{Isp}(:,:,I00(Isp)).*v{Isp}.^4,2)...
        ); %#ok<SAGROW> 
end

dT_dx{1}=zeros(dimx,1);
dT_dx{2}=zeros(dimx,1);
for Ix=1:dimx
    if Ix==1
        dT_dx{1}(1)=(T_{1}(2)-T_{1}(dimx))/2/dx;
        dT_dx{2}(1)=(T_{2}(2)-T_{2}(dimx))/2/dx;
    elseif Ix==dimx
        dT_dx{1}(dimx)=(T_{1}(1)-T_{1}(dimx-1))/2/dx;
        dT_dx{1}(dimx)=(T_{2}(1)-T_{2}(dimx-1))/2/dx;
    else
        dT_dx{1}(Ix-1)=(T_{1}(Ix+1)-T_{1}(Ix-1))/2/dx;
        dT_dx{2}(Ix-1)=(T_{2}(Ix+1)-T_{2}(Ix-1))/2/dx;
    end
end
kappa_simu=-q{1}./dT_dx{1};
kappa_simu_ion=-q{2}./dT_dx{2};
kappa_simu(q{1}==0|dT_dx{1}==0)=NaN;
kappa_simu_ion(q{2}==0|dT_dx{2}==0)=NaN;

n_e_cm3=n_{1}*ne;
n_i_cm3=n_{2}*ne;
omega_pe_sec=5.64e4*sqrt(n_e_cm3);
T_e_eV=T_{1}*5.11e5;
T_i_eV=T_{2}*5.11e5;
coeff1=[3.16 4.9  6.1  6.9];
if Z(2)==1
coeff2=3.906;
end
lnL_ei=(T_e_eV>=10*Z(2)^2).*(24-log(n_e_cm3.^(1/2)./T_e_eV))+...
            (T_e_eV<10*Z(2)^2).*(23-log(n_e_cm3.^(1/2)*Z(2).*T_e_eV.^(-3/2)));
tau_ei_sec=3.44e5*T_e_eV.^(3/2)./(n_{1}*ne)./lnL_ei/Z(2);
tau_ei=omega_pe_sec.*tau_ei_sec;
lnL_ii=23-log(Z(2)^2./T_i_eV.*sqrt(2*n_i_cm3*Z(2)^2./T_i_eV));
tau_i_sec=2.09e7*T_i_eV.^(3/2)./(n_{2}*ne)./lnL_ii*sqrt(mass(2)/1836)/Z(2)^2;
tau_i=omega_pe_sec.*tau_i_sec;
kappa_theo_e=coeff1(Z(2))*n_{1}.*T_{1}.*tau_ei;
if Z(2)==1
kappa_theo_i=coeff2(Z(2))*n_{2}.*T_{2}.*tau_i/mass(2);
end