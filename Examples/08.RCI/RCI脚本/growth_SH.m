load('Info.mat')
load('f020000.mat')
x_range=(x>2&x<73)|(x>77&x<148);

n_e=trapz(v{1},4*pi*f{1}(:,:,1).*v{1}.^2,2);
T_e=4*pi/3*trapz(v{1},f{1}(:,:,1).*v{1}.^4,2)./n_e;%-16*pi^2/27*trapz(v{1},f{1}(:,:,2).*v{1}.^3,2).^2./n_e.^2;
vT_e=sqrt(T_e);
n_e_cm_3=n_e*ne;
T_e_eV=T_e*5.1100e05;
lnL_ei=24-log(n_e_cm_3.^(1/2).*T_e_eV.^-1);
tau_ei=1.9410e10*sqrt(ne)*T_e_eV.^(3/2)./n_e_cm_3/Z(2)./lnL_ei;
lambda_ei=vT_e.*tau_ei;
    LT_e=nan(dimx,1);
    for Ix=2:dimx-1
LT_e(Ix)=(abs((log(T_e(Ix+1))-log(T_e(Ix-1)))/2/dx))^-1;
    end
delta_T=lambda_ei./LT_e;
lambda_De=sqrt(T_e./n_e);
n_i=trapz(v{2},4*pi*f{2}(:,:,1).*v{2}.^2,2);
omega_pi=Z(2)*sqrt(n_i/mass(2));
omega_pe=sqrt(n_e);
gamma0=sqrt(pi/8)*omega_pi.^2./omega_pe;

k=0:0.01:40;
vT_i=dist_ini{2}{2}   %vT_i=dist_ini{1}{5}
T_i=vT_i^2*mass(2);
g=Z(2)*T_e/T_i;
r=g./(1+k.^2.*lambda_De.^2);
lnL_ii=23-log(Z(2)^2./T_i.*(2*n_i*Z(2)^2./T_i).^(1/2));
lambda_ii=lambda_ei*sqrt(2).*(n_e./n_i)/Z(2).*(lnL_ei./lnL_ii).*(T_i./T_e).^2;
k_i=k.*lambda_ii;
Q=(r.^(3/2).*k_i.^2+k_i.*r.^(1/2))./(r.^(3/2).*k_i.^2+3*k_i.*r.^(1/2)+10);
G=(3*r.^3+11*r.^2+12)./(r.^3+7*r);
omega=k.*vT_i.*sqrt(r+5/3+Q.*(G-5/3));
c_s=sqrt(Z(2)*T_e/mass(2));
gamma_s=sqrt(pi/8)*omega.^4./k.^3./c_s.^3.*omega_pi./omega_pe;
ksi=(Z(2)+0.5)/(Z(2)+2.12);
p_T_SH=k.*vT_e./omega*3.2.*ksi.*delta_T;
gamma_e_SH=gamma_s.*(-1+p_T_SH);

gamma_i_H=k.*vT_i.*k_i.*(r+3.02)./(r+1.67).*(0.8*r.*k_i.^2+1.49)./(r.^2.*k_i.^4+4.05*r.*k_i.^2+2.33);
gamma_i_L=sqrt(pi/8)*r.^2.*exp(-r/2-G/2).*(10+21*r+r.^3)./(2*r.^2+r.^3);
R=1./(1+1./(r.*k_i.^2.*(0.05*r+0.04)));
gamma_i=gamma_i_H+R.*gamma_i_L;
gamma_SH=gamma_e_SH-gamma_i;

gamma_e_SH_range=gamma_e_SH(x_range,:);
gamma_SH_range=gamma_SH(x_range,:);

figure;plot(x(x_range),max(gamma_SH(x_range,:),[],2),x(x_range),max(gamma_e_SH(x_range,:),[],2));xlabel('x');legend('\gamma','\gamma_e');
figure;[~,Ix]=max(max(gamma_e_SH_range,[],2));plot(k,gamma_SH_range(Ix,:),k,gamma_e_SH_range(Ix,:));legend('\gamma','\gamma_e');ylim([0,max(gamma_e_SH_range(Ix,:))]);xlabel('k');title(['x=',num2str(x(Ix))])
figure;plot(x(x_range),delta_T(x_range,:));xlabel('x');ylabel('\delta_T^{SH}')
[~,max_k_Index_SH]=max(gamma_e_SH(x_range,:),[],2);figure;plot(x(x_range),k(max_k_Index_SH));xlabel('x');ylabel('k_{max}')