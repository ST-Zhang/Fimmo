function f=Distribution_TwoStream2(x,v,Isp,vT,vd,amp_noise,xmax)
global dimx dimv dimlm I10 I20
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,X]=meshgrid(v,x);
f_ini(:,:,1)=exp(-(V-vd).^2/2/vT^2)/(2*pi)^(3/2);
f_ini(:,:,I20(Isp))=exp(-(V-vd).^2/2/vT^2)/(2*pi)^(3/2)*2;
f_ini(:,:,I10(Isp))=-exp(-(V-vd).^2/2/vT^2)/(2*pi)^(3/2)*amp_noise.*sin(2*pi*X/xmax);
f=f_ini;