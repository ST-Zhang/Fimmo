function f=Distribution_Maxwell_linear(x,v,Isp,vTe)
global dimx dimv dimlm
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,X]=meshgrid(v,x);
n=(X<5)*0+...
    (X>=5&X<=35).*7/30.*(X-5)+...
    (X>35)*0;
f_ini(:,:,1)=n.*exp(-V.^2/2/vTe^2)/(2*pi)^(3/2)/vTe^3;
f=f_ini;