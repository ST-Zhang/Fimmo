function f=Distribution_Maxwell_WKB_linear4(x,v,Isp,vTe,L)
global dimx dimv dimlm
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,X]=meshgrid(v,x);
n=(X<5)*0+...
    (X>=5&X<=L+10).*1/L.*(X-5)+...
    (X>L+10)*(L+5)/L;
f_ini(:,:,1)=n.*exp(-V.^2/2/vTe^2)/(2*pi)^(3/2)/vTe^3;
f=f_ini;