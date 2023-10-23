function f=Distribution_heat2(x,v,Isp,xmax,n1,vT1,vT2)
global dimx dimv dimlm
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,X]=meshgrid(v,x);
vT=(vT1+vT2)/2+(vT1-vT2)/2*sin(2*pi/xmax*X);
n=n1;
f_ini(:,:,1)=n./(2*pi*vT.^2).^(3/2).*exp(-V.^2/2./vT.^2);
f=f_ini;