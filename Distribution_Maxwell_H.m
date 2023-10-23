function f=Distribution_Maxwell_H(x,v,Isp,xmax,vTe)
global dimx dimv dimlm
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,~]=meshgrid(v,x);
f_ini(:,:,1)=exp(-V.^2/2/vTe^2)/(2*pi)^(3/2)/vTe^3;
f_ini(x<xmax/2,:,1)=0;
f=f_ini;