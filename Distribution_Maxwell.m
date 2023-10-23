function f=Distribution_Maxwell(x,v,Isp,vT)
global dimx dimv dimlm
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,~]=meshgrid(v,x);
f_ini(:,:,1)=exp(-V.^2/2/vT^2)/(2*pi)^(3/2)/vT^3;
f=f_ini;