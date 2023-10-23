function f=Distribution_Bump(x,v,Isp,vTe,v_bump,amp_bump)
global dimx dimv dimlm
f_ini=zeros(dimx,dimv(Isp),dimlm(Isp));
[V,~]=meshgrid(v,x);
f_ini(:,:,1)=1/vTe^3*exp(-V.^2/2/vTe^2)+amp_bump/vTe^3*exp(-(V-v_bump).^2/2/vTe^2);
f=f_ini;