% folder_name={...
%     '07.1.Nx=1024';
%     '07.2.续算'
%     '07.3.续算'
%     };
% file_space={...
%     0:500:990000;
%     0:500:990000;
%     0:500:375500
%     };
folder_name={...
    '02.0加离子碰撞'
    };
file_space={...
    0:200:140000
    };
xrange=[0,xmax];%[45,55];%[0,xmax]
k0=[8,14];
kmin=4.0;
kmax=50;
%-----------------------------------------------
N_folder=length(folder_name);
dt_space=zeros(1,N_folder);
Ex_t=cell(1,N_folder);
time=cell(1,N_folder);
time0=0;
for Ifolder=1:N_folder
    cd(folder_name{Ifolder});
    load('Info.mat')
    dt_space(Ifolder)=dt;
    time{Ifolder}=(file_space{Ifolder}-file_space{Ifolder}(1))*dt+time0;
    time0=time{Ifolder}(end);
    Ex_t{Ifolder}=zeros(dimx,length(file_space{Ifolder}));
    for Ifile=1:length(file_space{Ifolder})
        load(['Field',num2str(file_space{Ifolder}(Ifile),'%.6d'),'.mat']);
        Ex_t{Ifolder}(:,Ifile)=Ex;
    end
    cd('..')
end
Ex_t=cell2mat(Ex_t);
time=cell2mat(time);
[~,xrange_left]=min(abs(xrange(1)-x));
[~,xrange_right]=min(abs(xrange(end)-x));
Ek_t=zeros(xrange_right-xrange_left+1,size(Ex_t,2));
for It=1:size(Ex_t,2)
    Ek_t(:,It)=abs(fft(Ex_t(xrange_left:xrange_right,It)))/size(Ex_t,1);
end
k_space=2*pi/dx*(0:size(Ek_t,1)/2-1)/size(Ek_t,1);
[~,Ik0_left]=min(abs(k0(1)-k_space));
[~,Ik0_right]=min(abs(k0(end)-k_space));
Ex_k0_t=mean(Ek_t(Ik0_left:Ik0_right,:),1);

figure('Color',[1 1 1]);
plot(time,Ex_k0_t);
figure('Color',[1 1 1]);
pcolor(k_space(k_space>=kmin&k_space<=kmax),time,Ek_t(k_space>=kmin&k_space<=kmax,:)');
cmap=colormap;
cmap(1,:)=[1 1 1];
colormap(cmap)
xlabel('k');ylabel('t');colormap;
set(gca,'ticklength',[0.03,0.06],'xminortick','on','yminortick','on','linewidth',2.0,'fontsize',20,'fontweight','bold');
shading interp