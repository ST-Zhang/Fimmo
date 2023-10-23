clear;clear global;
global dt lmax mmax dimsp flag bdy alg mass Z ne%#ok<*GVMIS>
%global x;;
restoredefaultpath;
% INPUT
ExecutePath='07.IB\01_Te=1keV,Z=41';
addpath(['E:\Fimmo\',ExecutePath])
INPUT;

%----------%
% INITIALIZATION
mmax=mmax.*(mmax<=lmax)+lmax.*(mmax>lmax);
global dimlm list_l list_m I00 I10 I11 I20 I21 I22 IlmSwitchTable
[list_l,list_m]=getList_lm(lmax,mmax);
dimlm=zeros(1,dimsp);
for Isp=1:dimsp
    dimlm(Isp)=length(list_l{Isp});
end
I00=ones(1,dimsp);
I10=zeros(1,dimsp);
I11=zeros(1,dimsp);
I20=zeros(1,dimsp);
I21=zeros(1,dimsp);
I22=zeros(1,dimsp);
for Isp=1:dimsp
    if lmax(Isp)==0
        [I10(Isp),I11(Isp),I20(Isp),I21(Isp),I22(Isp)]=deal(-1,-1,-1,-1,-1);
    elseif lmax(Isp)==1&&mmax(Isp)==0
        [I10(Isp),I11(Isp),I20(Isp),I21(Isp),I22(Isp)]=deal(2,-1,-1,-1,-1);
    elseif lmax(Isp)==1&&mmax(Isp)==1
        [I10(Isp),I11(Isp),I20(Isp),I21(Isp),I22(Isp)]=deal(2,3,-1,-1,-1);
    elseif mmax(Isp)==0
        [I10(Isp),I11(Isp),I20(Isp),I21(Isp),I22(Isp)]=deal(2,-1,3,-1,-1);
    elseif mmax(Isp)==1
        [I10(Isp),I11(Isp),I20(Isp),I21(Isp),I22(Isp)]=deal(2,3,4,5,-1);
    else
        [I10(Isp),I11(Isp),I20(Isp),I21(Isp),I22(Isp)]=deal(2,3,4,5,6);
    end
end
if flag.collision_on
    IlmSwitchTable=getIlmSwitchTable;
end

tend=Nt*dt;

global dx dimx dv dimv DifMat v_lmfiltered_2dx ar
if flag.debug_info
    disp('Initializing grids...')
end
[x,dx,dimx]=xGrid_Initialize(Nx,xmax,alg.xgridrbdy);
[v,dv,dimv,DifMat]=vGrid_Initialize(Nv,vmax);
v_lmfiltered_2dx=getvlmfiltered_2dx(v,list_l,dx);
for Isp=1:dimsp
    ar{Isp}=[0.75,ones(1,dimv(Isp)-2),0.5];
end
t=tGrid_Initialize(dt,tend);


Bx=Bx0*ones(dimx,1);
Ft=zeros(dimx,4);
if exist('Ey0','var')
    Ey0_=Ey0(x);
else
    Ey0_=0;
end
if exist('Ez0','var')
    Ez0_=Ez0(x);
else
    Ez0_=0;
end
if exist('By0','var')
    By0_=By0(x);
else
    By0_=0;
end
if exist('Bz0','var')
    Bz0_=Bz0(x);
else
    Bz0_=0;
end
Ft(:,1)=Ey0_+Bz0_;
Ft(:,2)=Ez0_-By0_;
Ft(:,3)=Ey0_-Bz0_;
Ft(:,4)=Ez0_+By0_;
clear Ey0_ Ez0_ By0_ Bz0_

global Ft_ext Ex0_ext
syms t_sym
if flag.extenal_fields
    if exist('Ey0_extenal','var')
        Ey0_ext=Ey0_extenal(x,t_sym).';
    else
        Ey0_ext=zeros(dimx,1);
    end
    if exist('Ez0_extenal','var')
        Ez0_ext=Ez0_extenal(x,t_sym).';
    else
        Ez0_ext=zeros(dimx,1);
    end
    if exist('By0_extenal','var')
        By0_ext=By0_extenal(x,t_sym).';
    else
        By0_ext=zeros(dimx,1);
    end
    if exist('Bz0_extenal','var')
        Bz0_ext=Bz0_extenal(x,t_sym).';
    else
        Bz0_ext=zeros(dimx,1);
    end
    if exist('Ex0_extenal','var')
        Ex0_ext=Ex0_extenal(x,t_sym).';
    else
        Ex0_ext=zeros(dimx,1);
    end
    Ft_ext=[Ey0_ext+1i*Ez0_ext,...
        Ey0_ext-1i*Ez0_ext,...
        Bz0_ext+1i*By0_ext,...
        Bz0_ext-1i*By0_ext];
    if isa(Ft_ext,'sym')
        Ft_ext=matlabFunction(Ft_ext);
    else
        Ft_ext=@(t)Ft_ext;
    end
    if isa(Ex0_ext,'sym')
        Ex0_ext=matlabFunction(Ex0_ext);
    else
        Ex0_ext=@(t)Ex0_ext;
    end
end

global FieldTransMat laser FieldSpreadMat1 FieldSpreadMat2
FieldTransMat=0.5*[1,1,1,1;1i,-1i,-1i,1i;1,1,-1,-1;1i,-1i,1i,-1i];
if contains(bdy.Ft,'laser')&&~flag.laser_on
    error('flag.laser==0 but the boundaries of transverse fields contain ''laser''.');
elseif ~contains(bdy.Ft,'laser')&&flag.laser_on
    error('flag.laser==1 but the boundaries of transverse fields do not contain ''laser''.');
end
if flag.laser_on
    laser.amp_L_y=laser_L(1,1);
    laser.freq_L_y=laser_L(1,2);
    laser.phase_L_y=laser_L(1,3);
    laser.duration_L_y=laser_L(1,4);
    laser.amp_L_z=laser_L(2,1);
    laser.freq_L_z=laser_L(2,2);
    laser.phase_L_z=laser_L(2,3);
    laser.duration_L_z=laser_L(2,4);
    laser.amp_R_y=laser_R(1,1);
    laser.freq_R_y=laser_R(1,2);
    laser.phase_R_y=laser_R(1,3);
    laser.duration_R_y=laser_R(1,4);
    laser.amp_R_z=laser_R(2,1);
    laser.freq_R_z=laser_R(2,2);
    laser.phase_R_z=laser_R(2,3);
    laser.duration_R_z=laser_R(2,4);
end
switch bdy.Ft
    case {'laser','outflow'}
        FieldSpreadMat1=full(spdiags([-ones(dimx,1),ones(dimx,1)],[-1,1],dimx,dimx));
        FieldSpreadMat1(end,end-2:end)=[1,-4,3];
        FieldSpreadMat2=full(spdiags([-ones(dimx,1),ones(dimx,1)],[-1,1],dimx,dimx));
        FieldSpreadMat2(1,1:3)=[-3,4,-1];
    case 'periodic'
        FieldSpreadMat1=full(spdiags([-ones(dimx,1),ones(dimx,1)],[-1,1],dimx,dimx));
        FieldSpreadMat1(1,dimx)=-1;
        FieldSpreadMat1(dimx,1)=1;
        FieldSpreadMat2=FieldSpreadMat1;
    otherwise
        error('Invalid boundary type');
end

global DetFactorMat_0 DetFactorMat_1
DetFactorMat_0=[3675/2816,-1225/2816,441/2816,-75/2816;-1733/1056,877/352,-365/352,197/1056;10/11,-18/11,10/11,-2/11;-71/528,47/176,-31/176,23/528]';
DetFactorMat_1=[35/8,-35/24,21/40,-5/56;-71/12,47/12,-31/20,23/84;5/2,-13/6,11/10,-3/14;-1/3,1/3,-1/5,1/21]';

global v_OSC_Mat
if flag.IB_heating
    v_OSC_Mat=v_OSC(x)';
else
    v_OSC_Mat=1;
end
if flag.heating
    Te_hot=(Te_hot_eV(x)/511000)';
    f0_hot=1./(2*pi*Te_hot).^1.5.*exp(-v{1}.^2/2./Te_hot);
    f0_hot=4*pi*dv(1)*f0_hot./(4*pi*dv(1)*sum(ar{1}.*f0_hot.*v{1}.^2,2));
end

if file_start_number==-1
    save([ExecutePath,'\Info.mat']);
end

global Gamma
if file_start_number==-1
    if flag.debug_info
        disp('Initializing distribution functions...')
    end
    f=Distribution_Initialize(x,v,dist_ini,density_profile_ini,xmax);
    if flag.debug_info
        disp('Normalizing distribution functions...')
    end
    f=Distribution_Normalize(x,v,f,density_profile_ini,number_ratio,bdy.f,xmax);
    f=Distribution_Neutralize(v,x,f,dist_neu_ini,Z,bdy.f);
    N=zeros(1,dimsp);
    for Isp=1:dimsp
        N(Isp)=trapz([x,xmax],tpzInt(dv(Isp),4*pi*f{Isp}([1:end,1],:,1).*v{Isp}.^2));
    end
    N0=N(1);
    save([ExecutePath,'\Info.mat'],'N0','-append');
    if max_multiple>0
        max_f=[max(real(f{1}),[],'all'),max(real(f{2}),[],'all')];
        save([ExecutePath,'\Info.mat'],'max_f','-append');
    end
    if flag.collision_on
        [Gamma,coll_medianinfo]=get_Gamma_from_f(v,f,0.5);
        save([ExecutePath,'\Info.mat'],'coll_medianinfo','-append')
        if flag.debug_info
            disp('Average value of collision parameters:')
            disp(coll_medianinfo);
        end
    end
    if flag.debug_info
        disp('Initializing Ex...')
    end
    Ex=Exfield_Initialize(f,x,v,Z,bdy.f);
    save([ExecutePath,'\f',num2str(0,'%.7d'),'.mat'],'f','N');
    save([ExecutePath,'\Field',num2str(0,'%.7d'),'.mat'],'Ex','Ft');
else
    f=load([ExecutePath,'\f',num2str(file_start_number,'%.7d'),'.mat']).f;
    N=load([ExecutePath,'\f',num2str(file_start_number,'%.7d'),'.mat']).N;
    Ex=load([ExecutePath,'\Field',num2str(file_start_number,'%.7d'),'.mat']).Ex;
    Ft=load([ExecutePath,'\Field',num2str(file_start_number,'%.7d'),'.mat']).Ft;
end
if max_multiple>0
    max_f=[max(real(f{1}),[],'all'),max(real(f{2}),[],'all')];
end

if flag.probe_on
    if file_start_number==-1
        disp('Probe is selected. If the program is interrupted, please type: save([ExecutePath,''\probe.mat''],''probe_x'',''probe'')')
        probe_xIndex=zeros(size(probe_x));
        for Iprobe=1:length(probe_x)
            [~,probe_xIndex(Iprobe)]=min(abs(probe_x(Iprobe)-x));
        end
        cell_empty=cell(size(probe_x));
        [cell_empty{:,:}]=deal(nan(1,Nt));
        probe=struct('Ex',cell_empty,'Ey',cell_empty,'Ez',cell_empty,...
            'Bx',cell_empty,'By',cell_empty,'Bz',cell_empty);
    elseif exist('probe.mat','file')
        load('probe.mat');
        probe_xIndex=zeros(size(probe_x));
        for Iprobe=1:length(probe_x)
            [~,probe_xIndex(Iprobe)]=min(abs(probe_x(Iprobe)-x));
        end
        if Nt>length(probe(1).Ex)
            cell_empty=cell(size(probe_x));
            [cell_empty{:,:}]=deal(nan(1,Nt-length(probe(1).Ex)));
            for Iprobe=1:length(probe_x)
                probe(Iprobe).Ex=[probe(Iprobe).Ex,cell_empty{Iprobe}]; %#ok<SAGROW> 
                probe(Iprobe).Ey=[probe(Iprobe).Ey,cell_empty{Iprobe}]; %#ok<SAGROW> 
                probe(Iprobe).Ez=[probe(Iprobe).Ez,cell_empty{Iprobe}]; %#ok<SAGROW> 
                probe(Iprobe).Bx=[probe(Iprobe).Bx,cell_empty{Iprobe}]; %#ok<SAGROW> 
                probe(Iprobe).By=[probe(Iprobe).By,cell_empty{Iprobe}]; %#ok<SAGROW> 
                probe(Iprobe).Bz=[probe(Iprobe).Bz,cell_empty{Iprobe}]; %#ok<SAGROW> 
            end
        end
    else
        warning('Missing file ''probe.mat''. It is recreated.')
        probe_xIndex=zeros(size(probe_x));
        for Iprobe=1:length(probe_x)
            [~,probe_xIndex(Iprobe)]=min(abs(probe_x(Iprobe)-x));
        end
        cell_empty=cell(size(probe_x));
        [cell_empty{:,:}]=deal(nan(1,Nt));
        probe=struct('Ex',cell_empty,'Ey',cell_empty,'Ez',cell_empty,...
            'Bx',cell_empty,'By',cell_empty,'Bz',cell_empty);
    end
end

if file_start_number~=-1
    save([ExecutePath,'\Info.mat'],'Nt','dt','tend','file_start_number','file_write','-append');
end
if flag.init_pause
disp('Initialization compelte. Press Enter to Start.');pause;
end

% TENSORS
global tensor JacR_ JacDL JacRDL_
if flag.read_tensor&&exist('tensor.mat','file')
    if flag.debug_info
        disp('Load tensors from files.')
    end
    load([ExecutePath,'\tensor.mat']);
else
    if flag.debug_info
        disp('Spawning tensors...')
    end
    for Isp=1:dimsp
        if flag.debug_info
            disp(['For Isp=',num2str(Isp),':'])
            disp('Spawning A...')
        end
        tensor.A{Isp}=Advection_x(Isp,bdy.f);
        if strcmp(alg.convection,'TVD')
            if flag.debug_info
                disp('Spawning H...')
            end
            tensor.H{Isp}=TVD_mm(Isp,bdy.f);
        end
        if flag.debug_info
            disp('Spawning field tensors...')
        end
        tensor.Ex{Isp}=Electric_x(v{Isp},Ivc(Isp),Isp,-Z(Isp)/mass(Isp));
        tensor.Eyz0{Isp}=Electric_yz_Iso(v{Isp},Ivc(Isp),Isp,-Z(Isp)/mass(Isp));
        if mmax(Isp)~=0
            tensor.EyzM{Isp}=Electric_yz_Ani_M(v{Isp},Isp,-Z(Isp)/mass(Isp));
            tensor.EyzP{Isp}=Electric_yz_Ani_P(v{Isp},Isp,-Z(Isp)/mass(Isp));
        end
        tensor.Byz0{Isp}=Magnetic_yz_Iso(Isp,-Z(Isp)/mass(Isp));
        if mmax(Isp)~=0
            tensor.Bx{Isp}=Magnetic_x(Isp,-Z(Isp)/mass(Isp));
            tensor.ByzM{Isp}=Magnetic_yz_Ani_M(Isp,-Z(Isp)/mass(Isp));
            tensor.ByzP{Isp}=Magnetic_yz_Ani_P(Isp,-Z(Isp)/mass(Isp));
        end
    end
    if flag.collision_on
        if flag.debug_info
            disp('Spawning collision tensors...')
        end
        for Isp1=1:dimsp
            for Isp2=1:dimsp
                if flag.collision_on_ie&&Isp1>1&&Isp2==1
                    tensor.C_Iso{Isp1,Isp2}=Collision_Iso(v,Ivc,[Isp1,Isp2],mass);
                    tensor.C_Ani_00{Isp1,Isp2}=Collision_Ani_00(v,Ivc,[Isp1,Isp2],mass);
                    tensor.C_Ani_lm{Isp1,Isp2}=Collision_Ani_lm(v,Ivc,[Isp1,Isp2],mass);
                end
                if flag.collision_on_ee&&Isp1==1&&Isp2==1
                    tensor.C_Iso{1,1}=Collision_Iso(v,Ivc,[1,1],mass);
                    tensor.C_Ani_00{1,1}=Collision_Ani_00(v,Ivc,[1,1],mass);
                    tensor.C_Ani_lm{1,1}=Collision_Ani_lm(v,Ivc,[1,1],mass);
                end
                if flag.collision_on_ii&&Isp1>1&&Isp2>1
                    tensor.C_Iso{Isp1,Isp2}=Collision_Iso(v,Ivc,[Isp1,Isp2],mass);
                    tensor.C_Ani_00{Isp1,Isp2}=Collision_Ani_00(v,Ivc,[Isp1,Isp2],mass);
                    tensor.C_Ani_lm{Isp1,Isp2}=Collision_Ani_lm(v,Ivc,[Isp1,Isp2],mass);
                end
            end
        end
        if flag.collision_on_ei
            if ischar('Te')&&strcmp(Te,'auto')
                Te=4*pi/3*sum(dv(1)*sum(ar{1}.*f{1}(:,:,1).*v{1}.^4,2))*dx*511000;
            end
            for Iion=1:dimsp-1
                tensor.C_ei{Iion}=Collision_ei(v,Iion,Te);
            end
        end
    end
    if flag.IB_heating
        for Iion=1:dimsp-1
            tensor.CIB{Iion}=Collision_IB(v,Ivc,Iion,Te);
        end
    end
    save([ExecutePath,'\tensor.mat'],'tensor','JacR_','JacDL','JacRDL_','-v7.3');
end
[tensor.Jx,tensor.Jy,tensor.Jz]=Current(v,Z);

% INVARIABLE BOUNDARY
global f_L f_R df_dt_bdy H_bdy
for Isp=1:dimsp
    df_dt_bdy{Isp}=0;
    H_bdy{Isp}=0;
end
if strcmp(bdy.f,'invariable')
    if file_start_number==-1||file_start_number==0
        if flag.debug_info
            disp('Initializing invariable boundaries of distribution functions...')
        end
        [df_dt_bdy,f_L,f_R,H_bdy]=InvarBdy_Initialize(f,alg.convection);
        save([ExecutePath,'\InvBdy.mat'],'f_L','f_R','df_dt_bdy','H_bdy');
    elseif exist('InvBdy.mat','file')
        load([ExecutePath,'\InvBdy.mat']);
    else
        error('InvBdy.mat is not detected.')
    end
end

% MAIN LOOP
fnd=attach(f,Ex,Ft);
if file_start_number==-1
    file_start_number=0;
else
    load([ExecutePath,'\Info.mat'],'N0');
end
Ifile=file_start_number;
disp('Main Loop:');
tic;
for It=file_start_number+1:Nt
    fnd=time_push(v,fnd,Bx,t(It),dt,alg.timepush);
    if flag.heating
        f_hot=sum(ar{1}.*reshape(fnd(1:dimx*dimv(1)),[dimx,dimv(1)]).*v{1}.^2,2).*f0_hot;
        fnd(1:dimx*dimv(1))=alpha*f_hot(:)+(1-alpha)*fnd(1:dimx*dimv(1));
    end
    if number_tol>0
        f{1}=reshape(fnd(1:dimx*dimv(1)*dimlm(1)),[dimx,dimv(1),dimlm(1)]);
    end
    if any(isnan(fnd),'all')...
            ||...
            (number_tol>0&&abs(trapz([x,xmax],tpzInt(dv(1),4*pi*f{1}([1:end,1],:,1).*v{1}.^2))-N0)/N0>number_tol)...
            ||...
            (max_multiple>0&&max(real(fnd(1:dimx*dimv(1)*dimlm(1))),[],'all')>max_f(1)*max_multiple)
        Ifile=Ifile+1;
        for Isp=1:dimsp
            N(Isp)=trapz([x,xmax],tpzInt(dv(Isp),4*pi*f{Isp}([1:end,1],:,1).*v{Isp}.^2));
        end
        if ~flag.notransfield && ~flag.isuniform
            [f,Ex,Ft]=detach(fnd);
        elseif flag.notransfield && ~flag.isuniform
            [f,Ex]=detach(fnd);
        elseif flag.notransfield && flag.isuniform
            f=detach(fnd);
        elseif ~flag.notransfield && flag.isuniform
            [f,Ft]=detach(fnd);
        end
        save([ExecutePath,'\f',num2str(Ifile,'%.7d'),'_error.mat'],'f','N');
        save([ExecutePath,'\Field',num2str(Ifile,'%.7d'),'_error.mat'],'Ex','Ft');
        toc;
        restoredefaultpath;
        error('Program interrupted.');
    end
    if flag.loop_pause
        disp('One Loop Finished.');pause;
    end
    Ifile=Ifile+1;
    if ~mod(Ifile,file_write)
        for Isp=1:dimsp
            N(Isp)=trapz([x,xmax],tpzInt(dv(Isp),4*pi*f{Isp}([1:end,1],:,1).*v{Isp}.^2));
        end
        if ~flag.notransfield && ~flag.isuniform
            [f,Ex,Ft]=detach(fnd);
        elseif flag.notransfield && ~flag.isuniform
            [f,Ex]=detach(fnd);
        elseif flag.notransfield && flag.isuniform
            f=detach(fnd);
        elseif ~flag.notransfield && flag.isuniform
            [f,Ft]=detach(fnd);
        end
        save([ExecutePath,'\f',num2str(Ifile,'%.7d'),'.mat'],'f','N');
        save([ExecutePath,'\Field',num2str(Ifile,'%.7d'),'.mat'],'Ex','Ft');
        disp([num2str(Ifile,'%.7d'),' saved.']);toc;
        if Ifile~=Nt
            tic;
        end
        if flag.write_pause
            disp('One Write Finished.');pause;
        end
    end
    if flag.probe_on
        if ~flag.notransfield && ~flag.isuniform
            [f,Ex,Ft]=detach(fnd);
        elseif flag.notransfield && ~flag.isuniform
            [f,Ex]=detach(fnd);
        elseif flag.notransfield && flag.isuniform
            f=detach(fnd);
        elseif ~flag.notransfield && flag.isuniform
            [f,Ft]=detach(fnd);
        end
        for Iprobe=1:length(probe_x)
            probe(Iprobe).Ex(Ifile)=Ex(probe_xIndex(Iprobe),1);
            probe(Iprobe).Ey(Ifile)=0.5*(Ft(probe_xIndex(Iprobe),1)+Ft(probe_xIndex(Iprobe),3));
            probe(Iprobe).Ez(Ifile)=0.5*(Ft(probe_xIndex(Iprobe),4)+Ft(probe_xIndex(Iprobe),2));
            probe(Iprobe).Bx(Ifile)=Bx(probe_xIndex(Iprobe),1);
            probe(Iprobe).By(Ifile)=0.5*(Ft(probe_xIndex(Iprobe),4)-Ft(probe_xIndex(Iprobe),2));
            probe(Iprobe).Bz(Ifile)=0.5*(Ft(probe_xIndex(Iprobe),1)-Ft(probe_xIndex(Iprobe),3));
        end
        if ~mod(Ifile,file_write)
            save([ExecutePath,'\probe.mat'],'probe_x','probe_xIndex','probe');
        end
    end
end
% if flag.probe_on
%     save([ExecutePath,'\probe.mat'],'probe_x','probe_xIndex','probe')
% end
restoredefaultpath;
disp('The code is running to the end.')


%----------%
% FUNCTIONS
function [list_l,list_m]=getList_lm(lmax,mmax)
global dimsp
for Isp=1:dimsp
    list_l_sp=[];
    list_m_sp=[];
    for l=0:lmax(Isp)
        for m=0:min([mmax(Isp),l])
            list_l_sp=[list_l_sp,l];  %#ok<AGROW>
            list_m_sp=[list_m_sp,m];  %#ok<AGROW>
        end
    end
    list_l{Isp}=list_l_sp; %#ok<AGROW> 
    list_m{Isp}=list_m_sp; %#ok<AGROW> 
end
end
function [Gamma,coll_medianinfo]=get_Gamma_from_f(v,f,lnL_min)
global dimx dimsp I10 I11 mass Z ne dv
n=zeros(dimx,dimsp);
T=zeros(dimx,dimsp);
for Isp=1:dimsp
    n(:,Isp)=4*pi*tpzInt(dv(Isp),f{Isp}(:,:,1).*v{Isp}.^2);
    if I10(Isp)<0
        T(:,Isp)=4*pi/3*mass(Isp)*tpzInt(dv(Isp),f{Isp}(:,:,1).*v{Isp}.^4)./n(:,Isp);
    elseif I11(Isp)<0
        T(:,Isp)=4*pi/3*mass(Isp)*tpzInt(dv(Isp),f{Isp}(:,:,1).*v{Isp}.^4)./n(:,Isp)...
            -16*pi^2/27*mass(Isp)*tpzInt(dv(Isp),f{Isp}(:,:,2).*v{Isp}.^3).^2./n(:,Isp).^2;
    else
        T(:,Isp)=4*pi/3*mass(Isp)*tpzInt(dv(Isp),f{Isp}(:,:,1).*v{Isp}.^4)./n(:,Isp)...
            -16*pi^2/27*mass(Isp)*tpzInt(dv(Isp),f{Isp}(:,:,2).*v{Isp}.^3).^2./n(:,Isp).^2 ...
            -64*pi^2/27*mass(Isp)*tpzInt(dv(Isp),real(f{Isp}(:,:,3)).*v{Isp}.^3).^2./n(:,Isp).^2 ...
            -64*pi^2/27*mass(Isp)*tpzInt(dv(Isp),imag(f{Isp}(:,:,3)).*v{Isp}.^3).^2./n(:,Isp).^2;
    end
end
n(n<0)=0;
T(n<=0|T<=0)=NaN;
n=n*ne;
T=T*5.1100e05;
Gamma=zeros(dimx,dimsp,dimsp);
for Isp1=1:dimsp
    for Isp2=1:dimsp
        if Isp1==1&&Isp2==1 %e-e
            lnL=23.5-log(n(:,1).^(1/2).*T(:,1).^(-5/4))-(1e-5+(log(T(:,1))-2).^2/16).^(1/2);
        elseif Isp1==1||Isp2==1 %e-i,i-e
            lnL=(T(:,1)>=10*Z(Isp2)^2).*(24-log(n(:,1).^(1/2)./T(:,1)))+...
                (T(:,1)<10*Z(Isp2)^2).*(23-log(n(:,1).^(1/2)*Z(Isp2).*T(:,1).^(-3/2)));
        elseif Isp1>1&&Isp2>1 %i-i
            lnL=23-log(...
                Z(Isp1)*Z(Isp2)*(mass(Isp1)+mass(Isp2))./(mass(Isp1)*T(:,Isp2)+mass(Isp2)*T(:,Isp1))...
                .*sqrt(n(:,Isp1)*Z(Isp1)^2./T(:,Isp1)+n(:,Isp2)*Z(Isp2)^2./T(:,Isp2))...
                );
        end
        lnL(lnL<lnL_min|isnan(lnL)|logical(imag(lnL)))=lnL_min;
        Gamma(:,Isp1,Isp2)=real(5.303e-19*sqrt(n(:,1))*Z(Isp1)^2*Z(Isp2)^2/mass(Isp1)^2.*lnL);
    end
end
if nargout==2
    n_=n;n_(n_<=0)=NaN;
    T_=T;T_(T_<=0)=NaN;
    coll_medianinfo.n_cm_3=median(n_,'omitnan');
    coll_medianinfo.T_eV=median(T_,'omitnan');
    omega_pe_sec=5.64e4*sqrt(coll_medianinfo.n_cm_3(1));
    coll_medianinfo.lnL_ee=23.5-log(coll_medianinfo.n_cm_3(1)^(1/2)*coll_medianinfo.T_eV(1)^(-5/4))-(1e-5+(log(coll_medianinfo.T_eV(1))-2)^2/16)^(1/2);
    coll_medianinfo.tau_ee=omega_pe_sec*3.44e5*coll_medianinfo.T_eV(1)^(3/2)/coll_medianinfo.n_cm_3(1)/coll_medianinfo.lnL_ee;
    coll_medianinfo.lnL_ei=zeros(1,dimsp-1);
    coll_medianinfo.tau_ei=zeros(1,dimsp-1);
    for Isp=2:dimsp
        if coll_medianinfo.T_eV(1)<10*Z(Isp)^2
            coll_medianinfo.lnL_ei(Isp-1)=23-log(coll_medianinfo.n_cm_3(1).^(1/2)*Z(Isp)*coll_medianinfo.T_eV(1).^(-3/2));
        elseif coll_medianinfo.T_eV(1)>10*Z(Isp)^2
            coll_medianinfo.lnL_ei(Isp-1)=24-log(coll_medianinfo.n_cm_3(1).^(1/2)*coll_medianinfo.T_eV(1).^-1);
        end
        coll_medianinfo.tau_ei(Isp-1)=coll_medianinfo.tau_ee/Z(Isp)*coll_medianinfo.lnL_ee/coll_medianinfo.lnL_ei(Isp-1);
    end
    coll_medianinfo.lnL_ii=zeros(dimsp-1,dimsp-1);
    coll_medianinfo.tau_ii=zeros(1,dimsp-1);
    for Isp1=2:dimsp
        for Isp2=2:dimsp
            mu=mass([Isp1,Isp2])/1836;
            Ti=coll_medianinfo.T_eV([Isp1,Isp2]);
            ni=coll_medianinfo.n_cm_3([Isp1,Isp2]);
            coll_medianinfo.lnL_ii(Isp1-1,Isp2-1)=23-log(...
                Z(Isp1)*Z(Isp2)*(mu(1)+mu(2))/(mu(1)*Ti(2)+mu(2)*Ti(1))*...
                sqrt(ni(1)*Z(Isp1)^2/Ti(1)+ni(2)*Z(Isp2)^2/Ti(2)));
        end
        coll_medianinfo.tau_ii(Isp1-1)=omega_pe_sec*2.09e7*Ti(1)^(3/2)/ni(1)/coll_medianinfo.lnL_ii(Isp1-1,Isp1-1)*sqrt(mass(Isp1)/1836)/Z(Isp1)^2;
    end
end
end
function [x,dx,dimx]=xGrid_Initialize(Nx,xmax,alg_xgridrbdy)
global bdy
switch alg_xgridrbdy
    case 'exclude'
        dx=xmax/Nx;
        x=(0:Nx-1)*dx;
    case 'include'
        x=linspace(0,xmax,Nx);
        dx=mean(diff(x));
    otherwise
        switch bdy.f
            case {'periodic'}
                disp('x-grid excludes x=xmax because boundary is ''periodic''.')
                dx=xmax/Nx;
                x=(0:Nx-1)*dx;
            otherwise
                disp('x-grid includes x=xmax because boundary is not ''periodic''.')
                x=linspace(0,xmax,Nx);
                dx=mean(diff(x));
        end 
end
dimx=Nx;
end
function [v,dv,dimv,DifMat]=vGrid_Initialize(Nv,vmax)
global dimsp
v=cell(1,dimsp);
dv=vmax./Nv;
dimv=Nv;
DifMat=cell(1,dimsp);
for Isp=1:dimsp
    v{Isp}=(0:Nv(Isp)-1)*dv(Isp)+dv(Isp)/2;
    DifMat{Isp}=spdiags([1/2/dv(Isp)*ones(Nv(Isp),1),-1/2/dv(Isp)*ones(Nv(Isp),1)],[1,-1],Nv(Isp),Nv(Isp));
    DifMat{Isp}(1,:)=0;
    DifMat{Isp}(end,end-2:end)=[1,-4,3]/2/dv(Isp);
end
end
function v_filtered=getvlmfiltered_2dx(v,list_l,dx)
global dimv dimlm dimsp
v_filtered=cell(1,dimsp);
for Isp=1:dimsp
    v_filtered{Isp}=zeros(1,dimv(Isp),dimlm(Isp));
    for Ilm=1:dimlm(Isp)
        l=list_l{Isp}(Ilm);
        v_filtered{Isp}(1,max(1,l):end,Ilm)=v{Isp}(max(1,l):end)/2/dx;
    end
end
end
function t=tGrid_Initialize(dt,tend)
t=0:dt:tend;
end
function f=Distribution_Initialize(x,v,dist_ini,density_profile_ini,xmax)
global dimsp list_l dimlm
for Isp=1:dimsp
    Distribution_Function=str2func(['Distribution_',dist_ini{Isp}{1}]);
    f{Isp}=Distribution_Function(x,v{Isp},Isp,dist_ini{Isp}{2:end}); %#ok<AGROW> 
    switch density_profile_ini{Isp}{1}
        case 'trapezoid'
            DPV=density_profile_ini{Isp}{2};
            [~,X]=meshgrid(v{Isp},x);
            nx=(X-DPV(1))/(DPV(2)-DPV(1)).*(X>=DPV(1)&X<DPV(2))+...
                (X>=DPV(2)&X<=DPV(3))+...
                (DPV(4)-X)/(DPV(4)-DPV(3)).*(X>DPV(3)&X<=DPV(4));
            f{Isp}=f{Isp}.*nx; %#ok<AGROW>
        case 'trapezoid-nonuniform'
            DPVx=density_profile_ini{Isp}{2};
            DPVy=density_profile_ini{Isp}{3};
            [~,X]=meshgrid(v{Isp},x);
            nx=(X>=DPVx(1)&X<=DPVx(2)).*(((DPVy(2)-DPVy(1))/(DPVx(2)-DPVx(1))*(X-DPVx(2))+DPVy(2)));
            f{Isp}=f{Isp}.*nx; %#ok<AGROW>
        case {'constant','const'}
            f{Isp}=f{Isp}*density_profile_ini{Isp}{2}; %#ok<AGROW> 
        case 'quad'
            DPVx=density_profile_ini{Isp}{2};
            DPVy=density_profile_ini{Isp}{3};
            [~,X]=meshgrid(v{Isp},x);
            nx= (X>=DPVx(1)&X<DPVx(2)) .* (DPVy(1)*(X-DPVx(1))/(DPVx(2)-DPVx(1)))+...
                (X>=DPVx(2)&X<=DPVx(3)).* (((DPVy(2)-DPVy(1))/(DPVx(3)-DPVx(2))*(X-DPVx(3))+DPVy(2)))+...
                (X>DPVx(3)&X<=DPVx(4)) .* (-DPVy(2)*(X-DPVx(4))/(DPVx(4)-DPVx(3)));
            f{Isp}=f{Isp}.*nx; %#ok<AGROW>
        case 'smooth'
            xL=density_profile_ini{Isp}{2}(1)*xmax;
            if length(density_profile_ini{Isp}{2})==1
                xR=xmax-xL;
            elseif length(density_profile_ini{Isp}{2})==2
                xR=density_profile_ini{Isp}{2}(2)*xmax;
            end
            expo=density_profile_ini{Isp}{3};
            [~,X]=meshgrid(v{Isp},x);
            nx= (X<=xL) .* sin(X/xL*pi/2).^expo+...
                (X>xL&X<xR) .* 1+...
                (X>=xR) .* cos((X-xR)/(xmax-xR)*pi/2).^expo;
            f{Isp}=f{Isp}.*nx; %#ok<AGROW>
        case {'default','neutrality'}
        otherwise
            error('Invalid Density_Profile_Type.')
    end
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);
    f{Isp}(:,1:end<l,Ilm)=0;
end
end
end
function f=Distribution_Normalize(x,v,f,density_profile,number_ratio,bdy_f,xmax)
global dimsp dx dv flag
density_profile_e=density_profile{1};
switch density_profile_e{1}
    case 'trapezoid'
        DPV=density_profile_e{2};
        Int_ne=(DPV(3)-DPV(2)+DPV(4)-DPV(1))/2;
    case 'trapezoid-nonuniform'
        DPVx=density_profile_e{2};
        DPVy=density_profile_e{3};
        Int_ne=(DPVy(1)+DPVy(2))*(DPVx(2)-DPVx(1))/2;
    case {'constant','const'}
        Int_ne=density_profile_e{2}*(x(end)+dx);
    case 'quad'
        DPVx=density_profile_e{2};
        DPVy=density_profile_e{3};
        Int_ne=((DPVx(2)-DPVx(1))*DPVy(1)+(DPVx(3)-DPVx(2))*(DPVy(1)+DPVy(2))+(DPVx(4)-DPVx(3))*DPVy(2))/2;
    case 'smooth'
        xL=density_profile_e{2}(1)*xmax;
        if length(density_profile_e{2})==1
            xR=xmax-xL;
        elseif length(density_profile_e{2})==2
            xR=density_profile_e{2}(2)*xmax;
        end
        expo=density_profile_e{3};
        Int_ne=(xL+xmax-xR)/sqrt(pi)*gamma(1/2+expo/2)/gamma(1+expo/2)+(xR-xL);
    case 'default'
        switch bdy_f
            case 'periodic'
                Int_ne=trapz([x,x(end)+dx],tpzInt(dv(1),4*pi*f{1}([1:end,1],:,1).*v{1}.^2));
            otherwise
                Int_ne=trapz(x,tpzInt(dv(1),4*pi*f{1}(:,:,1).*v{1}.^2));
        end
    otherwise
        error('Invalid Density_Profile_Type.');
end
for Isp=1:dimsp
    switch bdy_f
        case 'periodic'
            Integrated_f00=trapz([x,x(end)+dx],tpzInt(dv(Isp),4*pi*f{Isp}([1:end,1],:,1).*v{Isp}.^2));
        otherwise
            Integrated_f00=trapz(x,tpzInt(dv(Isp),4*pi*f{Isp}(:,:,1).*v{Isp}.^2));
    end
    f{Isp}=f{Isp}/Integrated_f00*(number_ratio(Isp)/number_ratio(1))*Int_ne;
    MaxImagf=abs(max(imag(f{Isp}),[],'all'));
    if MaxImagf<1e-9&&MaxImagf>0
        f{Isp}=real(f{Isp});
        if flag.debug_info
            disp(['The imagary part of f{',num2str(Isp),'} is omitted because the maximal value is only ',num2str(MaxImagf),'.']);
        end
    end
end
end
function f=Distribution_Neutralize(v,x,f,dist_neu_ini,Z,bdy_f)
global dimsp dx dv
switch dist_neu_ini{1}
    case 'none'
    case 'everywhere'
        neuIsp=dist_neu_ini{2};
        rho_other=0;
        for Isp=[1:neuIsp-1,neuIsp+1:dimsp]
            rho_other=rho_other+tpzInt(dv(Isp),Z(Isp)*4*pi*f{Isp}(:,:,1).*v{Isp}.^2);
        end
        rho_itself=tpzInt(dv(neuIsp),Z(neuIsp)*4*pi*f{neuIsp}(:,:,1).*v{neuIsp}.^2);
        f{neuIsp}=-f{neuIsp}.*rho_other./rho_itself;
        f{neuIsp}(isnan(f{neuIsp}))=0;
    case 'overall'
        neuIsp=dist_neu_ini{2};
        q_other=0;
        for Isp=[1:neuIsp-1,neuIsp+1:dimsp]
            switch bdy_f
                case 'periodic'
                    q_other=q_other+trapz([x,x(end)+dx],tpzInt(dv(Isp),Z(Isp)*4*pi*f{Isp}([1:end,1],:,1).*v{Isp}.^2));
                otherwise
                    q_other=q_other+trapz(x,tpzInt(dv(Isp),Z(Isp)*4*pi*f{Isp}(:,:,1).*v{Isp}.^2));
            end
        end
        switch bdy_f
            case 'periodic'
                q_itself=trapz([x,x(end)+dx],tpzInt(dv(neuIsp),Z(neuIsp)*4*pi*f{neuIsp}([1:end,1],:,1).*v{neuIsp}.^2));
            otherwise
                q_itself=trapz(x,tpzInt(dv(neuIsp),Z(neuIsp)*4*pi*f{neuIsp}(:,:,1).*v{neuIsp}.^2));
        end
        f{neuIsp}=-f{neuIsp}*q_other/q_itself;
    otherwise
        error('Invalid neutralize type.')
end
end
function [df_dt_bdy,f_L,f_R,H_bdy]=InvarBdy_Initialize(f0,alg_convection)
global v_lmfiltered_2dx dimsp dimx dimv dimlm list_l list_m JacRDL_ JacDL
df_dt_bdy=cell(1,dimsp);
f_L=cell(1,dimsp);
f_R=cell(1,dimsp);
H_bdy=cell(1,dimsp);
for Isp=1:dimsp
    df_dt_bdy{Isp}=zeros(size(f0{Isp}));
    f_L{Isp}=f0{Isp}(1,:,:);
    f_R{Isp}=f0{Isp}(end,:,:);
    for Ilm=1:dimlm(Isp)
        l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
        K1=(l-m)/(2*l-1);
        K2=(l+m+1)/(2*l+3);
        Ilm_lM=LM(Ilm,'l-',Isp);
        Ilm_lP=LM(Ilm,'l+',Isp);
        f_L_lM=emp20(f_L{Isp}(1,:,Ilm_lM));
        f_R_lM=emp20(f_R{Isp}(1,:,Ilm_lM));
        f_L_lP=emp20(f_L{Isp}(1,:,Ilm_lP));
        f_R_lP=emp20(f_R{Isp}(1,:,Ilm_lP));
        df_dt_bdy{Isp}(1,:,Ilm)=-v_lmfiltered_2dx{Isp}(1,:,Ilm).*(-K1*f_L_lM-K2*f_L_lP);
        df_dt_bdy{Isp}(end,:,Ilm)=-v_lmfiltered_2dx{Isp}(1,:,Ilm).*(K1*f_R_lM+K2*f_R_lP);
    end
    if strcmp(alg_convection,'TVD')
        df_dt_bdy{Isp}(1,:,:)=df_dt_bdy{Isp}(1,:,:)+v_lmfiltered_2dx{Isp}.*permute((permute(f_L{Isp},[2,3,1])*JacRDL_{Isp}),[3,1,2]);
        df_dt_bdy{Isp}(end,:,:)=df_dt_bdy{Isp}(end,:,:)+v_lmfiltered_2dx{Isp}.*permute((permute(f_R{Isp},[2,3,1])*JacRDL_{Isp}),[3,1,2]);
        H_bdy{Isp}=zeros([dimx+3,dimv(Isp),dimlm(Isp)]);
        H_bdy{Isp}(2,:,:)=permute(squeeze(-f_L{Isp})*JacDL{Isp}',[3,1,2]);
        H_bdy{Isp}(dimx+2,:,:)=permute(squeeze(f_R{Isp})*JacDL{Isp}',[3,1,2]);
    end
end
end
function Ex=Exfield_Initialize(f,x,v,Z,bdy_f)
global dimx dimsp dx dv
rho=zeros(dimx,1);
for Isp=1:dimsp
    rho=rho+Z(Isp)*4*pi*tpzInt(dv(Isp),v{Isp}.^2.*f{Isp}(:,:,1));
end
rho=rho.*(abs(rho)>1e-14);
switch bdy_f
    case 'periodic'
        e=ones(dimx,1);
        D2=full(spdiags([-e,16*e,-30*e,16*e,-e],[2,1,0,-1,-2],dimx,dimx));
        D2(1,end-1:end)=[-1,16];
        D2(2,end)=-1;
        D2(end-1,1)=-1;
        D2(end,1:2)=[16,-1];
        D2(end+1,1)=1;
        phi=(D2'*D2)\(D2'*[-rho;0])*12*dx^2;
        D1=full(spdiags([-e,8*e,-8*e,e],[2,1,-1,-2],dimx,dimx));
        D1(1,end-1:end)=[1,-8];
        D1(2,end)=1;
        D1(end-1,1)=-1;
        D1(end,1:2)=[8,-1];
        Ex=-D1*phi/12/dx;
        if max(abs(-D2(1:end-1,:)*phi/12/dx^2-rho))>1e-8
            error(['The density obtained from phi may be far from that integrated by distribution. ',...
                'Check whether neutrality is satisfied when bdy.f is ''periodic''.'])
        end
    case {'vaccum','even','invariable'}
        e=ones(dimx+4,1);
        D2=full(spdiags([-e,16*e,-30*e,16*e,-e],[2,1,0,-1,-2],dimx+4,dimx+4));
        D2(1:2,:)=zeros(2,dimx+4);
        D2(end-1:end,:)=zeros(2,dimx+4);
        D2(1,1:2)=[-1,1];
        D2(2,2:3)=[-1,1];
        D2(end-1,end-2:end-1)=[-1,1];
        D2(end,end-1:end)=[-1,1];
        D2(end+1,1)=1;
        Q=trapz(x,rho);
        phi=(D2'*D2)\(D2'*[-Q*dx/2;-Q*dx/2;-rho*12*dx^2;Q*dx/2;Q*dx/2;0]);
        D1=full(spdiags([-e,8*e,-8*e,e],[2,1,-1,-2],dimx+4,dimx+4));
        D1=D1(3:end-2,:);
        Ex=-D1*phi/12/dx;
    case 'vaccum|reflect'
        e=ones(dimx+2,1);
        D2=full(spdiags([-e,16*e,-30*e,16*e,-e],[2,1,0,-1,-2],dimx+2,dimx+2));
        D2(1:2,:)=zeros(2,dimx+2);
        D2(1,1:2)=[-1,1];
        D2(2,2:3)=[-1,1];
        D2(end-1,end-3)=-2;
        D2(end,end-2:end-1)=[-2,32];
        D2(end+1,1)=1;
        Q=trapz(x,rho);
        phi=(D2'*D2)\(D2'*[-Q*dx;-Q*dx;-rho*12*dx^2;0]);
        D1=full(spdiags([-e,8*e,-8*e,e],[2,1,-1,-2],dimx+2,dimx+2));
        D1(end-1,end-3)=0;
        D1(end,end-2:end-1)=[0,0];
        D1=D1(3:end,:);
        Ex=-D1*phi/12/dx;
    case 'reflect|vaccum'
        e=ones(dimx+2,1);
        D2=full(spdiags([-e,16*e,-30*e,16*e,-e],[2,1,0,-1,-2],dimx+2,dimx+2));
        D2(end-1:end,:)=zeros(2,dimx+2);
        D2(end-1,end-2:end-1)=[-1,1];
        D2(end,end-1:end)=[-1,1];
        D2(2,4)=-2;
        D2(1,2:3)=[32,-2];
        D2(end+1,1)=1;
        Q=trapz(x,rho);
        phi=(D2'*D2)\(D2'*[-rho*12*dx^2;Q*dx;Q*dx;0]);
        D1=full(spdiags([-e,8*e,-8*e,e],[2,1,-1,-2],dimx+2,dimx+2));
        D1(2,4)=0;
        D1(1,[2,3])=[0,0];
        D1=D1(1:end-2,:);
        Ex=-D1*phi/12/dx;
    case 'reflect'
        e=ones(dimx,1);
        D2=full(spdiags([-e,16*e,-30*e,16*e,-e],[2,1,0,-1,-2],dimx,dimx));
        D2(1,2:3)=[32,-2];
        D2(2,4)=-2;
        D2(end-1,end-3)=-2;
        D2(end,end-2:end-1)=[-2,32];
        D2(end+1,1)=1;
        phi=(D2'*D2)\(D2'*[-rho*12*dx^2;0]);
        D1=full(spdiags([-e,8*e,-8*e,e],[2,1,-1,-2],dimx,dimx));
        D1(1,[2,3])=[0,0];
        D1(2,4)=0;
        D1(end-1,end-3)=0;
        D1(end,end-2:end-1)=[0,0];
        Ex=-D1*phi/12/dx;
    otherwise
        error('Invalid boundary type.');
end
end
function fnd=attach(f,Ex,Ft)
global flag dimsp
for Isp=1:dimsp
    f{Isp}=f{Isp}(:);
end
if ~flag.notransfield && ~flag.isuniform
    fnd=[cell2mat(f');Ex;Ft(:)];
elseif flag.notransfield && ~flag.isuniform
    fnd=[cell2mat(f');Ex];
elseif flag.notransfield && flag.isuniform
    fnd=cell2mat(f');
elseif ~flag.notransfield && flag.isuniform
    fnd=[cell2mat(f');Ft(:)];
end
end
function [f,F1,F2]=detach(fnd)
global flag dimsp dimx dimv dimlm
fStartIndex=1;
f=cell(1,dimsp);
for Isp=1:dimsp
    f{Isp}=reshape(fnd(fStartIndex:fStartIndex+dimx*dimv(Isp)*dimlm(Isp)-1),[dimx,dimv(Isp),dimlm(Isp)]);
    fStartIndex=fStartIndex+dimx*dimv(Isp)*dimlm(Isp);
end
if ~flag.notransfield && ~flag.isuniform
    F1=fnd(fStartIndex:fStartIndex+dimx-1);
    F2=reshape(fnd(fStartIndex+dimx:end),[dimx,4]);
elseif flag.notransfield && ~flag.isuniform
    F1=fnd(fStartIndex:fStartIndex+dimx-1);
elseif ~flag.notransfield && flag.isuniform
    F1=reshape(fnd(1:4*dimx),[dimx,4]);
end
end
function fnd=time_push(v,fnd,Bx,t,dt,alg_timepush)
%global x DifMat;;
switch alg_timepush
    case {'midpoint','rk2'}
        k_one_fnd=FP_Discretization(v,fnd,Bx,t);
        k_two_fnd=FP_Discretization(v,fnd+dt/2*k_one_fnd,Bx,t+dt/2);
        fnd=fnd+dt*k_two_fnd;
    case 'Heun2'
        k_one_fnd=FP_Discretization(v,fnd,Bx,t);
        k_two_fnd=FP_Discretization(v,fnd+dt*k_one_fnd,Bx,t+dt);
        fnd=fnd+dt/2*(k_one_fnd+k_two_fnd);
    case 'Heun3'
        k_one_fnd=FP_Discretization(v,fnd,Bx,t);
        k_two_fnd=FP_Discretization(v,fnd+dt/3*k_one_fnd,Bx,t+dt/3);
        k_three_fnd=FP_Discretization(v,fnd+dt*(2/3)*k_two_fnd,Bx,t+dt*(2/3));
        fnd=fnd+dt/4*(k_one_fnd+3*k_three_fnd);
    case 'rk4'
        k_one_fnd=FP_Discretization(v,fnd,Bx,t);
        k_two_fnd=FP_Discretization(v,fnd+dt/2*k_one_fnd,Bx,t+dt/2);
        k_three_fnd=FP_Discretization(v,fnd+dt/2*k_two_fnd,Bx,t+dt/2);
        k_four_fnd=FP_Discretization(v,fnd+dt*k_three_fnd,Bx,t+dt);
        fnd=fnd+dt/6*(k_one_fnd+2*k_two_fnd+2*k_three_fnd+k_four_fnd);
%     case 'trapezoidal'
%         impfun=@(fnd_out)fnd_out-fnd-dt/2*(FP_Discretization(v,fnd,Bx,t)+FP_Discretization(v,fnd_out,Bx,t+dt));
%         fnd=fsolve(impfun,fnd,'Display','none');
    otherwise
        error('Invalid alg_timepush type')
end
end
function k_fnd=FP_Discretization(v,fnd,Bx,t)
global tensor lmax mmax FieldTransMat I10 I11 dx flag dimsp bdy IlmSwitchTable FieldSpreadMat1 FieldSpreadMat2 dimx dimlm Ft_ext Ex0_ext v_OSC_Mat
if ~flag.notransfield && ~flag.isuniform
    [f,Ex,Ft]=detach(fnd);
    if flag.extenal_fields
        Ft_PM=Ft*FieldTransMat+Ft_ext(t);
        Ex=Ex+Ex0_ext(t);
    else
        Ft_PM=Ft*FieldTransMat;
    end
    k_f=cell(1,dimsp);
    for Isp=1:dimsp
        if mmax(Isp)~=0
            if flag.cons_E_N
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +consN(v,contract_E(Ex,tensor.Ex{Isp},f{Isp})...
                    +real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp}))...
                    +contract_E(Ft_PM(:,2),tensor.EyzM{Isp},f{Isp})...
                    +contract_E(Ft_PM(:,1),tensor.EyzP{Isp},f{Isp}),Isp)...
                    +Bx.*(tensor.Bx{Isp}.*f{Isp})...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}))...
                    +contract_B(Ft_PM(:,4),tensor.ByzM{Isp},f{Isp})...
                    +contract_B(Ft_PM(:,3),tensor.ByzP{Isp},f{Isp});
            else
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +contract_E(Ex,tensor.Ex{Isp},f{Isp})...
                    +real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp}))...
                    +contract_E(Ft_PM(:,2),tensor.EyzM{Isp},f{Isp})...
                    +contract_E(Ft_PM(:,1),tensor.EyzP{Isp},f{Isp})...
                    +Bx.*(tensor.Bx{Isp}.*f{Isp})...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}))...
                    +contract_B(Ft_PM(:,4),tensor.ByzM{Isp},f{Isp})...
                    +contract_B(Ft_PM(:,3),tensor.ByzP{Isp},f{Isp});
            end
        else
            if flag.cons_E_N
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +consN(v,contract_E(Ex,tensor.Ex{Isp},f{Isp})...
                    +real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp})),Isp)...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}));
            else
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +contract_E(Ex,tensor.Ex{Isp},f{Isp})...
                    +real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp}))...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}));
            end
        end
    end
elseif flag.notransfield && ~flag.isuniform
    [f,Ex]=detach(fnd);
    if flag.extenal_fields
        Ex=Ex+Ex0_ext(t);
    end
    k_f=cell(1,dimsp);
    for Isp=1:dimsp
        if mmax(Isp)~=0
            if flag.cons_E_N
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +consN(v,contract_E(Ex,tensor.Ex{Isp},f{Isp}),Isp)...
                    +Bx.*(tensor.Bx{Isp}.*f{Isp});
            else
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +contract_E(Ex,tensor.Ex{Isp},f{Isp})...
                    +Bx.*(tensor.Bx{Isp}.*f{Isp});
            end
        else
            if flag.cons_E_N
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +consN(v,contract_E(Ex,tensor.Ex{Isp},f{Isp}),Isp);
            else
                k_f{Isp}=contract_A(tensor,f{Isp},Isp)...
                    +contract_E(Ex,tensor.Ex{Isp},f{Isp});
            end
        end
    end
elseif flag.notransfield && flag.isuniform
    f=detach(fnd);
    k_f=cell(1,dimsp);
    for Isp=1:dimsp
        if mmax(Isp)~=0
            k_f{Isp}=Bx.*(tensor.Bx{Isp}.*f{Isp});
        else
            k_f{Isp}=zeros(size(f{Isp}));
        end
    end
elseif ~flag.notransfield && flag.isuniform
    [f,Ft]=detach(fnd);
    if flag.extenal_fields
        Ft_PM=Ft*FieldTransMat+Ft_ext(t);
    else
        Ft_PM=Ft*FieldTransMat;
    end
    k_f=cell(1,dimsp);
    for Isp=1:dimsp
        if mmax(Isp)~=0
            if flag.cons_E_N
                k_f{Isp}=consN(v,real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp}))...
                    +contract_E(Ft_PM(:,2),tensor.EyzM{Isp},f{Isp})...
                    +contract_E(Ft_PM(:,1),tensor.EyzP{Isp},f{Isp}),Isp)...
                    +Bx.*(tensor.Bx{Isp}.*f{Isp})...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}))...
                    +contract_B(Ft_PM(:,4),tensor.ByzM{Isp},f{Isp})...
                    +contract_B(Ft_PM(:,3),tensor.ByzP{Isp},f{Isp});
            else
                k_f{Isp}=real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp}))...
                    +contract_E(Ft_PM(:,2),tensor.EyzM{Isp},f{Isp})...
                    +contract_E(Ft_PM(:,1),tensor.EyzP{Isp},f{Isp})...
                    +Bx.*(tensor.Bx{Isp}.*f{Isp})...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}))...
                    +contract_B(Ft_PM(:,4),tensor.ByzM{Isp},f{Isp})...
                    +contract_B(Ft_PM(:,3),tensor.ByzP{Isp},f{Isp});
            end
        else
            if flag.cons_E_N
                k_f{Isp}=consN(v,real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp})),Isp)...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}));
            else
                k_f{Isp}=real(contract_E(Ft_PM(:,1),tensor.Eyz0{Isp},f{Isp}))...
                    +real(contract_B(Ft_PM(:,4),tensor.Byz0{Isp},f{Isp}));
            end
        end
    end
end
if flag.collision_on
    Gamma=get_Gamma_from_f(v,f,0.5);
    if flag.collision_on_ee
        if flag.cons_ee_NE %e-e,Isotropic
            k_f{1}(:,:,1)=k_f{1}(:,:,1)+Gamma(:,1,1).*consNE(v,contract_C0(tensor.C_Iso{1,1},f{1}(:,:,1),f{1}(:,:,1)),1);
        else
            k_f{1}(:,:,1)=k_f{1}(:,:,1)+Gamma(:,1,1).*contract_C0(tensor.C_Iso{1,1},f{1}(:,:,1),f{1}(:,:,1)); 
        end
        if lmax(1)>0 %e-e,Anisotrpic
            if flag.cons_ee_P&&I10(1)>0
                k_f{1}=k_f{1}...
                    +Gamma(:,1,1).*consP(v,...
                    (contract_Clm00(tensor.C_Ani_00{1,1},f{1}(:,:,1),f{1},dimlm(1))...
                    +contract_Clmlm(tensor.C_Ani_lm{1,1},f{1},f{1}(:,:,1),1)),1);
            else
                k_f{1}=k_f{1}...
                    +Gamma(:,1,1).*...
                    (contract_Clm00(tensor.C_Ani_00{1,1},f{1}(:,:,1),f{1},dimlm(1))...
                    +contract_Clmlm(tensor.C_Ani_lm{1,1},f{1},f{1}(:,:,1),1));
            end
        end
    end
    if flag.collision_on_ei % e-i
        for Iion=1:dimsp-1
            if flag.collision_on_ei
                k_f{1}=k_f{1}+Gamma(:,1,Iion+1).*contract_Cei(tensor.C_ei{Iion},f{1},f{Iion+1});
            end
        end
    end
    if flag.collision_on_ii
        for Isp=2:dimsp
            if flag.cons_ii_NE % i1-i1,Isotropic
                k_f{Isp}(:,:,1)=k_f{Isp}(:,:,1)+Gamma(:,Isp,Isp).*consNE(v,contract_C0(tensor.C_Iso{Isp,Isp},f{Isp}(:,:,1),f{Isp}(:,:,1)),Isp); %same i-i,Isotropic,Conservation
            else
                k_f{Isp}(:,:,1)=k_f{Isp}(:,:,1)+Gamma(:,Isp,Isp).*contract_C0(tensor.C_Iso{Isp,Isp},f{Isp}(:,:,1),f{Isp}(:,:,1)); %same i-i,Isotropic
            end
            if lmax(Isp)>0
                if flag.cons_ii_P % i1-i1,Anisotropic
                    k_f{Isp}=k_f{Isp}...
                        +Gamma(:,Isp,Isp).*consP(v,...
                        (contract_Clm00(tensor.C_Ani_00{Isp,Isp},f{Isp}(:,:,1),f{Isp},dimlm(Isp))...
                        +contract_Clmlm(tensor.C_Ani_lm{Isp,Isp},f{Isp},f{Isp}(:,:,1),Isp)),Isp);
                else
                    k_f{Isp}=k_f{Isp}...
                        +Gamma(:,Isp,Isp).*...
                        (contract_Clm00(tensor.C_Ani_00{Isp,Isp},f{Isp}(:,:,1),f{Isp},dimlm(Isp))...
                        +contract_Clmlm(tensor.C_Ani_lm{Isp,Isp},f{Isp},f{Isp}(:,:,1),Isp));
                end
            end
            for Isp0=Isp+1:dimsp % i1-i2,Isotropic
                if flag.cons_i1i2_NE
                    [Ci1i2,Ci2i1]=consNE2(v,...
                        Gamma(:,Isp,Isp0).*contract_C0(tensor.C_Iso{Isp,Isp0},f{Isp}(:,:,1),f{Isp0}(:,:,1)),...
                        Gamma(:,Isp0,Isp).*contract_C0(tensor.C_Iso{Isp0,Isp},f{Isp0}(:,:,1),f{Isp}(:,:,1)),...
                        Isp,Isp0);
                    k_f{Isp}(:,:,1)=k_f{Isp}(:,:,1)+Ci1i2;
                    k_f{Isp0}(:,:,1)=k_f{Isp0}(:,:,1)+Ci2i1;
                else
                    k_f{Isp}(:,:,1)=k_f{Isp}(:,:,1)+Gamma(:,Isp,Isp0).*contract_C0(tensor.C_Iso{Isp,Isp0},f{Isp}(:,:,1),f{Isp0}(:,:,1));
                    k_f{Isp0}(:,:,1)=k_f{Isp0}(:,:,1)+Gamma(:,Isp0,Isp).*contract_C0(tensor.C_Iso{Isp0,Isp},f{Isp0}(:,:,1),f{Isp}(:,:,1));
                end
            end
            if lmax(Isp)>0 %i1-i2,Anisotropic
                for Isp0=Isp+1:dimsp
                    k_f{Isp}(:,:,IlmSwitchTable{Isp,Isp0}(1,:))=k_f{Isp}(:,:,IlmSwitchTable{Isp,Isp0}(1,:))...
                        +Gamma(:,Isp,Isp0).*contract_Clm00(tensor.C_Ani_00{Isp,Isp0}(:,:,IlmSwitchTable{Isp,Isp0}(2,:)),f{Isp}(:,:,1),f{Isp0}(:,:,IlmSwitchTable{Isp,Isp0}(2,:)),length(IlmSwitchTable{Isp,Isp0}(2,:)));
                    k_f{Isp}=k_f{Isp}...
                        +Gamma(:,Isp,Isp0).*contract_Clmlm(tensor.C_Ani_lm{Isp,Isp0},f{Isp},f{Isp0}(:,:,1),Isp);
                    if lmax(Isp0)>0
                        k_f{Isp0}(:,:,IlmSwitchTable{Isp0,Isp}(1,:))=k_f{Isp0}(:,:,IlmSwitchTable{Isp0,Isp}(1,:))...
                            +Gamma(:,Isp0,Isp).*contract_Clm00(tensor.C_Ani_00{Isp0,Isp}(:,:,IlmSwitchTable{Isp0,Isp}(2,:)),f{Isp0}(:,:,1),f{Isp}(:,:,IlmSwitchTable{Isp0,Isp}(2,:)),length(IlmSwitchTable{Isp,Isp0}(2,:)));
                        k_f{Isp0}=k_f{Isp0}...
                            +Gamma(:,Isp0,Isp).*contract_Clmlm(tensor.C_Ani_lm{Isp0,Isp},f{Isp0},f{Isp}(:,:,1),Isp0);
                    end
                end
            end
        end
    end
end
if flag.collision_on_ie %i-e
    for Isp=2:dimsp
        k_f{Isp}(:,:,1)=k_f{Isp}(:,:,1)+Gamma(:,Isp,1).*contract_C0(tensor.C_Iso{Isp,1},f{Isp}(:,:,1),f{1}(:,:,1));% i-e Isotropic
        if lmax(Isp)>0 %i-e,Anisotropic
            k_f{Isp}(:,:,IlmSwitchTable{Isp,1}(1,:))=k_f{Isp}(:,:,IlmSwitchTable{Isp,1}(1,:))...
                +Gamma(:,Isp,1).*contract_Clm00(tensor.C_Ani_00{Isp,1}(:,:,IlmSwitchTable{Isp,1}(2,:)),f{Isp}(:,:,1),f{1}(:,:,IlmSwitchTable{Isp,1}(2,:)),length(IlmSwitchTable{Isp,1}(1,:))); %i-i,Anisotropic
            k_f{Isp}=k_f{Isp}...
                +Gamma(:,Isp,1).*contract_Clmlm(tensor.C_Ani_lm{Isp,1},f{Isp},f{1}(:,:,1),Isp);
        end
    end
end
if flag.IB_heating
    for Iion=1:dimsp-1
        k_f{1}(:,:,1)=k_f{1}(:,:,1)+Gamma(:,1,Iion+1).*v_OSC_Mat.*contract_C0(tensor.CIB{Iion},f{1}(:,:,1),f{Iion+1}(:,:,1));
    end
end
if ~flag.isuniform && any(I10>0)
    k_Ex=contract_Jx(tensor.Jx,f,dimsp,I10);
else
    k_Ex=0;
end
if ~flag.notransfield
    if ~flag.isuniform
        if any(I11>0)
            Jyz=contract_Jyz(tensor.Jy,tensor.Jz,f,dimsp,I11);
            switch bdy.Ft
                case 'laser'
                    Ft_laser=laser_source(t);
                    k_Ft(:,1:2)=-Jyz-(FieldSpreadMat1*Ft(:,1:2)-[Ft_laser(1:2);zeros(dimx-1,2)])/2/dx;
                    k_Ft(:,3:4)=-Jyz+(FieldSpreadMat2*Ft(:,3:4)+[zeros(dimx-1,2);Ft_laser(3:4)])/2/dx;
                case {'periodic','outflow'}
                    k_Ft(:,1:2)=-Jyz-FieldSpreadMat1*Ft(:,1:2)/2/dx;
                    k_Ft(:,3:4)=-Jyz+FieldSpreadMat2*Ft(:,3:4)/2/dx;
                otherwise
                    error('Invalid boundary type.');
            end
        else
            switch bdy.Ft
                case 'laser'
                    Ft_laser=laser_source(t);
                    k_Ft(:,1:2)=-(FieldSpreadMat1*Ft(:,1:2)-[Ft_laser(1:2);zeros(dimx-1,2)])/2/dx;
                    k_Ft(:,3:4)=(FieldSpreadMat2*Ft(:,3:4)+[zeros(dimx-1,2);Ft_laser(3:4)])/2/dx;
                case {'periodic','outflow'}
                    k_Ft(:,1:2)=-FieldSpreadMat1*Ft(:,1:2)/2/dx;
                    k_Ft(:,3:4)=FieldSpreadMat2*Ft(:,3:4)/2/dx;
                otherwise
                    error('Invalid boundary type.');
            end
        end
    elseif any(I11>0)
        k_Ft(:,1:4)=-contract_Jyz(tensor.Jy,tensor.Jz,f,dimsp,I11);
    else
        error('When the transverse field is added, at least one I11 of the species should be defined.(有纵向场时，至少定义一种粒子的I11)')
    end
else
    k_Ft=0;
end
k_fnd=attach(k_f,k_Ex,k_Ft);
end
function Ft=laser_source(t)
global laser
Ft=zeros(1,4);
if t<laser.duration_L_y
    Ft(1)=2*laser.amp_L_y*cos(laser.freq_L_y*t+laser.phase_L_y);
end
if t<laser.duration_L_z
    Ft(2)=2*laser.amp_L_z*cos(laser.freq_L_z*t+laser.phase_L_z);
end
if t<laser.duration_R_y
    Ft(3)=2*laser.amp_R_y*cos(laser.freq_R_z*t+laser.phase_R_z);
end
if t<laser.duration_R_z
    Ft(4)=2*laser.amp_R_z*cos(laser.freq_R_y*t+laser.phase_R_y);
end
end
% function f1=cellpro(f0,dt,k_f)
% global dimsp
% f1=cell(1,dimsp);
% for Isp=1:dimsp
%     f1{Isp}=f0{Isp}+dt*k_f{Isp};
% end
% end
% function f_out=getrkf(k_fs,b,num,dt,f_in)
% global dimsp
% f_out=cell(1,dimsp);
% k_f=cell(1,dimsp);
% for Isp=1:dimsp
%     k_f{Isp}=b(1)*k_fs{1}{Isp};
%     for Inum=2:num
%         k_f{Isp}=k_f{Isp}+b(Inum)*k_fs{Inum}{Isp};
%     end
%     f_out{Isp}=f_in{Isp}+dt*k_f{Isp};
% end
% end
function result=emp20(M)
% Convert [] to 0
if isempty(M)
    result=0;
else
    result=M;
end
end
% function Iv_in_2=Iv_switch(v1,Iv_in_1,v2)
% %find the index Iv2 to make v2(Iv2) nearest to v1(Iv1)  
% [~,Iv_in_2]=min(abs(v2-v1(Iv_in_1)));
% end
function I=tpzInt(dv,F)
I=dv*(sum(F,2)-0.25*F(:,1,:)-0.5*F(:,end,:));
end
function IntegralFactor=IF(Iv0,IvL,IvR,De_M,De_P)
IntegralFactor=0.5*De_M(Iv0)*(Iv0>IvL&Iv0<=IvR)+0.5*De_P(Iv0)*(Iv0>=IvL&Iv0<IvR);
end
function Nlm_=LM(Nlm,shift,Isp)
% 返回Ilm中l或m加减1后的索引，超出范围的索引返回为[]。
global list_l list_m dimlm
N=1:dimlm(Isp);
l=list_l{Isp}(Nlm);m=list_m{Isp}(Nlm);
switch shift
    case 'l-'
        Nlm_=N(list_l{Isp}==l-1&list_m{Isp}==m);
    case 'l+'
        Nlm_=N(list_l{Isp}==l+1&list_m{Isp}==m);
    case 'm-'
        Nlm_=N(list_l{Isp}==l&list_m{Isp}==m-1);
    case 'm+'
        Nlm_=N(list_l{Isp}==l&list_m{Isp}==m+1);
    case 'l-m-'
        Nlm_=N(list_l{Isp}==l-1&list_m{Isp}==m-1);
    case 'l-m+'
        Nlm_=N(list_l{Isp}==l-1&list_m{Isp}==m+1);
    case 'l+m-'
        Nlm_=N(list_l{Isp}==l+1&list_m{Isp}==m-1);
    case 'l+m+'
        Nlm_=N(list_l{Isp}==l+1&list_m{Isp}==m+1);
end
end
function Ix_=Xs(Ix,shift,bdy_f)
% 返回索引Ix向右移动shift次后相关的新索引Ix_
global dimx
Ix_=Ix+shift;
switch bdy_f
    case 'periodic'
        if Ix_<1
            Ix_=Ix_+dimx;
        elseif Ix_>dimx
            Ix_=Ix_-dimx;
        end
    case {'vaccum','invariable'}
        if Ix_<1||Ix_>dimx
            Ix_=[];
        end
    case 'vaccum|reflect'
        if Ix_<1
            Ix_=[];
        elseif Ix_>dimx
            Ix_=2*dimx-Ix_;
        end
    case 'reflect|vaccum'
        if Ix_<1
            Ix_=2-Ix_;
        elseif Ix_>dimx
            Ix_=[];
        end
    case 'reflect'
        if Ix_<1
            Ix_=2-Ix_;
        elseif Ix_>dimx
            Ix_=2*dimx-Ix_;
        end
    case 'even'
        if Ix==1||Ix==dimx||Ix_<1||Ix_>dimx
            Ix_=[];
        end
    otherwise
        error('Invaild bdy.f type.');
end
end
function Ix0__Ix_=Xe(Ix,shift,bdy_f,Ix0,l_plus_m)
% 返回 δ(Ix0,向右移动shift次后的Ix)。
global dimx
Ix_=Ix+shift;
Ix0__Ix_=Ix0==Ix_;
switch bdy_f
    case 'periodic'
        if Ix_<1
            Ix_=Ix_+dimx;
            Ix0__Ix_=Ix0==Ix_;
        elseif Ix_>dimx
            Ix_=Ix_-dimx;
            Ix0__Ix_=Ix0==Ix_;
        end
    case {'vaccum','invariable'}
        if Ix_<1||Ix_>dimx
            Ix0__Ix_=0;
        end
    case 'vaccum|reflect'
        if Ix_<1
            Ix0__Ix_=0;
        elseif Ix_>dimx
            Ix_=2*dimx-Ix_;
            Ix0__Ix_=(-1)^(l_plus_m)*(Ix0==Ix_);
        end
    case 'reflect|vaccum'
        if Ix_<1
            Ix_=2-Ix_;
            Ix0__Ix_=(-1)^(l_plus_m)*(Ix0==Ix_);
        elseif Ix_>dimx
            Ix0__Ix_=0;
        end
    case 'reflect'
        if Ix_<1
            Ix_=2-Ix_;
            Ix0__Ix_=(-1)^(l_plus_m)*(Ix0==Ix_);
        elseif Ix_>dimx
            Ix_=2*dimx-Ix_;
            Ix0__Ix_=(-1)^(l_plus_m)*(Ix0==Ix_);
        end
    case 'even'
        if Ix==1||Ix==dimx||Ix_<1||Ix_>dimx
            Ix0__Ix_=0;
        end
    otherwise
        error('Invaild bdy.f type.');
end
end
function IlmSwitchTable=getIlmSwitchTable
global list_l list_m dimlm dimsp
IlmSwitchTable=cell([dimsp,dimsp]);
for Isp1=1:dimsp
    for Isp2=1:dimsp
        if Isp1==Isp2
            IlmSwitchTable{Isp1,Isp2}=[1:dimlm(Isp1);1:dimlm(Isp2)];
        elseif Isp2>Isp1
            for Ilm=1:dimlm(Isp1)
                Seq=1:dimlm(Isp2);
                Ilm_=Seq(list_l{Isp2}==list_l{Isp1}(Ilm)&list_m{Isp2}==list_m{Isp1}(Ilm));
                if ~isempty(Ilm_)
                    IlmSwitchTable{Isp1,Isp2}=[IlmSwitchTable{Isp1,Isp2},[Ilm;Ilm_]];
                end
            end
        else
            IlmSwitchTable{Isp1,Isp2}=IlmSwitchTable{Isp2,Isp1}([2,1],:);
        end
    end
end
end
function D_D=det0(m,n,dv)
% 对于4×4矩阵D1=[1,v1^2,v1^3,v1^4;1,v2^2,v2^3,v2^4;......v4^4],将第m行n列替换为1，非m行n列替换为0，得到D2，
% 输出D2的行列式除以D1的行列式
global DetFactorMat_0
if ~ismember(m,[1,2,3,4])
    D_D=0;
elseif n==1
    D_D=DetFactorMat_0(m,n);
else
    D_D=DetFactorMat_0(m,n)/dv^n;
end
end
function D_D=det1(m,n,dv)
% 对于4×4矩阵D1=[v1,v1^2,v1^3,v1^4;v2,v2^2,v2^3,v2^4;......v4^4],将第m行n列替换为1，非m行n列替换为0，得到D2，
% 输出D2的行列式除以D1的行列式
global DetFactorMat_1
if ~ismember(m,[1,2,3,4])
    D_D=0;
else
    D_D=DetFactorMat_1(m,n)/dv^n;
end
end
function A=Advection_x(Isp,bdy_f)
global dimx dimlm list_l list_m alg JacR_ JacDL JacRDL_
if strcmp(alg.convection,'TVD')
    JacA=zeros(dimlm(Isp),dimlm(Isp));
    for Ilm=1:dimlm(Isp)
        l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
        K1=(l-m)/(2*l-1);
        K2=(l+m+1)/(2*l+3);
        Ilm_lM=LM(Ilm,'l-',Isp);
        Ilm_lP=LM(Ilm,'l+',Isp);
        if ~isempty(Ilm_lM)
            JacA(Ilm,Ilm_lM)=K1;
        end
        if ~isempty(Ilm_lP)
            JacA(Ilm,Ilm_lP)=K2;
        end
    end
    [JacR,JacD]=eig(JacA);
    JacR_{Isp}=JacR';
    JacDL{Isp}=abs(JacD)/JacR;
    JacRDL=JacR*abs(JacD)/JacR;
    JacRDL(abs(JacRDL)<eps(10))=0;
    JacRDL(abs(imag(JacRDL))<eps(10))=real(JacRDL(abs(imag(JacRDL))<eps(10)));
    JacRDL_{Isp}=JacRDL';
    if ~isreal(JacRDL)
        warning('The convection tensor contains complex numbers.')
    end
     lengthlm0=nnz(JacRDL(1,:));
end
switch alg.convection
    case 'central'
        Ixlm_space=zeros(dimx*dimlm(Isp)*2*2,1);
        Ixlm0_space=zeros(dimx*dimlm(Isp)*2*2,1);
        Value_space=zeros(dimx*dimlm(Isp)*2*2,1);
    case 'TVD'
        try
            Ixlm_space=zeros(dimx*dimlm(Isp)*lengthlm0*3,1);
            Ixlm0_space=zeros(dimx*dimlm(Isp)*lengthlm0*3,1);
            Value_space=zeros(dimx*dimlm(Isp)*lengthlm0*3,1);
        catch
            disp('The max size of initial Advection tensor is reduced.');
            Ixlm_space=zeros(dimx*dimlm(Isp)*2*2,1);
            Ixlm0_space=zeros(dimx*dimlm(Isp)*2*2,1);
            Value_space=zeros(dimx*dimlm(Isp)*2*2,1);
        end
    otherwise
        error('Invalid alg.convection type.')
end
ISparse=1;
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    K1=(l-m)/(2*l-1);
    K2=(l+m+1)/(2*l+3);
    Ilm_lM=LM(Ilm,'l-',Isp);
    Ilm_lP=LM(Ilm,'l+',Isp);
        for Ix=1:dimx
            switch alg.convection
                case 'central'
                    Ix0_space=[Xs(Ix,-1,bdy_f),Xs(Ix,1,bdy_f)];
                    Ilm0_space=[Ilm_lM,Ilm_lP];
                case 'TVD'
                    Ix0_space=[Xs(Ix,-1,bdy_f),Ix,Xs(Ix,1,bdy_f)];
                    Ilm0_space=1:dimlm(Isp);
                otherwise
                    error('Invalid alg.convection type.')
            end
            for Ilm0=Ilm0_space
                for Ix0=Ix0_space
                    switch alg.convection
                        case 'central'
                            Ixlm_space(ISparse)=Ix+(Ilm-1)*dimx;
                            Ixlm0_space(ISparse)=Ix0+(Ilm0-1)*dimx;
                            Value_space(ISparse)...
                                =-1 ...
                                *(K1*emp20(Ilm0==Ilm_lM)+K2*emp20(Ilm0==Ilm_lP))...
                                *(Xe(Ix,1,bdy_f,Ix0,l+m)-Xe(Ix,-1,bdy_f,Ix0,l+m));
                            ISparse=ISparse+1;
                        case 'TVD'
                            Value=...
                                -1 ...
                                *(K1*emp20(Ilm0==Ilm_lM)+K2*emp20(Ilm0==Ilm_lP))...
                                *(Xe(Ix,1,bdy_f,Ix0,l+m)-Xe(Ix,-1,bdy_f,Ix0,l+m))...
                                +1 ...
                                *JacRDL(Ilm,Ilm0)...
                                *(Xe(Ix,1,bdy_f,Ix0,l+m)-2*(Ix0==Ix)+Xe(Ix,-1,bdy_f,Ix0,l+m));
                            if Value~=0
                                Ixlm_space(ISparse)=Ix+(Ilm-1)*dimx;
                                Ixlm0_space(ISparse)=Ix0+(Ilm0-1)*dimx;
                                Value_space(ISparse)=Value;
                                ISparse=ISparse+1;
                            end
                    end
                end
            end
        end
end
Ixlm_space(ISparse:end)=[];
Ixlm0_space(ISparse:end)=[];
Value_space(ISparse:end)=[];
A=sparse(Ixlm0_space,Ixlm_space,Value_space,dimx*dimlm(Isp),dimx*dimlm(Isp));
end
function H=TVD_mm(Isp,bdy_f)
global dimx dimlm list_l list_m JacDL
Ixlm_space=zeros((dimx+3)*dimlm(Isp)^2*2,1);
Ixlm0_space=zeros((dimx+3)*dimlm(Isp)^2*2,1);
Value_space=zeros((dimx+3)*dimlm(Isp)^2*2,1);
ISparse=1;
for k=1:dimlm(Isp)
    l=list_l{Isp}(k);m=list_m{Isp}(k);
    for Ix=1:dimx+3
        Ix0_space=[Xs(Ix,-1,bdy_f),Xs(Ix,-2,bdy_f)];
        for Ilm0=1:dimlm(Isp)
            for Ix0=Ix0_space
                Ixlm_space(ISparse)=Ix+(k-1)*(dimx+3);
                Ixlm0_space(ISparse)=Ix0+(Ilm0-1)*dimx;
                Value_space(ISparse)...
                    =JacDL{Isp}(k,Ilm0)...
                    *(Xe(Ix,-1,bdy_f,Ix0,l+m)-Xe(Ix,-2,bdy_f,Ix0,l+m));
                ISparse=ISparse+1;
            end
        end
    end
end
Ixlm_space(ISparse:end)=[];
Ixlm0_space(ISparse:end)=[];
Value_space(ISparse:end)=[];
H=sparse(Ixlm0_space,Ixlm_space,Value_space,dimx*dimlm(Isp),(dimx+3)*dimlm(Isp));
end
function Ex=Electric_x(v,Ivc,Isp,cmr)
global dimlm dimv list_l list_m I10 DifMat dv
Ivlm_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Ivlm0_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Value_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
ISparse=1;
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    K1=(l-m)/(2*l-1);
    K2=(l+m+1)/(2*l+3);
    Ilm_lM=LM(Ilm,'l-',Isp);
    Ilm_lP=LM(Ilm,'l+',Isp);
    for Iv=max(l,1):dimv(Isp)
        if (l>=1&&Iv>=2)||(l==0&&Iv>Ivc)
            if Iv<dimv(Isp)
                Iv0_space=[Iv-1,Iv+1];
%                 Iv0_space=[Iv-1,Iv,Iv+1];
            elseif Iv==dimv(Isp)
                Iv0_space=[dimv(Isp)-2,dimv(Isp)-1,dimv(Isp)];
            end
            for Iv0=Iv0_space
                for Ilm0=[Ilm_lM,Ilm_lP]
%                     Value=K1*v(Iv)^(l-1)*emp20(Ilm0==Ilm_lM)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(-l+1)...
%                         +K2*v(Iv)^(-l-2)*emp20(Ilm0==Ilm_lP)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(l+2);
                    Value=K1*emp20(Ilm0==Ilm_lM)*DifMat{Isp}(Iv,Iv0)*(v(Iv0)/v(Iv))^(-l+1)...
                        +K2*emp20(Ilm0==Ilm_lP)*DifMat{Isp}(Iv,Iv0)*(v(Iv0)/v(Iv))^(l+2);
%                     Value=K1*emp20(Ilm0==Ilm_lM)*(DifMat{Isp}(Iv,Iv0)+(-l-1)*(Iv0==Iv)/v(Iv))...
%                         +K2*emp20(Ilm0==Ilm_dvlP)*(DifMat{Isp}(Iv,Iv0)+(l+2)*(Iv0==Iv)/v(Iv));
                    if Value
                        Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                        Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                        Value_space(ISparse)=Value;
                        ISparse=ISparse+1;
                    end
                end
            end
        elseif l>=1&&Iv==1
            B1=1/(v(2)-v(1))/(1+(v(2)-v(1))/2/v(1));
            B2=(2*l+3)*v(1)^l/v(2)^(l+1);
            for Iv0=[1,2]
                for Ilm0=union(1,Ilm_lP)
                    if Ilm==I10(Isp)&&Ilm0==1
                        Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                        Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                        Value_space(ISparse)...
                            =K1*B1*((Iv0==2)-(Iv0==1))...
                            +K2*B2*(Iv0==2)*emp20(Ilm0==Ilm_lP);
                        ISparse=ISparse+1;
                    else
                        Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                        Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                        Value_space(ISparse)...
                            =K2*B2*(Iv0==2)*emp20(Ilm0==Ilm_lP);
                        ISparse=ISparse+1;
                    end
                end
            end
        elseif l==0&&Iv<=Ivc&&I10(Isp)>=0
            for Iv0=1:4
                d1=det1(Iv0,1,dv(Isp));
                d2=det1(Iv0,2,dv(Isp));
                d3=det1(Iv0,3,dv(Isp));
                d4=det1(Iv0,4,dv(Isp));
                Ivlm_space(ISparse)=Iv;
                Ivlm0_space(ISparse)=Iv0+(I10(Isp)-1)*dimv(Isp);
                Value_space(ISparse)...
                    =1/3*(3*d1+4*d2*v(Iv)+5*d3*v(Iv)^2+6*d4*v(Iv)^3);
                ISparse=ISparse+1;
            end
        end
    end
end
Ivlm_space(ISparse:end)=[];
Ivlm0_space(ISparse:end)=[];
Value_space(ISparse:end)=[];
Ex=sparse(Ivlm0_space,Ivlm_space,cmr*Value_space,dimv(Isp)*dimlm(Isp),dimv(Isp)*dimlm(Isp));
end
function Eyz=Electric_yz_Iso(v,Ivc,Isp,cmr)
global dimlm dimv list_l list_m I11 DifMat dv
Ivlm_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Ivlm0_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Value_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
ISparse=1;
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    if m==0
        K3=-l*(l-1)/(2*l-1);
        K4=(l+1)*(l+2)/(2*l+3);
        Ilm_lMmP=LM(Ilm,'l-m+',Isp);
        Ilm_lPmP=LM(Ilm,'l+m+',Isp);
        for Iv=max(l,1):dimv(Isp)
            if (l>=1&&Iv>=2)||(l==0&&Iv>Ivc)
                for Iv0=intersect([Iv-1,Iv+1],1:dimv(Isp))
                    for Ilm0=[Ilm_lMmP,Ilm_lPmP]
                        Value=K3*v(Iv)^(l-1)*emp20(Ilm0==Ilm_lMmP)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(-l+1)...
                            +K4*v(Iv)^(-l-2)*emp20(Ilm0==Ilm_lPmP)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(l+2);
                        if Value
                            Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                            Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                            Value_space(ISparse)=Value;
                            ISparse=ISparse+1;
                        end
                    end
                end
            elseif l>=1&&Iv==1
                B2=(2*l+3)*v(1)^l/v(2)^(l+1);
                for Ilm0=Ilm_lPmP
                    Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                    Ivlm0_space(ISparse)=2+(Ilm0-1)*dimv(Isp);
                    Value_space(ISparse)...
                        =K4*B2;
                    ISparse=ISparse+1;
                end
            elseif l==0&&Iv<=Ivc&&I11(Isp)>0
                for Iv0=1:4
                    d1=det1(Iv0,1,dv(Isp));
                    d2=det1(Iv0,2,dv(Isp));
                    d3=det1(Iv0,3,dv(Isp));
                    d4=det1(Iv0,4,dv(Isp));
                    Ivlm_space(ISparse)=Iv;
                    Ivlm0_space(ISparse)=Iv0+(I11(Isp)-1)*dimv(Isp);
                    Value_space(ISparse)...
                        =2/3*(3*d1+4*d2*v(Iv)+5*d3*v(Iv)^2+6*d4*v(Iv)^3);
                    ISparse=ISparse+1;
                end
            end
        end
    end
end
Ivlm_space(ISparse:end)=[];
Ivlm0_space(ISparse:end)=[];
Value_space(ISparse:end)=[];
Eyz=sparse(Ivlm0_space,Ivlm_space,cmr*Value_space,dimv(Isp)*dimlm(Isp),dimv(Isp)*dimlm(Isp));
end
function Eyz=Electric_yz_Ani_M(v,Isp,cmr)
global dimlm dimv list_l list_m I11 DifMat
Ivlm_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Ivlm0_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Value_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
ISparse=1;
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    if m~=0
        K5=1/(2*l-1);
        K6=-1/(2*l+3);
        Ilm_lMmM=LM(Ilm,'l-m-',Isp);
        Ilm_lPmM=LM(Ilm,'l+m-',Isp);
        for Iv=max(l,1):dimv(Isp)
            if Iv>=2
                for Iv0=intersect([Iv-1,Iv+1],1:dimv(Isp))
                    for Ilm0=[Ilm_lMmM,Ilm_lPmM] 
                        Value=K5*v(Iv)^(l-1)*emp20(Ilm0==Ilm_lMmM)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(-l+1)...
                            +K6*v(Iv)^(-l-2)*emp20(Ilm0==Ilm_lPmM)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(l+2);
                        if Value
                            Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                            Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                            Value_space(ISparse)=Value;
                            ISparse=ISparse+1;
                        end
                    end
                end
            elseif Iv==1
                B1=1/(v(2)-v(1))/(1+(v(2)-v(1))/2/v(1));
                B2=(2*l+3)*v(1)^l/v(2)^(l+1);
                for Iv0=[1,2]
                    for Ilm0=union(1,Ilm_lPmM)
                        if Ilm==I11(Isp)&&Ilm0==1
                            Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                            Ivlm0_space(ISparse)=Iv0;
                            Value_space(ISparse)...
                                =K5*B1*((Iv0==2)-(Iv0==1))...
                                +K6*B2*(Iv0==2)*emp20(Ilm0==Ilm_lPmM);
                            ISparse=ISparse+1;
                        else
                            Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                            Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                            Value_space(ISparse)...
                                =K6*B2*(Iv0==2)*emp20(Ilm0==Ilm_lPmM);
                            ISparse=ISparse+1;
                        end
                    end
                end
            end
        end
    end
end
Ivlm_space(ISparse:end)=[];
Ivlm0_space(ISparse:end)=[];
Value_space(ISparse:end)=[];
Eyz=sparse(Ivlm0_space,Ivlm_space,cmr*Value_space/2,dimv(Isp)*dimlm(Isp),dimv(Isp)*dimlm(Isp));
end
function Eyz=Electric_yz_Ani_P(v,Isp,cmr)
global dimlm dimv list_l list_m DifMat
Ivlm_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Ivlm0_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
Value_space=zeros(dimv(Isp)*dimlm(Isp)*3*2,1);
ISparse=1;
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    if m~=0
        K7=-(l-m)*(l-m-1)/(2*l-1);
        K8=(l+m+1)*(l+m+2)/(2*l+3);
        Ilm_lMmP=LM(Ilm,'l-m+',Isp);
        Ilm_lPmP=LM(Ilm,'l+m+',Isp);
        for Iv=max(l,1):dimv(Isp)
            if Iv>=2
                for Iv0=intersect([Iv-1,Iv+1],1:dimv(Isp))
                    for Ilm0=[Ilm_lMmP,Ilm_lPmP]
                        Value=K7*v(Iv)^(l-1)*emp20(Ilm0==Ilm_lMmP)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(-l+1)...
                            +K8*v(Iv)^(-l-2)*emp20(Ilm0==Ilm_lPmP)*DifMat{Isp}(Iv,Iv0)*v(Iv0)^(l+2);
                        if Value
                            Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                            Ivlm0_space(ISparse)=Iv0+(Ilm0-1)*dimv(Isp);
                            Value_space(ISparse)=Value;
                            ISparse=ISparse+1;
                        end
                    end
                end
            elseif Iv==1
                B2=(2*l+3)*v(1)^l/v(2)^(l+1);
                for Ilm0=Ilm_lPmP
                    Ivlm_space(ISparse)=Iv+(Ilm-1)*dimv(Isp);
                    Ivlm0_space(ISparse)=2+(Ilm0-1)*dimv(Isp);
                    Value_space(ISparse)...
                        =K8*B2;
                    ISparse=ISparse+1;
                end
            end
        end
    end
end
Ivlm_space(ISparse:end)=[];
Ivlm0_space(ISparse:end)=[];
Value_space(ISparse:end)=[];
Eyz=sparse(Ivlm0_space,Ivlm_space,cmr*Value_space/2,dimv(Isp)*dimlm(Isp),dimv(Isp)*dimlm(Isp));
end
function Bx=Magnetic_x(Isp,cmr)
global dimlm list_m
Bx=cmr*reshape(-1i*list_m{Isp},1,1,dimlm(Isp));
end
function Byz=Magnetic_yz_Iso(Isp,cmr)
global dimlm list_l list_m
Byz=zeros(dimlm(Isp),dimlm(Isp));
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    if m==0
        Byz(LM(Ilm,'m+',Isp),Ilm)=cmr*l*(l+1);
    end
end
end
function Byz=Magnetic_yz_Ani_M(Isp,cmr)
global dimlm list_l list_m
Byz=zeros(dimlm(Isp),dimlm(Isp));
for Ilm=1:dimlm(Isp)
    l=list_l{Isp}(Ilm);m=list_m{Isp}(Ilm);
    K9=(l-m)*(l+m+1)/2;
    if m~=0
        Byz(LM(Ilm,'m+',Isp),Ilm)=cmr*K9;
    end
end
end
function Byz=Magnetic_yz_Ani_P(Isp,cmr)
global dimlm list_m
Byz=zeros(dimlm(Isp),dimlm(Isp));
for Ilm=1:dimlm(Isp)
    m=list_m{Isp}(Ilm);
    if m~=0
        Byz(LM(Ilm,'m-',Isp),Ilm)=-1/2*cmr;
    end
end
end
function C00_00=Collision_Iso(v,Ivc,Isp12,mass)
global dimv flag dv
Isp1=Isp12(1);Isp2=Isp12(2);
if Isp1>1&&Isp2==1&&~flag.collision_on_ie % 忽略i-e情况
    C00_00=0;
elseif Isp1==1&&Isp2>1 % e-i
    C00_00=0; %各向同性近似为零
else % e-e;i1-i1;i2-i2;i-e(同种或不同种)
    C00_00=zeros([dimv(Isp1),dimv(Isp2),dimv(Isp1)]);
    mu=mass(Isp2)/mass(Isp1);
    for Iv=1:dimv(Isp1)
        v_Iv=v{Isp1}(Iv);
        Iv_p=find(v{Isp2}>=v_Iv,1);
        Iv_m=find(v{Isp2}<=v_Iv,1,'last');
        v_p=v{Isp2}(Iv_p);
        v_m=v{Isp2}(Iv_m);
        for Iv1=1:dimv(Isp1)
            if Iv<=Ivc(Isp1)
                if Iv1<=4
                    d1=det0(Iv1,2,dv(Isp1));
                    d2=det0(Iv1,3,dv(Isp1));
                    d3=det0(Iv1,4,dv(Isp1));
                    dIv1__Iv=2*d1*v_Iv+3*d2*v_Iv^2+4*d3*v_Iv^3;
                    ddIv1__Iv=2*d1+6*d2*v_Iv+12*d3*v_Iv^2;
                else
                    dIv1__Iv=0;
                    ddIv1__Iv=0;
                end
            else
                dIv1__Iv=((Iv1==Iv+1)-(Iv1==Iv-1))/2/dv(Isp1);
                ddIv1__Iv=((Iv1==Iv+1)-2*(Iv1==Iv)+(Iv1==Iv-1))/dv(Isp1)^2;
            end
            for Iv2=1:dimv(Isp2)
                v_Iv2=v{Isp2}(Iv2);
                if isempty(Iv_m)||Iv_m<=Ivc(Isp2)
                    D0=det0(Iv2,1,dv(Isp2));
                    D1=det0(Iv2,2,dv(Isp2));
                    D2=det0(Iv2,3,dv(Isp2));
                    D3=det0(Iv2,4,dv(Isp2));
                    I_2=1/3*D0*v_Iv^3+1/5*D1*v_Iv^5+1/6*D2*v_Iv^6+1/7*D3*v_Iv^7;
                    I_3=1/5*D0*v_Iv^5+1/7*D1*v_Iv^7+1/8*D2*v_Iv^8+1/9*D3*v_Iv^9;
                    Iv2__Iv=D0+D1*v_Iv^2+D2*v_Iv^3+D3*v_Iv^4;
                elseif Iv_m>Ivc(Isp2)
                    g2=1*(Iv2==1||Iv2==Iv_m)+2*(Iv2>1&&Iv2<Iv_m);
                    if v_Iv<v{Isp2}(end)
                        I_2=v_Iv2^2*dv(Isp2)/2*g2+...
                            v_Iv2^2/2*v{Isp2}(1)*(Iv2==1)+...
                            v_Iv2^2/2*(v_Iv-v_m)*(...
                            1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v_m)*(Iv2==Iv_m+1))...
                            +(Iv2==Iv_m));
                        I_3=v_Iv2^4*dv(Isp2)/2*g2+...
                            v_Iv2^4/2*v{Isp2}(1)*(Iv2==1)+...
                            v_Iv2^4/2*(v_Iv-v_m)*(...
                            1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v_m)*(Iv2==Iv_m+1))...
                            +(Iv2==Iv_m));
                        Iv2__Iv=1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v{Isp2}(Iv_m))*(Iv2==Iv_m+1));
                    elseif v_Iv>=v{Isp2}(end)
                        I_2=v_Iv2^2*dv(Isp2)/2*g2+...
                            v_Iv2^2/2*v{Isp2}(1)*(Iv2==1);
                        I_3=v_Iv2^4*dv(Isp2)/2*g2+...
                            v_Iv2^4/2*v{Isp2}(1)*(Iv2==1);
                        Iv2__Iv=0;
                    end
                end
                g1=1*(Iv2==Iv_p)+2*(Iv2>Iv_p);
                if v_Iv>v{Isp2}(end)
                    I_1=0;
                elseif v_Iv>v{Isp2}(1)
                    I_1=v_Iv2*dv(Isp2)/2*g1+...
                        v_Iv2/2*(v_p-v_Iv)*(...
                        1/dv(Isp2)*((v_p-v_Iv)*(Iv2==Iv_p-1)+(v_Iv-v{Isp2}(Iv_p-1))*(Iv2==Iv_p))...
                        +(Iv2==Iv_p));
                elseif v_Iv<=v{Isp2}(1)
                    I_1=v_Iv2*dv(Isp2)/2*g1+...
                        v_Iv2/2*(v_p-v_Iv)*(...
                        v_Iv/v{Isp2}(1)*(Iv2==1)...
                        +(Iv2==Iv_p));
                end
                C00_00(Iv1,Iv2,Iv)=4*pi/3*(...
                    dIv1__Iv*(3/mu/v_Iv^2*I_2-1/v_Iv^4*I_3+2/v_Iv*I_1)+...
                    ddIv1__Iv*(1/v_Iv^3*I_3+I_1)...
                    +3/mu*(Iv1==Iv)*Iv2__Iv);
            end
        end
    end
end
C00_00=reshape(C00_00,[dimv(Isp1)*dimv(Isp2),dimv(Isp1)]);
end
function Clm_ei=Collision_ei(v,Iion,Te)
global dimv dimlm list_l ar dv
Clm_ei=zeros([dimv(Iion+1),dimv(1),dimlm(1)]);
ve=v{1};
vi=v{Iion+1};
v_th=sqrt(Te/1000/511);
v_bar=100^(-1/3)*v_th;
v_c=6*v_th;
v_max=7*v_th;
for Ilm=1:dimlm(1)
    l=list_l{1}(Ilm);
    for Iv=max(l,1):dimv(1)
        if ve(Iv)<=v_bar
            nu=1/v_bar^3;
        elseif ve(Iv)<=v_c
            nu=1/ve(Iv)^3;
        elseif ve(Iv)<=v_max
            nu=1/ve(Iv)^3*(1-sin(pi/2*(ve(Iv)-v_c)/(v_max-v_c))^2);
        else
            nu=0;
        end
        for Iv0=1:dimv(Iion+1)
            if l>0
                Clm_ei(Iv0,Iv,Ilm)=-1/2*l*(l+1)*nu*4*pi*ar{Iion+1}(Iv0)*vi(Iv0)^2*dv(Iion+1);
            end
        end
    end
end
end
function CIB=Collision_IB(v,Ivc,Iion,Te)
global dimv alg dv ar
CIB=zeros([dimv(1),dimv(Iion+1),dimv(1)]);
ve=v{1};
dve=dv(1);
vi=v{Iion+1};
De_M_i=[vi(1),diff(vi)];
De_P_i=[diff(vi),De_M_i(end)];
De_M_e=[ve(1),diff(ve)];
De_P_e=[diff(ve),De_M_e(end)];
De_e=(De_M_e+De_P_e)/2;
v_th=sqrt(Te/1000/511);
%v_bar=100^(-1/3)*v_th;
for Iv=1:dimv(1)
    for Iv1=1:dimv(1)
        for Iv2=1:dimv(Iion+1)
            switch alg.IBheating
                case 'Langdon'
                    if Iv>=Ivc(1)&&Iv<length(ve)
                        CIB(Iv1,Iv2,Iv)=1/6/ve(Iv)^2/dve^2*...
                            (((Iv1==Iv+1)-(Iv1==Iv))*2/(ve(Iv+1)+ve(Iv))-((Iv1==Iv)-(Iv1==Iv-1))*2/(ve(Iv)+ve(Iv-1)))...
                            *4*pi*ar{Iion+1}(Iv2)*vi(Iv2)^2*dv(Iion+1);
                    elseif Iv<Ivc(1)
                        CIB(Iv1,Iv2,Iv)=1/6/ve(Ivc(1))^2/dve^2*...
                            (((Iv1==Ivc(1)+1)-(Iv1==Ivc(1)))*2/(ve(Ivc(1)+1)+ve(Ivc(1)))-((Iv1==Ivc(1))-(Iv1==Ivc(1)-1))*2/(ve(Ivc(1))+ve(Ivc(1)-1)))...
                            *4*pi*ar{Iion+1}(Iv2)*vi(Iv2)^2*dv(Iion+1);
                    end
                case 'sc'
                    if Iv>=Ivc(1)&&Iv<length(ve)
                        v_P=(ve(Iv)+ve(Iv+1))/2;
                        v_M=(ve(Iv)+ve(Iv-1))/2;
                        CIB(Iv1,Iv2,Iv)=-1/3/ve(Iv)^2/v_th^3/De_e(Iv)*...
                            (((Iv+1==Iv1)-(Iv==Iv1))*v_P^2/De_P_e(Iv)-((Iv==Iv1)-(Iv-1==Iv1))*v_M^2/De_M_e(Iv))...
                            *4*pi*vi(Iv2)^2*IF(Iv2,1,dimv(Iion+1),De_M_i,De_P_i);
                    elseif Iv<Ivc(1)
                        Iv1__0=((Iv1==1)-(Iv1==2)*ve(1)^2/ve(2)^2)/(1-ve(1)^2/ve(2)^2);
                        Iv1__d=((Iv1==2)-Iv1__0)/ve(2)^2;
                        CIB(Iv1,Iv2,Iv)=-2/v_th^3*Iv1__d...
                            *4*pi*vi(Iv2)^2*IF(Iv2,1,dimv(Iion+1),De_M_i,De_P_i);
%                         CIB(Iv1,Iv2,Iv)=-1/3/v_th^3 ...
%                            *8/(ve(2)-ve(1))^2*((Iv1==2)-(Iv1==1)) ...
%                            *4*pi*vi(Iv2)^2*IF(Iv2,1,dimv(Iion+1),De_M_i,De_P_i);不稳定
                    end
                otherwise
                    error('Invalid alg.IBheating type.');
            end
        end
    end
end
CIB=reshape(CIB,[dimv(1)*dimv(Iion+1),dimv(1)]);
end
function Clm_00=Collision_Ani_00(v,Ivc,Isp12,mass)
global dimv dimlm list_l dv
Isp1=Isp12(1);Isp2=Isp12(2);
Clm_00=zeros([dimv(Isp1),dimv(Isp2),dimv(Isp1),dimlm(Isp2)]);
mu=mass(Isp2)/mass(Isp1);
for Ilm=1:dimlm(Isp2)
    if Ilm>1
        l=list_l{Isp1}(Ilm);
        C1=(l+1)*(l+2)/(2*l+1)/(2*l+3);
        C2=-(l-1)*l/(2*l+1)/(2*l-1);
        C3=(-l*(l+1)/2-(l+1))/(2*l+1)/(2*l+3);
        C4=(-l*(l+1)/2+(l+2))/(2*l+1)/(2*l+3);
        C5=(l*(l+1)/2+(l-1))/(2*l+1)/(2*l-1);
        C6=(l*(l+1)/2-l)/(2*l+1)/(2*l-1);
        C7=(l+1)/(2*l+1);
        C8=l/(2*l+1);
        for Iv=max(l,1):dimv(Isp1)
            v_Iv=v{Isp1}(Iv);
            Iv_p=find(v{Isp2}>=v_Iv,1);
            Iv_m=find(v{Isp2}<=v_Iv,1,'last');
            v_p=v{Isp2}(Iv_p);
            v_m=v{Isp2}(Iv_m);
            for Iv1=1:dimv(Isp1)
                if Iv<=Ivc(Isp1)
                    if Iv1<=4
                        d1=det0(Iv1,2,dv(Isp1));
                        d2=det0(Iv1,3,dv(Isp1));
                        d3=det0(Iv1,4,dv(Isp1));
                        dIv1__Iv=2*d1*v_Iv+3*d2*v_Iv^2+4*d3*v_Iv^3;
                        ddIv1__Iv=2*d1+6*d2*v_Iv+12*d3*v_Iv^2;
                    else
                        dIv1__Iv=0;
                        ddIv1__Iv=0;
                    end
                else
                    dIv1__Iv=((Iv1==Iv+1)-(Iv1==Iv-1))/2/dv(Isp1);
                    ddIv1__Iv=((Iv1==Iv+1)-2*(Iv1==Iv)+(Iv1==Iv-1))/dv(Isp1)^2;
                end
                for Iv2=max(l,1):dimv(Isp2)
                    v_Iv2=v{Isp2}(Iv2);
                    if isempty(Iv_m)
                        I_1=0.5*v_Iv^2/v{Isp2}(1)*(Iv2==1)*v{Isp2}(1)^(l+2);
                        I_2=0.5*v_Iv^2/v{Isp2}(1)*(Iv2==1)*v{Isp2}(1)^(l+4);
                        Iv2__Iv=v_Iv/v{Isp2}(1)*(Iv2==1);
                    else
                        g2=1*(Iv2==1||Iv2==Iv_m)+2*(Iv2>1&&Iv2<Iv_m);
                        if v_Iv<v{Isp2}(end)
                            I_1=v_Iv2^(l+2)*dv(Isp2)/2*g2+...
                                v_Iv2^(l+2)/2*v{Isp2}(1)*(Iv2==1)+...
                                v_Iv2^(l+2)/2*(v_Iv-v_m)*(...
                                1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v_m)*(Iv2==Iv_m+1))...
                                +(Iv2==Iv_m));
                            I_2=v_Iv2^(l+4)*dv(Isp2)/2*g2+...
                                v_Iv2^(l+4)/2*v{Isp2}(1)*(Iv2==1)+...
                                v_Iv2^(l+4)/2*(v_Iv-v_m)*(...
                                1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v_m)*(Iv2==Iv_m+1))...
                                +(Iv2==Iv_m));
                            Iv2__Iv=1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v{Isp2}(Iv_m))*(Iv2==Iv_m+1));
                        elseif v_Iv>=v{Isp2}(end)
                            I_1=v_Iv2^(l+2)*dv(Isp2)/2*g2+...
                                v_Iv2^(l+2)/2*v{Isp2}(1)*(Iv2==1);
                            I_2=v_Iv2^(l+4)*dv(Isp2)/2*g2+...
                                v_Iv2^(l+4)/2*v{Isp2}(1)*(Iv2==1);
                            Iv2__Iv=0;
                        end
                    end
                    g1=1*(Iv2==Iv_p)+2*(Iv2>Iv_p);
                    if v_Iv>v{Isp2}(end)
                        J_1=0;
                        J_2=0;
                    elseif v_Iv>v{Isp2}(1)
                        J_1=v_Iv2^(-l+1)*dv(Isp2)/2*g1+...
                            v_Iv2^(-l+1)/2*(v_p-v_Iv)*(...
                            1/dv(Isp2)*((v_p-v_Iv)*(Iv2==Iv_p-1)+(v_Iv-v{Isp2}(Iv_p-1))*(Iv2==Iv_p))...
                            +(Iv2==Iv_p));
                        J_2=v_Iv2^(-l+3)*dv(Isp2)/2*g1+...
                            v_Iv2^(-l+3)/2*(v_p-v_Iv)*(...
                            1/dv(Isp2)*((v_p-v_Iv)*(Iv2==Iv_p-1)+(v_Iv-v{Isp2}(Iv_p-1))*(Iv2==Iv_p))...
                            +(Iv2==Iv_p));
                    elseif v_Iv<=v{Isp2}(1)
                        J_1=v_Iv2^(-l+1)*dv(Isp2)/2*g1+...
                            v_Iv2^(-l+1)/2*(v_p-v_Iv)*(...
                            v_Iv/v{Isp2}(1)*(Iv2==1)...
                            +(Iv2==Iv_p));
                        J_2=v_Iv2^(-l+3)*dv(Isp2)/2*g1+...
                            v_Iv2^(-l+3)/2*(v_p-v_Iv)*(...
                            v_Iv/v{Isp2}(1)*(Iv2==1)...
                            +(Iv2==Iv_p));
                    end
                    Clm_00(Iv1,Iv2,Iv,Ilm)=4*pi/mu*(Iv1==Iv)*Iv2__Iv...
                        -4*pi*(mu-1)/mu*dIv1__Iv*(C7*v_Iv^(-l-2)*I_1-C8*v_Iv^(l-1)*J_1)...
                        +2*pi*ddIv1__Iv*(C1*v_Iv^(-l-3)*I_2+C1*v_Iv^l*J_1+C2*v_Iv^(-l-1)*I_1+C2*v_Iv^(l-2)*J_2)...
                        +4*pi*dIv1__Iv*(C3*v_Iv^(-l-4)*I_2+C4*v_Iv^(l-1)*J_1+C5*v_Iv^(-l-2)*I_1+C6*v_Iv^(l-3)*J_2);
                end
            end
        end
    end
end
Clm_00=reshape(Clm_00,[dimv(Isp1)*dimv(Isp2),dimv(Isp1),dimlm(Isp2)]);
end
function Clm_lm=Collision_Ani_lm(v,Ivc,Isp12,mass)
global dimv dimlm list_l dv
Isp1=Isp12(1);Isp2=Isp12(2);
Clm_lm=zeros([dimv(Isp1),dimv(Isp2),dimv(Isp1),dimlm(Isp1)]);
mu=mass(Isp2)/mass(Isp1);
for Ilm=1:dimlm(Isp1)
    if Ilm>1
        l=list_l{Isp1}(Ilm);
        for Iv=max(l,1):dimv(Isp1)
            v_Iv=v{Isp1}(Iv);
            Iv_p=find(v{Isp2}>=v_Iv,1);
            Iv_m=find(v{Isp2}<=v_Iv,1,'last');
            v_p=v{Isp2}(Iv_p);
            v_m=v{Isp2}(Iv_m);
            for Iv1=max(l,1):dimv(Isp1)
                if l==1&&Iv==1
                    dIv1__Iv=(3*(Iv1==1)+(Iv1==2))/3/dv(Isp1);
                    ddIv1__Iv=(4*(Iv1==2)-12*(Iv1==1))/3/dv(Isp1)^2;
                else
                    dIv1__Iv=((Iv1==Iv+1)-(Iv1==Iv-1))/2/dv(Isp1);
                    ddIv1__Iv=((Iv1==Iv+1)-2*(Iv1==Iv)+(Iv1==Iv-1))/dv(Isp1)^2;
                end
                for Iv2=1:dimv(Isp2)
                    v_Iv2=v{Isp2}(Iv2);
                    if isempty(Iv_m)||Iv_m<=Ivc(Isp2)
                        D0=det0(Iv2,1,dv(Isp2));
                        D1=det0(Iv2,2,dv(Isp2));
                        D2=det0(Iv2,3,dv(Isp2));
                        D3=det0(Iv2,4,dv(Isp2));
                        I_1=1/3*D0*v_Iv^3+1/5*D1*v_Iv^5+1/6*D2*v_Iv^6+1/7*D3*v_Iv^7;
                        I_2=1/5*D0*v_Iv^5+1/7*D1*v_Iv^7+1/8*D2*v_Iv^8+1/9*D3*v_Iv^9;
                        Iv2__Iv=D0+D1*v_Iv^2+D2*v_Iv^3+D3*v_Iv^4;
                        %                             I_1=0.5*v_Iv^2/v{Isp2}(1)*(Iv2==1)*v{Isp2}(1)^2;
                        %                             I_2=0.5*v_Iv^2/v{Isp2}(1)*(Iv2==1)*v{Isp2}(1)^4;
                        %                             Iv2__Iv=v_Iv/v{Isp2}(1)*(Iv2==1);
                    else
                        g2=1*(Iv2==1||Iv2==Iv_m)+2*(Iv2>1&&Iv2<Iv_m);
                        if v_Iv<v{Isp2}(end)
                            I_1=v_Iv2^2*dv(Isp2)/2*g2+...
                                v_Iv2^2/2*v{Isp2}(1)*(Iv2==1)+...
                                v_Iv2^2/2*(v_Iv-v_m)*(...
                                1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v_m)*(Iv2==Iv_m+1))...
                                +(Iv2==Iv_m));
                            I_2=v_Iv2^4*dv(Isp2)/2*g2+...
                                v_Iv2^4/2*v{Isp2}(1)*(Iv2==1)+...
                                v_Iv2^4/2*(v_Iv-v_m)*(...
                                1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v_m)*(Iv2==Iv_m+1))...
                                +(Iv2==Iv_m));
                            Iv2__Iv=1/dv(Isp2)*((v{Isp2}(Iv_m+1)-v_Iv)*(Iv2==Iv_m)+(v_Iv-v{Isp2}(Iv_m))*(Iv2==Iv_m+1));
                        elseif v_Iv>=v{Isp2}(end)
                            I_1=v_Iv2^2*dv(Isp2)/2*g2+...
                                v_Iv2^2/2*v{Isp2}(1)*(Iv2==1);
                            I_2=v_Iv2^4*dv(Isp2)/2*g2+...
                                v_Iv2^4/2*v{Isp2}(1)*(Iv2==1);
                            Iv2__Iv=0;
                        end
                    end
                    g1=1*(Iv2==Iv_p)+2*(Iv2>Iv_p);
                    if v_Iv>v{Isp2}(end)
                        J_1=0;
                    elseif v_Iv>v{Isp2}(1)
                        J_1=v_Iv2*dv(Isp2)/2*g1+...
                            v_Iv2/2*(v_p-v_Iv)*(...
                            1/dv(Isp2)*((v_p-v_Iv)*(Iv2==Iv_p-1)+(v_Iv-v{Isp2}(Iv_p-1))*(Iv2==Iv_p))...
                            +(Iv2==Iv_p));
                    elseif v_Iv<=v{Isp2}(1)
                        J_1=v_Iv2*dv(Isp2)/2*g1+...
                            v_Iv2/2*(v_p-v_Iv)*(...
                            v_Iv/v{Isp2}(1)*(Iv2==1)...
                            +(Iv2==Iv_p));
                    end
                    Clm_lm(Iv1,Iv2,Iv,Ilm)=4*pi/mu*(Iv1==Iv)*Iv2__Iv-4*pi*(mu-1)/mu/v_Iv^2*I_1*dIv1__Iv...
                        +4*pi/3*(1/v_Iv^3*I_2+J_1)*ddIv1__Iv...
                        +4*pi/3*(-1/v_Iv^4*I_2+2/v_Iv*J_1+3/v_Iv^2*I_1)*dIv1__Iv...
                        -2*pi/3*l*(l+1)*(-1/v_Iv^5*I_2+2/v_Iv^2*J_1+3/v_Iv^3*I_1)*(Iv1==Iv);
                end
            end
        end
    end
end
Clm_lm=reshape(Clm_lm,[dimv(Isp1)*dimv(Isp2),dimv(Isp1),dimlm(Isp1)]);
end
function [Jx,Jy,Jz]=Current(v,Z)
global ar dv dimv dimsp
Jx=cell(1,dimsp);
Jy=cell(1,dimsp);
Jz=cell(1,dimsp);
for Isp=1:dimsp
    Jx{Isp}=zeros(dimv(Isp),1);
    Jy{Isp}=zeros(dimv(Isp),1);
    Jz{Isp}=zeros(dimv(Isp),1);
    for Iv0=1:dimv(Isp)
        Jx{Isp}(Iv0)=Z(Isp)*4*pi/3*ar{Isp}(Iv0)*v{Isp}(Iv0)^3*dv(Isp);
        Jy{Isp}(Iv0)=Z(Isp)*8*pi/3*ar{Isp}(Iv0)*v{Isp}(Iv0)^3*dv(Isp);
        Jz{Isp}(Iv0)=-Z(Isp)*8*pi/3*ar{Isp}(Iv0)*v{Isp}(Iv0)^3*dv(Isp);
    end
end
end
function df_dt=contract_A(TA,f,Isp)
global dimx dimv dimlm alg JacR_ v_lmfiltered_2dx df_dt_bdy H_bdy
switch alg.convection
    case 'central'
        df_dt=v_lmfiltered_2dx{Isp}.*permute(reshape(reshape(permute(f,[2,1,3]),dimv(Isp),[])*TA.A{Isp},dimv(Isp),dimx,dimlm(Isp)),[2,1,3])...
            +df_dt_bdy{Isp};
    case 'TVD'
        H_=permute(reshape(reshape(permute(f,[2,1,3]),dimv(Isp),[])*TA.H{Isp},dimv(Isp),dimx+3,dimlm(Isp)),[2,1,3])...
            +H_bdy{Isp};
        Hm = H_(1:dimx,:,:);
        H  = H_(2:dimx+1,:,:);
        Hp = H_(3:dimx+2,:,:);
        Hpp= H_(4:dimx+3,:,:);
        sH=sign(H);
        mmH=(sH==sign(Hp)).*sH.*((sH==sign(Hpp)).*min(abs(cat(4,H,Hp,Hpp)),[],4)-(sH==sign(Hm)).*min(abs(cat(4,H,Hp,Hm)),[],4));
        df_dt=v_lmfiltered_2dx{Isp}.*(...
            permute(reshape(reshape(permute(f,[2,1,3]),dimv(Isp),[])*TA.A{Isp},dimv(Isp),dimx,dimlm(Isp)),[2,1,3])...
            -reshape(reshape(mmH,dimx*dimv(Isp),[])*JacR_{Isp},size(mmH)))...
            +df_dt_bdy{Isp};
end
end
function df_dt=contract_E(E,TE,f)
df_dt=reshape(E.*(f(:,:)*TE),size(f));
end
function df_dt=contract_B(B,TB,f)
df_dt=B.*tensorprod(f,TB,3,1);
end
function df_dt=contract_C0(TC0,f0_1,f0_2)
ff=f0_1.*permute(f0_2,[1,3,2]);
df_dt=ff(:,:)*TC0;
end
function df_dt=contract_Cei(TCei,fe,fi)
df_dt=fe.*pagemtimes(fi(:,:,1),TCei);
end
function df_dt=contract_Clm00(TC,f0,Flm,lengthlm)
global dimx
df_dt=pagemtimes(reshape(f0.*permute(Flm,[1,4,2,3]),dimx,[],lengthlm),TC);
end
function df_dt=contract_Clmlm(TC,flm,F0,Isplm)
global dimx dimlm
df_dt=pagemtimes(reshape(permute(flm,[1,2,4,3]).*permute(F0,[1,3,2]),dimx,[],dimlm(Isplm)),TC);
end
function dEx_dt=contract_Jx(TJx,f,dimsp,I10)
dEx_dt=0;
for Isp=1:dimsp
    if I10(Isp)>0
        dEx_dt=dEx_dt-f{Isp}(:,:,I10(Isp))*TJx{Isp};
    end
end
end
function Jyz=contract_Jyz(TJy,TJz,f,dimsp,I11)
Jyz=0;
for Isp=1:dimsp
    if I11(Isp)>0
        Jyz=Jyz+[real(f{Isp}(:,:,I11(Isp))*TJy{Isp}),imag(f{Isp}(:,:,I11(Isp))*TJz{Isp})];
    end
end
end
function C_out=consNE(v,C_in,Isp)
global ar
A1=0.5*sum(ar{Isp}.^2.*v{Isp}.^4.*C_in.^2,2);
A2=0.5*sum(ar{Isp}.^2.*v{Isp}.^6.*C_in.^2,2);
A4=0.5*sum(ar{Isp}.^2.*v{Isp}.^8.*C_in.^2,2);
B1=sum(ar{Isp}.*v{Isp}.^2.*C_in,2);
B2=sum(ar{Isp}.*v{Isp}.^4.*C_in,2);
D=(A1.*A4-A2.^2);
lN=(A4.*B1-A2.*B2)./D;lN(isnan(lN)|isinf(lN))=0;
lE=(A1.*B2-A2.*B1)./D;lE(isnan(lE)|isinf(lE))=0;
C_out=C_in-1/2*C_in.^2.*ar{Isp}.*v{Isp}.^2.*(lN+v{Isp}.^2.*lE);
end
function [C1_out,C2_out]=consNE2(v,C1_in,C2_in,Isp1,Isp2)
global ar mass dv
A11=sum(ar{Isp1}.^2.*v{Isp1}.^4.*C1_in.^2,2)*2*pi*dv(Isp1)^2;
A22=sum(ar{Isp2}.^2.*v{Isp2}.^4.*C2_in.^2,2)*2*pi*dv(Isp2)^2;
A31=sum(ar{Isp1}.^2.*v{Isp1}.^6.*C1_in.^2,2)*pi*dv(Isp1)^2*mass(Isp1);
A32=sum(ar{Isp2}.^2.*v{Isp2}.^6.*C2_in.^2,2)*pi*dv(Isp2)^2*mass(Isp2);
A33=sum(ar{Isp1}.^2.*v{Isp1}.^8.*C1_in.^2,2)*pi/2*dv(Isp1)^2*mass(Isp1)^2+sum(ar{Isp2}.^2.*v{Isp2}.^8.*C2_in.^2,2)*pi/2*dv(Isp2)^2*mass(Isp2)^2;
B1=sum(ar{Isp1}.*v{Isp1}.^2.*C1_in,2)*dv(Isp1);
B2=sum(ar{Isp2}.*v{Isp2}.^2.*C2_in,2)*dv(Isp2);
B3=sum(ar{Isp1}.*v{Isp1}.^4.*C1_in,2)/2*dv(Isp1)*mass(Isp1)+sum(ar{Isp2}.*v{Isp2}.^4.*C2_in,2)/2*dv(Isp2)*mass(Isp2);
D=A11.*A22.*A33-A22.*A31.^2-A11.*A32.^2;
lN1=(-B1.*A32.^2+A31.*B2.*A32+A22.*A33.*B1-A22.*A31.*B3)./D;lN1(isnan(lN1)|isinf(lN1))=0;
lN2=(-B2.*A31.^2+A32.*B1.*A31+A11.*A33.*B2-A11.*A32.*B3)./D;lN2(isnan(lN2)|isinf(lN2))=0;
lE=(A11.*A22.*B3-A11.*A32.*B2-A22.*A31.*B1)./D;lE(isnan(lE)|isinf(lE))=0;
C1_out=C1_in-2*pi*dv(Isp1)*lN1.*C1_in.^2.*ar{Isp1}.*v{Isp1}.^2-pi*dv(Isp1)*mass(Isp1)*lE.*C1_in.^2.*ar{Isp1}.*v{Isp1}.^4;
C2_out=C2_in-2*pi*dv(Isp2)*lN2.*C2_in.^2.*ar{Isp2}.*v{Isp2}.^2-pi*dv(Isp2)*mass(Isp2)*lE.*C2_in.^2.*ar{Isp2}.*v{Isp2}.^4;
end
function C_in=consP(v,C_in,Isp)
global ar
lP=sum(ar{Isp}.*v{Isp}.^3.*C_in(:,:,2),2)./sum(ar{Isp}.^2.*v{Isp}.^6.*C_in(:,:,2).^2,2);lP(isnan(lP)|isinf(lP))=0;
C_in(:,:,2)=C_in(:,:,2)-lP.*C_in(:,:,2).^2.*ar{Isp}.*v{Isp}.^3;
end
function C_in=consN(v,C_in,Isp)
global ar
lN=sum(ar{Isp}.*v{Isp}.^2.*C_in(:,:,1),2)./sum(ar{Isp}.^2.*v{Isp}.^4.*C_in(:,:,1).^2,2);lN(isnan(lN)|isinf(lN))=0;
C_in(:,:,1)=C_in(:,:,1)-lN.*C_in(:,:,1).^2.*ar{Isp}.*v{Isp}.^2;
end