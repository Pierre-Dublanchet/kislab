%----------------------------------------------------------------%
%                                                                %
%                           KISLAB                               %
%                                                                %
% Copyright (c) (2024) Pierre Dublanchet                         %
% Author:  Pierre Dublanchet (Mines Paris PSL/Armines)           %
% Contact: pierre.dublanchet@minesparis.psl.eu                   %
% License: MIT License                                           %
%                                                                %
%----------------------------------------------------------------%
% KIS_D: deterministic inversion step
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

par.l_smooth=ksipar.REGUL*ones(1,3);       %--smoothing parameter
lx=2*ksipar.R;                             %--x length scale (m)
ly=ksipar.R;                               %--y length scale (m)
th=ksipar.TH;                              %--fault angle

%--load observations---------%
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);
par.tobs=tobs;
[nj,nt]=size(dobs);
par.nt=nt;

%--normalize observations----%
t0=max(par.tobs);        %--normalized time scale (s)
s0=max(abs(dobs'));      %--strain scale (strain unit)
u0=max(abs(dobsms'));    %--slip scale (mm)

[nj,nt]=size(dobs);
par.tobs=par.tobs/t0; %--normalized observation time
for j=1:nj
    %par.dobs(j,:)=dobs(j,:)/s0(j); %--normalize strain observations by the max of each time series
    par.dobs(j,:)=dobs(j,:)/max(s0); %--normalize strain observations by the max of all time series
end
par.dobsms=dobsms/u0; %--normalized mean slip


%--compute error covariance matrix variance (diagonal elements)----%
par.cov_strain=((max(s0)./ksipar.STD_STRAIN).^2)*ones(nj,1);
par.cov_slip=(u0/ksipar.STD_SLIP).^2;

par.cov_strain=par.cov_strain.*ksipar.GAUGE_CONF_LEV;

%--load Green's function (axial strain e_zz)-----%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);

for j=1:nj
    %gf_ezz(j,:)=gf_ezz(j,:)*u0/s0(j); %--normalized green's function
    gf_ezz(j,:)=gf_ezz(j,:)*u0/max(s0); %--normalized green's function
end
[nj0,nu]=size(gf_ezz);
par.nu=nu;
par.gf_ezz=gf_ezz;

%--normalized initial vector---%
x0=zeros(3*nu,1);

x0(1:nu)=log(ksipar.INIT_DET_DU*ones(nu,1));
x0(nu+1:2*nu)=log(ksipar.INIT_DET_T0*ones(nu,1));
x0(2*nu+1:3*nu)=log(ksipar.INIT_DET_T*ones(nu,1));

%--load matrices to compute slip gradient--%
file_gradient=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/GHmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_gradient)
par.G=G(:,:)*lx;  %--normalized gradient x operator
par.H=H(:,:)*ly;  %--normalized gradient y operator
par.P=G*G;
par.Q=H*H;
par.S=par.P+par.Q;

%--load matrix to compute spatial averages--%
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_mv)
par.MV=MV;

%--switch for data taken into account-------%
par.a_data=ksipar.SWITCH_STRAIN_DATA;         %--strain data
par.am_data=ksipar.SWITCH_SLIP_DATA;          %--average slip


%--Perform inversion (LBFGS algorithm)------------------------%
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','lbfgs','InitialHessType','identity','GoalsExactAchieve',0,'MaxIter',ksipar.NITERINV_D,'Tolx',ksipar.TOLINV_D);
[x1,fval2,exitflag,output,grad] = fminlbfgs(@(x0)slipsawcut_smooth(x0,par),x0,options);

x1(1:nu)=u0*exp(x1(1:nu));           %--slip step (mm)
x1(nu+1:3*nu)=t0*exp(x1(nu+1:3*nu)); %--onset of slip (t0, s) and ramp duration (T, s)
Jobs=fval2;                          %--value of the objective function

%--Save inversion results (vector x1)-------%
filesave=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/inv_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_ls',num2str(ksipar.REGUL),'.mat'];
save(filesave,'x1','Jobs','par');


%----------------------------------------------------------%
%--Objective function and gradient of  objective function--%
%----------------------------------------------------------%
function [f,g]=slipsawcut_smooth(x,par)



if ( nargout > 1 )
    f=0;
    g=zeros(3*par.nu,1);
    for i=1:par.nt
        [psi,phi,lambda]=basis_fn_sc_smooth(par.tobs(i),exp(x(par.nu+1:2*par.nu)),exp(x(2*par.nu+1:3*par.nu)));
        u=psi.*exp(x(1:par.nu));
        y=exp(x);

        %-----------------------------------------------%
        %--Smoothing with slip gradient norm------------%
        %-----------------------------------------------%

        %         f=f+0.5*(((par.gf_ezz*u-dobs(:,i)).*a_data)')*((par.gf_ezz*u-dobs(:,i)).*a_data)+l_smooth*(((par.G*u)')*(par.G*u)+((H*u)')*(H*u))+0.5*am_data*((MV*u-dobsms(i)).^2);
        %
        %         g(1:nu)=g(1:nu)+psi.*((par.gf_ezz')*(a_data.*(a_data.*(par.gf_ezz*u-dobs(:,i)))))...
        %                  +2*l_smooth*(psi.*((par.G')*(par.G*u))+psi.*((par.H')*(par.H*u)))...
        %                  +am_data*(MV'.*psi)*(MV*u-dobsms(i));
        %         g(nu+1:2*nu)=g(nu+1:2*nu)+(phi.*exp(x(1:nu))).*((par.gf_ezz')*(a_data.*(a_data.*(par.gf_ezz*u-dobs(:,i)))))...
        %                      +2*l_smooth*((phi.*exp(x(1:nu))).*((par.G')*(par.G*u))+(phi.*exp(x(1:nu))).*((par.H')*(par.H*u)))...
        %                      +am_data*(MV'.*(phi.*exp(x(1:nu))))*(MV*u-dobsms(i));
        %         g(2*nu+1:3*nu)=g(2*nu+1:3*nu)+(lambda.*exp(x(1:nu))).*((par.gf_ezz')*(a_data.*(a_data.*(par.gf_ezz*u-dobs(:,i)))))...
        %                      +2*l_smooth*((lambda.*exp(x(1:nu))).*((par.G')*(par.G*u))+(lambda.*exp(x(1:nu))).*((par.H')*(par.H*u)))...
        %                      +am_data*(MV'.*(lambda.*exp(x(1:nu))))*(MV*u-dobsms(i));
        %------------------------------------------%
        %--Smoothing with x gradient norm------------%
        %------------------------------------------%

        f=f+0.5*(((par.gf_ezz*u-par.dobs(:,i)).*par.a_data)')*((par.gf_ezz*u-par.dobs(:,i)).*par.a_data.*par.cov_strain)+0.5*par.am_data*par.cov_slip*((par.MV*u-par.dobsms(i)).^2)...
            +par.l_smooth(1)*(((par.G*y(1:par.nu))')*(par.G*y(1:par.nu))+((par.H*y(1:par.nu))')*(par.H*y(1:par.nu)))...
            +par.l_smooth(2)*(((par.G*y(par.nu+1:2*par.nu))')*(par.G*y(par.nu+1:2*par.nu))+((par.H*y(par.nu+1:2*par.nu))')*(par.H*y(par.nu+1:2*par.nu)))...
            +par.l_smooth(3)*(((par.G*y(2*par.nu+1:3*par.nu))')*(par.G*y(2*par.nu+1:3*par.nu))+((par.H*y(2*par.nu+1:3*par.nu))')*(par.H*y(2*par.nu+1:3*par.nu)));

        g(1:par.nu)=g(1:par.nu)+psi.*((par.gf_ezz')*(par.a_data.*(par.cov_strain.*par.a_data.*(par.gf_ezz*u-par.dobs(:,i)))))...
            +par.am_data*par.cov_slip*(par.MV'.*psi)*(par.MV*u-par.dobsms(i))...
            +2*par.l_smooth(1)*(((par.G')*(par.G*y(1:par.nu)))+((par.H')*(par.H*y(1:par.nu))));
        g(par.nu+1:2*par.nu)=g(par.nu+1:2*par.nu)+(phi.*exp(x(1:par.nu))).*((par.gf_ezz')*(par.a_data.*(par.cov_strain.*par.a_data.*(par.gf_ezz*u-par.dobs(:,i)))))...
            +par.am_data*par.cov_slip*(par.MV'.*(phi.*exp(x(1:par.nu))))*(par.MV*u-par.dobsms(i))...
            +2*par.l_smooth(2)*(((par.G')*(par.G*y(par.nu+1:2*par.nu)))+((par.H')*(par.H*y(par.nu+1:2*par.nu))));
        g(2*par.nu+1:3*par.nu)=g(2*par.nu+1:3*par.nu)+(lambda.*exp(x(1:par.nu))).*((par.gf_ezz')*(par.a_data.*(par.cov_strain.*par.a_data.*(par.gf_ezz*u-par.dobs(:,i)))))...
            +par.am_data*par.cov_slip*(par.MV'.*(lambda.*exp(x(1:par.nu))))*(par.MV*u-par.dobsms(i))...
            +2*par.l_smooth(3)*(((par.G')*(par.G*y(2*par.nu+1:3*par.nu)))+((par.H')*(par.H*y(2*par.nu+1:3*par.nu))));
        %------------------------------------------%
        %--Smoothing with laplacian norm-----------%
        %------------------------------------------%
        %         f=f+0.5*((par.gf_ezz*u-dobs(:,i))')*(par.gf_ezz*u-dobs(:,i))...
        %            +l_smooth*(((S*u)')*(S*u))...
        %            +0.5*((MV*u-dobsms(i)).^2);
        %
        %
        %
        %         g(1:nu)=g(1:nu)+psi.*((par.gf_ezz')*(par.gf_ezz*u-dobs(:,i)))...
        %                  +2*l_smooth*(psi.*((S')*S*u))...
        %                  +(MV'.*psi)*(MV*u-dobsms(i));
        %         g(nu+1:2*nu)=g(nu+1:2*nu)+(phi.*exp(x(1:nu))).*((par.gf_ezz')*(par.gf_ezz*u-dobs(:,i)))...
        %                      +2*l_smooth*((phi.*exp(x(1:nu))).*((S')*S*u))...
        %                      +(MV'.*(phi.*exp(x(1:nu))))*(MV*u-dobsms(i));
        %         g(2*nu+1:3*nu)=g(2*nu+1:3*nu)+(lambda.*exp(x(1:nu))).*((par.gf_ezz')*(par.gf_ezz*u-dobs(:,i)))...
        %                      +2*l_smooth*((lambda.*exp(x(1:nu))).*((S')*S*u))...
        %                      +(MV'.*(lambda.*exp(x(1:nu))))*(MV*u-dobsms(i));
    end
    g=g.*exp(x);
else
    f=0;
    for i=1:par.nt
        psi=slip_basis_fn_sc_smooth(par.tobs(i),exp(x(par.nu+1:2*par.nu)),exp(x(2*par.nu+1:3*par.nu)));
        u=psi.*exp(x(1:par.nu));
        y=exp(x);
        %-----------------------------------------------%
        %--Smoothing with slip gradient norm------------%
        %-----------------------------------------------%
        %f=f+0.5*((a_data.*(par.gf_ezz*u-dobs(:,i)))')*((par.gf_ezz*u-dobs(:,i)).*a_data)+l_smooth*(((par.G*u)')*(par.G*u)+((par.H*u)')*(par.H*u))+0.5*am_data*((MV*u-dobsms(i)).^2);
        %-----------------------------------------------%
        %--Smoothing with x gradient norm------------%
        %-----------------------------------------------%
        f=f+0.5*(((par.gf_ezz*u-par.dobs(:,i)).*par.a_data)')*((par.gf_ezz*u-par.dobs(:,i)).*par.a_data.*par.cov_strain)+0.5*par.am_data*par.cov_slip*((par.MV*u-par.dobsms(i)).^2)...
            +par.l_smooth(1)*(((par.G*y(1:par.nu))')*(par.G*y(1:par.nu))+((par.H*y(1:par.nu))')*(par.H*y(1:par.nu)))...
            +par.l_smooth(2)*(((par.G*y(par.nu+1:2*par.nu))')*(par.G*y(par.nu+1:2*par.nu))+((par.H*y(par.nu+1:2*par.nu))')*(par.H*y(par.nu+1:2*par.nu)))...
            +par.l_smooth(3)*(((par.G*y(2*par.nu+1:3*par.nu))')*(par.G*y(2*par.nu+1:3*par.nu))+((par.H*y(2*par.nu+1:3*par.nu))')*(par.H*y(2*par.nu+1:3*par.nu)));
        %------------------------------------------%
        %--Smoothing with laplacian norm-----------%
        %------------------------------------------%
        %         f=f+0.5*((par.gf_ezz*u-dobs(:,i))')*(par.gf_ezz*u-dobs(:,i))...
        %             +l_smooth*((S*u)')*(S*u)...
        %             +0.5*((MV*u-dobsms(i)).^2);
    end
end

end



function [psi,phi,lambda]=basis_fn_sc_smooth(t,t0,T)

psi=zeros(length(t0),1);
phi=zeros(length(t0),1);
lambda=zeros(length(t0),1);
%-------------------------%
%--Parametrisation num 1--%
%-------------------------%
% indp=find(t0<=t);
% if isempty(indp)==0
%     psi(indp)=1-exp((t0(indp)-t)./T(indp));
%     phi(indp)=(psi(indp)-1)./T(indp);
%     lambda(indp)=-phi(indp).*(t0(indp)-t)./T(indp);
% end
%-------------------------%
%--Parametrisation num 2--%
%-------------------------%
% psi=0.5*(tanh((t-t0-T)./T)+1);
% phi=-(2./T).*exp(-2*(t-t0-T)./T)./((1+exp(-2*(t-t0-T)./T)).^2);
% lambda=-(2*(t-t0)./(T.^2)).*exp(-2*(t-t0-T)./T)./((1+exp(-2*(t-t0-T)./T)).^2);
%-------------------------%
%--Parametrisation num 3--%
%-------------------------%
indp=find(t0<=t);
if isempty(indp)==0
    psi(indp)=0.5*(1-cos(2*pi*(t-t0(indp))./T(indp)));
    phi(indp)=-(pi./T(indp)).*sin(2*pi*(t-t0(indp))./T(indp));
    lambda(indp)=-(pi*(t-t0(indp))./(T(indp).^2)).*sin(2*pi*(t-t0(indp))./T(indp)); %--a 2 is missing here!!
end
indpp=find(t0+0.5*T<=t);
psi(indpp)=ones(length(indpp),1);
phi(indpp)=zeros(length(indpp),1);
lambda(indpp)=zeros(length(indpp),1);
end




