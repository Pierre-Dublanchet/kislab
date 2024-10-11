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
% KIS_MCMC: mcmc iterations (bayesian step).
%           Generates mcmc models chain
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--number of result files
nsimgroup=floor(ksipar.NITER_MCMC_TOTAL/ksipar.NITER_MCMC);


%--load deterministic inversion results (x1)
if isempty(ksipar.METH_INIT_MCMC)==0
    file_invd=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/inv_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.METH_INIT_MCMC,'_ls',num2str(ksipar.REGUL),'.mat'];
    load(file_invd)
end

%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);

%--undersample data------%
[par.nj,nsdraw]=size(dobs);
if length(tobs)>ksipar.MCMC_TIME_SAMP+2
    ttobs=[0:max(tobs)/(ksipar.MCMC_TIME_SAMP-1):max(tobs)];
    ddobsms=interp1(tobs,dobsms,ttobs);
    ddobs=zeros(par.nj,ksipar.MCMC_TIME_SAMP);
    for j=1:par.nj
        ddobs(j,:)=interp1(tobs,dobs(j,:),ttobs);
    end
    tobs=ttobs;
    dobs=ddobs;
    dobsms=ddobsms;
end

%-normalize data------------------%
par.nt=length(tobs);
par.tobs=tobs/max(tobs);
par.s0=max(abs(dobs'));   %--strain scale
u0=max(abs(dobsms));  %--slip scale
par.dobs=zeros(size(dobs));
par.dobsms=zeros(size(dobsms));
for j=1:par.nj
    %par.dobs(j,:)=dobs(j,:)/par.s0(j); %--normalize strain observations by the max of each time series
    par.dobs(j,:)=dobs(j,:)/max(par.s0); %--normalize strain observations by the max of all time series
end
par.dobsms=dobsms/u0; %--normalized mean slip

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,par.nu]=size(gf_ezz);
for j=1:par.nj
    %par.gf_ezz(j,:)=gf_ezz(j,:)*u0/par.s0(j); %--normalized stiffness
    par.gf_ezz(j,:)=gf_ezz(j,:)*u0/max(par.s0); %--normalized stiffness
end


%--average  operator matrix--%
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_mv)
par.MV=MV;

%--compute error covariance matrix variance (diagonal elements)----%
par.cov_strain=((max(par.s0)./ksipar.STD_STRAIN).^2)*ones(par.nj,1);
par.cov_slip=(u0/ksipar.STD_SLIP).^2;

par.cov_strain=par.cov_strain.*ksipar.GAUGE_CONF_LEV;

%--define std of exploration steps
dm=[ksipar.MCMC_STD_STEP_DU*ones(par.nu,1);ksipar.MCMC_STD_STEP_T0*ones(par.nu,1);ksipar.MCMC_STD_STEP_T*ones(par.nu,1)];

%--define min and max parameter values
par.mmin=[ksipar.MCMC_DU_MINMAX(1)*ones(par.nu,1);ksipar.MCMC_T0_MINMAX(1)*ones(par.nu,1);ksipar.MCMC_T_MINMAX(1)*ones(par.nu,1)];
par.mmax=[ksipar.MCMC_DU_MINMAX(2)*ones(par.nu,1);ksipar.MCMC_T0_MINMAX(2)*ones(par.nu,1);ksipar.MCMC_T_MINMAX(2)*ones(par.nu,1)];

%--MCMC iterations-------------------------------------------%
for ksimgroup=1:nsimgroup
    %--initial model for the first run
    if ksimgroup==1
        if isempty(ksipar.METH_INIT_MCMC)==0
            minit(1:par.nu)=x1(1:par.nu)/u0;
            minit(par.nu+1:3*par.nu)=x1(par.nu+1:3*par.nu)/max(tobs);
            clear x1
        else
            minit(1:par.nu)=ksipar.INIT_MCMC_DU*ones(par.nu,1);
            minit(par.nu+1:2*par.nu)=ksipar.INIT_MCMC_T0*ones(par.nu,1);
            minit(2*par.nu+1:3*par.nu)=ksipar.INIT_MCMC_T*ones(par.nu,1);
        end
    end


    %--conduct mcmc exploration
    [chain,logP,alpha]=mcmc(minit,@(m)loglike(m,par),@(m)logmodelprior(m,par),dm,ksipar.NITER_MCMC,ksipar.MCMC_CHAIN_SKIP);

    %--initial model for next simulation
    minit=chain(end,:);

    %--make results with physical dimensions for recording
    chain(:,1:par.nu)=u0*chain(:,1:par.nu);
    chain(:,par.nu+1:2*par.nu)=max(tobs)*chain(:,par.nu+1:2*par.nu);
    chain(:,2*par.nu+1:3*par.nu)=max(tobs)*chain(:,2*par.nu+1:3*par.nu);

    %--save results (models chain, rejection rate alpha, -loglikelihood (logP)
    filesave=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/mcmc_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_nit',num2str(ksimgroup*(ksipar.NITER_MCMC/1e6)),'e6.mat'];
    save(filesave,'chain','alpha','logP');

end


function ll=loglike(m,par)

%--initialize loglikelihood and strain
ll=0.0;
strain=zeros(par.nj,par.nt);

%--loop on times
for i=1:par.nt

    %--compute slip
    psi=slip_basis_fn_sc_smooth(par.tobs(i),m(par.nu+1:2*par.nu),m(2*par.nu+1:3*par.nu));
    u=m(1:par.nu)'.*psi;

    %--strain
    strain(:,i)=par.gf_ezz*u;

    %--objective function
    ll=ll-0.5*((strain(:,i)-par.dobs(:,i))')*((strain(:,i)-par.dobs(:,i)).*par.cov_strain)-0.5*par.cov_slip*((par.MV*u-par.dobsms(i)).^2);


end


end

function lmp=logmodelprior(m,par)

%--uniform
%lmp=0;

%--uniform between min and max
lmp=0;
for i=1:length(m)
    if ((m(i)<par.mmin(i)) | (m(i)>par.mmax(i)))
        lmp=-10;
    end
end

end




