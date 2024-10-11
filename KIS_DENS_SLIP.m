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
% KIS_DENS_SLIP: compute pdf of reconstructed slip history 
%                for each fault node
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);
par.nt=length(tobs);
par.s0=max(abs(dobs'));   %--stress scale (Pa)
par.u0=max(abs(dobsms));  %--slip scale (mm)

%--undersample data------%
[par.nj,nsdraw]=size(dobs);
ttobs=[0:max(tobs)/(ksipar.MCMC_TIME_SAMP-1):max(tobs)];
par.ntt=length(ttobs);

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,par.nu]=size(gf_ezz);
 
%--average  operator matrix--%
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_mv)



%--load mcmc exploration result and compute model predictions
nsimgroup=floor(ksipar.NITER_MCMC_TOTAL/ksipar.NITER_MCMC);
ksimgroup0=floor(ksipar.NITER_MCMC_BURNIN/ksipar.NITER_MCMC)+1;

iter=0;

n_model_considered=(ksipar.NITER_MCMC_TOTAL-ksipar.NITER_MCMC_BURNIN)/ksipar.MCMC_FRECU;


    for j=1:par.ntt
        mslipdens(j).v=zeros(par.nu,n_model_considered);
    end

for k=ksimgroup0:nsimgroup

    fileparam=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/mcmc_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_nit',num2str(k*(ksipar.NITER_MCMC/1e6)),'e6.mat']
    load(fileparam)

    for i=1:ksipar.NITER_MCMC/ksipar.MCMC_CHAIN_SKIP
        if mod(i,ksipar.MCMC_FRECU)==0
            iter=iter+1;
            xx=squeeze(chain(i,:));
            for j=1:par.ntt
                %--compute slip
                psi=slip_basis_fn_sc_smooth(ttobs(j),xx(par.nu+1:2*par.nu),xx(2*par.nu+1:3*par.nu));
                u=xx(1:par.nu)'.*psi;
                mslipdens(j).v(:,iter)=u;
            end
        end
    end
end

vslipbin=[0:ksipar.MCMC_DU_MINMAX(2)*par.u0/ksipar.NSBIN:ksipar.MCMC_DU_MINMAX(2)*par.u0];

for i=1:par.nu
    matslipdens(i).mm=zeros(ksipar.NSBIN,par.ntt);
    for j=1:par.ntt
        matslipdens(i).mm(:,j)=histcounts(squeeze(mslipdens(j).v(i,:)),vslipbin);
    end
    indz=find(matslipdens(i).mm==0);
    matslipdens(i).mm(indz)=NaN;
end

nnodes=par.nu;

%--save pdf (slip density maps)
filesave=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/slip_dens_map_',ksipar.INV_FILE_ID,'.mat']
save(filesave,'matslipdens','vslipbin','nnodes','ttobs');


