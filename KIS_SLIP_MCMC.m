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
% KIS_SLIP_MCMC: compute model predictions:
%              reconstruct mean slip/mean strain histories
%              and slip/strain history range at each node from
%              mcmc chain
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


%--load mcmc exploration results and compute model predictions
nsimgroup=floor(ksipar.NITER_MCMC_TOTAL/ksipar.NITER_MCMC);
ksimgroup0=floor(ksipar.NITER_MCMC_BURNIN/ksipar.NITER_MCMC)+1;

iter=0;
for k=ksimgroup0:nsimgroup

    fileparam=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/mcmc_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_nit',num2str(k*(ksipar.NITER_MCMC/1e6)),'e6.mat']
    load(fileparam)

    for i=1:ksipar.NITER_MCMC/ksipar.MCMC_CHAIN_SKIP
        if mod(i,ksipar.MCMC_FRECU)==0
            iter=iter+1;

            xx=squeeze(chain(i,:));
            strain=zeros(par.nj,par.ntt);
            strains=zeros(par.nj,par.ntt);
            meanslip=zeros(par.ntt,1);
            meanslips=zeros(par.ntt,1);
            matu=zeros(par.nu,par.ntt);
            matuu=zeros(par.nu,par.ntt);

            for j=1:par.ntt
                %--compute slip
                psi=slip_basis_fn_sc_smooth(ttobs(j),xx(par.nu+1:2*par.nu),xx(2*par.nu+1:3*par.nu));
                u=xx(1:par.nu)'.*psi;
                matu(:,j)=u;
                matuu(:,j)=u.^2;

                %--compute strain
                strain(:,j)=gf_ezz*u;
                strains(:,j)=strain(:,j).^2;

                %--compute meanslip
                meanslip(j)=MV*u;
                meanslips(j)=meanslips(j).^2;
            end

            if iter==1 
                strain_m=strain;
                strain_p=strain;
                strain_mean=strain;
                strains_mean=strains;
                meanslip_m=meanslip';
                meanslip_p=meanslip';
                meanslip_mean=meanslip';
                meanslips_mean=meanslips';
                matu_m=matu;
                matu_p=matu;
                matu_mean=matu;
                matuu_mean=matuu;
            else
                for j=1:par.nj
                    strain_m(j,:)=min([strain_m(j,:);strain(j,:)]);
                    strain_p(j,:)=max([strain_p(j,:);strain(j,:)]);
                    meanslip_m=min([meanslip_m;meanslip']);
                    meanslip_p=max([meanslip_p;meanslip']);
                end
                for j=1:par.nu
                    matu_m(j,:)=min([matu_m(j,:) matu(j,:)]);
                    matu_p(j,:)=max([matu_p(j,:) matu(j,:)]);
                end
                matu_mean=matu_mean+matu;
                matuu_mean=matuu_mean+matuu;
                strain_mean=strain_mean+strain;
                strains_mean=strains_mean+strains;
                meanslip_mean=meanslip_mean+meanslip;
                meanslips_mean=meanslips_mean+meanslips;
 
            end

        end
    end

end

matu_mean=matu_mean/iter;
matuu_mean=matuu_mean/iter;
strain_mean=strain_mean/iter;
strains_mean=strains_mean/iter;
meanslip_mean=meanslip_mean/iter;
meanslips_mean=meanslips_mean/iter;

matu_std=sqrt(matuu_mean-matu_mean.^2);
strain_std=sqrt(strains_mean-strain_mean.^2);
meanslip_std=sqrt(meanslips_mean-meanslip_mean.^2);

matu_mm=matu_mean-matu_std;ind=find(matu_mm<0);matu_mm(ind)=0;
matu_pp=matu_mean+matu_std;
strain_mm=strain_mean-strain_std;
strain_pp=strain_mean+strain_std;
meanslip_mm=meanslip_mean-meanslip_std;
meanslip_pp=meanslip_mean+meanslip_std;

ind0=find(matu_m<0);matu_m(ind0)=0;

%--interpolate slip--%
mu_mean=matu_mean;
mu_m=matu_m;
mu_p=matu_p;
mu_mm=matu_mm;
mu_pp=matu_pp;
matu_mean=zeros(par.nu,par.nt);
matu_m=zeros(par.nu,par.nt);
matu_p=zeros(par.nu,par.nt);
matu_mm=zeros(par.nu,par.nt);
matu_pp=zeros(par.nu,par.nt);
for j=1:par.nu
    matu_mean(j,:)=interp1(ttobs,squeeze(mu_mean(j,:)),tobs);
    matu_m(j,:)=interp1(ttobs,squeeze(mu_m(j,:)),tobs);
    matu_p(j,:)=interp1(ttobs,squeeze(mu_p(j,:)),tobs);
    matu_mm(j,:)=interp1(ttobs,squeeze(mu_mm(j,:)),tobs);
    matu_pp(j,:)=interp1(ttobs,squeeze(mu_pp(j,:)),tobs);
end

%--interpolate strains--%
sm=strain_m;
sp=strain_p;
smm=strain_mm;
spp=strain_pp;
strain_m=zeros(par.nj,par.nt);
strain_p=zeros(par.nj,par.nt);
strain_mm=zeros(par.nj,par.nt);
strain_pp=zeros(par.nj,par.nt);
for j=1:par.nj
    strain_m(j,:)=interp1(ttobs,squeeze(sm(j,:)),tobs);
    strain_p(j,:)=interp1(ttobs,squeeze(sp(j,:)),tobs);
    strain_mm(j,:)=interp1(ttobs,squeeze(smm(j,:)),tobs);
    strain_pp(j,:)=interp1(ttobs,squeeze(spp(j,:)),tobs);
end

%-interpolate average slip
msm=meanslip_m;
msp=meanslip_p;
msmm=meanslip_mm;
mspp=meanslip_pp;
meanslip_m=zeros(par.nt,1);
meanslip_p=zeros(par.nt,1);
meanslip_mm=zeros(par.nt,1);
meanslip_pp=zeros(par.nt,1);
meanslip_m=interp1(ttobs,msm,tobs);
meanslip_p=interp1(ttobs,msp,tobs);
meanslip_mm=interp1(ttobs,msmm,tobs);
meanslip_pp=interp1(ttobs,mspp,tobs);

%--save predictions
filesave=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/predictions_mcmc_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'.mat'];
save(filesave,'meanslip_m','meanslip_p','meanslip_mm','meanslip_pp',...
              'strain_m','strain_p','strain_mm','strain_pp',...
              'matu_m','matu_p','matu_mm','matu_pp','matu_mean');

