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
% KIS_RMS_DET_MCMC: plot rms of deterministic inversion
%                   and mcmc chain.
%                   plot acceptance rate of mcmc iteration
%                   
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat


%--total number of mcmc iterations
nsimu=floor(ksipar.NITER_MCMC_TOTAL/ksipar.NITER_MCMC);
isimu0=floor(ksipar.NITER_MCMC_BURNIN/ksipar.NITER_MCMC)+1;
isimu0=1;

rms=zeros((nsimu-isimu0+1)*ksipar.NITER_MCMC/ksipar.MCMC_CHAIN_SKIP,1);
valpha=zeros(nsimu-isimu0+1,2);


for isimu=isimu0:nsimu   %--loop on simulations

    %--load mcmc chain
    fileparam=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/mcmc_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_nit',num2str(isimu*(ksipar.NITER_MCMC/1e6)),'e6.mat']
    load(fileparam)

    rms((isimu-isimu0)*ksipar.NITER_MCMC/ksipar.MCMC_CHAIN_SKIP+1:(isimu-isimu0+1)*ksipar.NITER_MCMC/ksipar.MCMC_CHAIN_SKIP)=sqrt(-2*logP(:,2)/(ksipar.NGOK*ksipar.MCMC_TIME_SAMP));
    valpha(isimu-isimu0+1,:)=alpha;

    if isimu==1
        rmsdet=sqrt(-2*logP(1,2)/(ksipar.NGOK*ksipar.MCMC_TIME_SAMP));
    end
end

ind=find(rms>0);


[pdf_rms,rms_edges] = histcounts(rms(ind),'Normalization','pdf');
rms_c=0.5*(rms_edges(1:end-1)+rms_edges(2:end));

figure(1);clf;
subplot(2,3,[1 2],'align')
semilogx([1:length(rms(ind))]*ksipar.MCMC_CHAIN_SKIP,rms(ind),'-k')
hold on
semilogx(1,rmsdet,'p','MarkerFaceColor',[.7 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',15)
tx(1)=text(1.5,0.63,'(a)')
ylabel('$RMS_{i}=\sqrt{2J/N_g N_t}$','Interpreter','latex')
xlabel('number of MCMC iterations','Interpreter','latex')
set(gca,'Fontsize',15)

subplot(2,3,3,'align')
plot(pdf_rms,rms_c,'-k','Linewidth',1)
hold on;
plot(0,rmsdet,'p','MarkerFaceColor',[.7 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',15)
tx(2)=text(1,0.63,'(b)')
xlabel('pdf','Interpreter','latex')
leg=legend('MCMC chain','Deterministic')
set(leg,'Interpreter','latex')
set(tx,'Fontsize',20,'Interpreter','latex')
set(gca,'Fontsize',15)

subplot(2,3,[4 5],'align')
plot([isimu0:nsimu]*ksipar.NITER_MCMC,1-valpha(:,1),'-ok')
hold on
plot([isimu0:nsimu]*ksipar.NITER_MCMC,1-valpha(:,2),'--ok')
ylabel('acceptance rate $\alpha$','Interpreter','latex')
xlabel('number of MCMC iterations','Interpreter','latex')
set(gca,'Fontsize',15)

figure(1)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.2*vop(2) 1.8*vop(3) 1.2*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/rmsidetmcmc_',ksipar.INV_FILE_ID,'.png'];
print(1,'-dpng',filefig)

