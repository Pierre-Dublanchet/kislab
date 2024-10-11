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
% KIS_PLT_PDF: plot joint posterior pdf of parameters
%              from mcmc (Bayesian) exploration
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--Model resolution matrix
file_resolution=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/Resolution_matrix_r',num2str(ksipar.HMAXF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(file_resolution)

%--load adjacency matrix
file_adja=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat'];
load(file_adja)

%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);

par.s0=max(abs(dobs'));   %--strain scale
par.u0=max(abs(dobsms));  %--slip scale
par.nj=length(par.s0);

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,par.nu]=size(gf_ezz);
for j=1:par.nj
    %par.gf_ezz(j,:)=gf_ezz(j,:)*par.u0/par.s0(j); %--normalized stiffness
    par.gf_ezz(j,:)=gf_ezz(j,:)*par.u0/max(par.s0); %--normalized stiffness
end

%--fault nodes  coordinates--%
gcfault0=g_coord_fault;
xf0=gcfault0(:,1);
yf0=gcfault0(:,2);

%--strain gages coordinates--%
gcjauge=[x_jauge y_jauge z_jauge];
gcjauge(:,1)=gcjauge(:,1)/cosd(90-ksipar.TH);
gcjauge(:,3)=0;

nn=length(DR);    %--number of nodes



%--define window and grid to compute and draw pdf
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);
u0=max(abs(dobsms));
tobsmax=max(tobs);
dumax=ksipar.MCMC_DU_MINMAX(2)*u0*1e3;                     %--in microns
t0max=ksipar.MCMC_T0_MINMAX(2)*tobsmax;
Tmax=ksipar.MCMC_T_MINMAX(2)*tobsmax;

dumin=ksipar.MCMC_DU_MINMAX(1)*u0*1e3;                     %--in microns
t0min=ksipar.MCMC_T0_MINMAX(1)*tobsmax;
Tmin=ksipar.MCMC_T_MINMAX(1)*tobsmax;



dt0=(t0max-t0min)/ksipar.NPTS_PDF;
dT=(Tmax-Tmin)/ksipar.NPTS_PDF;
ddu=(dumax-dumin)/ksipar.NPTS_PDF;


vt0=[t0min:dt0:t0max];
vT=[Tmin:dT:Tmax];
vdu=[dumin:ddu:dumax];

vct0=0.5*(vt0(1:end-1)+vt0(2:end));
vcT=0.5*(vT(1:end-1)+vT(2:end));
vcdu=0.5*(vdu(1:end-1)+vdu(2:end));

%--compute joint pdf
for i=1:nn
    pdf(i).Tdu=zeros(ksipar.NPTS_PDF,ksipar.NPTS_PDF);
    pdf(i).Tt0=zeros(ksipar.NPTS_PDF,ksipar.NPTS_PDF);
    pdf(i).t0du=zeros(ksipar.NPTS_PDF,ksipar.NPTS_PDF);
end


%--total number of mcmc iterations
nsimu=floor(ksipar.NITER_MCMC_TOTAL/ksipar.NITER_MCMC);
isimu0=floor(ksipar.NITER_MCMC_BURNIN/ksipar.NITER_MCMC)+1;

for isimu=isimu0:nsimu   %--loop on simulations

    %--load mcmc chain
    fileparam=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/mcmc_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_nit',num2str(isimu*(ksipar.NITER_MCMC/1e6)),'e6.mat']
    load(fileparam)

    for i=1:nn   %--loop on nodes


        %--extract parameter vectors from chain
        du=chain(:,i)*1e3;  %--slip in micron
        t0=chain(:,nn+i);
        T=chain(:,2*nn+i);

        %--compute joint pdf
        [N1,Xedges,Yedges] = histcounts2(T,du,vT,vdu,'Normalization','pdf');
        [N2,Xedges,Yedges] = histcounts2(T,t0,vT,vt0,'Normalization','pdf');
        [N3,Xedges,Yedges] = histcounts2(t0,du,vt0,vdu,'Normalization','pdf');

        pdf(i).Tdu=pdf(i).Tdu+N1*ksipar.NITER_MCMC/ksipar.NITER_MCMC_TOTAL;
        pdf(i).Tt0=pdf(i).Tt0+N2*ksipar.NITER_MCMC/ksipar.NITER_MCMC_TOTAL;
        pdf(i).t0du=pdf(i).t0du+N3*ksipar.NITER_MCMC/ksipar.NITER_MCMC_TOTAL;
    end

    clear chain

end

%--Identify non explored regions
for i=1:nn
    indz=find(pdf(i).Tdu==0);
    pdf(i).Tdu(indz)=NaN;
    indz=find(pdf(i).Tt0==0);
    pdf(i).Tt0(indz)=NaN;
    indz=find(pdf(i).t0du==0);
    pdf(i).t0du(indz)=NaN;
end



figure(1);clf;
figure(2);clf;
figure(3);clf;

figure(1);tiledlayout(ksipar.NXFIG_NODES,ksipar.NYFIG_NODES,"TileSpacing","none");
figure(2);tiledlayout(ksipar.NXFIG_NODES,ksipar.NYFIG_NODES,"TileSpacing","none");
figure(3);tiledlayout(ksipar.NXFIG_NODES,ksipar.NYFIG_NODES,"TileSpacing","none");
%--plot joint pdfs
for i=1:nn


    %--position of the panel in the figure from node number
    fpos(i)=i;
    xf=floor((fpos(i)-1)/ksipar.NYFIG_NODES)+1;
    yf=fpos(i)-(xf-1)*ksipar.NYFIG_NODES;


    figure(1);
    ax(1)=nexttile(fpos(i))
    truc=imagesc(vcdu,vcT,pdf(i).Tdu/max(max(pdf(i).Tdu)))
    txnode(1)=text(0.02*dumax,0.1*Tmax,['node ',num2str(i),'(',num2str(floor(xf0(i)*1000)/10),',',num2str(floor(yf0(i)*1000)/10),')'])
    set(truc,'AlphaData',1-isnan(pdf(i).Tdu))
    caxis([0 1]);
    if xf==ksipar.NXFIG_NODES
        xlabel('$\Delta u$ ($\mu$m)','Interpreter','latex');
    else
        set(gca,'XTickLabel',[])
    end
    if yf==1
        ylabel('T (s)','Interpreter','latex');
    else
        set(gca,'YTickLabel',[])
    end
    if (fpos(i)==nn)
        cb=colorbar;
        ylabel(cb,'pdf/max(pdf)','Interpreter','latex');
    end
    xlim([dumin dumax]);ylim([Tmin Tmax]);
    set(gca,'Fontsize',12)


    figure(2);
    ax(2)=nexttile(fpos(i))
    truc=imagesc(vct0,vcT,pdf(i).Tt0/max(max(pdf(i).Tt0)));
    txnode(2)=text(0.02*t0max,0.1*Tmax,['node ',num2str(i),'(',num2str(floor(xf0(i)*1000)/10),',',num2str(floor(yf0(i)*1000)/10),')'])
    set(truc,'AlphaData',1-isnan(pdf(i).Tt0))
    caxis([0 1]);
    
    if xf==ksipar.NXFIG_NODES
        xlabel('$t_0$ (s)','Interpreter','latex');
   else
        set(gca,'XTickLabel',[])
    end
    if yf==1
        ylabel('T (s)','Interpreter','latex');
    else
        set(gca,'YTickLabel',[])
    end
    if (fpos(i)==nn)
        cb=colorbar;
        ylabel(cb,'pdf/max(pdf)','Interpreter','latex');
    end
    xlim([t0min t0max]);ylim([Tmin Tmax]);
    set(gca,'Fontsize',12)

    figure(3);
    ax(3)=nexttile(fpos(i))
    truc=imagesc(vcdu,vct0,pdf(i).t0du/max(max(pdf(i).t0du)))
    txnode(3)=text(0.02*dumax,0.1*t0max,['node ',num2str(i),'(',num2str(floor(xf0(i)*1000)/10),',',num2str(floor(yf0(i)*1000)/10),')'])
    set(truc,'AlphaData',1-isnan(pdf(i).t0du))
    caxis([0 1]);
    if xf==ksipar.NXFIG_NODES
        xlabel('$\Delta u$ ($\mu$m)','Interpreter','latex');
    else
        set(gca,'XTickLabel',[])
    end
    if yf==1
        ylabel('$t_0$ (s)','Interpreter','latex');
    else
        set(gca,'YTickLabel',[])
    end
    if (fpos(i)==nn)
        cb=colorbar;
        ylabel(cb,'pdf/max(pdf)','Interpreter','latex');
    end
    xlim([dumin dumax]);ylim([t0min t0max]);
    set(gca,'Fontsize',12)

    for isp=1:3
        colormap(ax(isp),'autumn');
    end
    set(txnode,'Interpreter','latex','Fontsize',12)

end



for i=1:3
    figure(i)
    vop=get(gcf,'outerposition');
    set(gcf,'outerposition',[0.01*vop(1) 0.01*vop(2) 2*vop(3) 2*vop(4)]);
end


print(1,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/pdf_Tdu_mcmc_',ksipar.INV_FILE_ID,'.png'])
print(2,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/pdf_Tto_mcmc_',ksipar.INV_FILE_ID,'.png'])
print(3,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/pdf_todu_mcmc_',ksipar.INV_FILE_ID,'.png'])






