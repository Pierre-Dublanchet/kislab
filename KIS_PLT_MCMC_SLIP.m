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
% KIS_PLT_MCMC_SLIP: plot mcmc (Bayesian) inversion results
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load deterministic inversion result
file_invd=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/inv_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_ls',num2str(ksipar.REGUL),'.mat'];
load(file_invd)

%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);

par.nt=length(tobs);
data.xdata=tobs;
par.s0=max(abs(dobs'));   %--strain scale
par.u0=max(abs(dobsms));  %--slip scale
par.nj=length(par.s0);
for j=1:par.nj
    %data.ydata(j,:)=dobs(j,:)/par.s0(j); %--normalized strain observations
    data.ydata(j,:)=dobs(j,:)/max(par.s0); %--normalized strain observations
end
data.ydata(par.nj+1,:)=dobsms/par.u0; %--normalized mean slip

%--data error stdev
std_dobs=ksipar.STD_STRAIN*ones(1,par.nj);
std_dobsms=ksipar.STD_SLIP;
std_dobs=std_dobs./sqrt(ksipar.GAUGE_CONF_LEV);

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,par.nu]=size(gf_ezz);
for j=1:par.nj
    %par.gf_ezz(j,:)=gf_ezz(j,:)*par.u0/par.s0(j); %--normalized stiffness
    par.gf_ezz(j,:)=gf_ezz(j,:)*par.u0/max(par.s0); %--normalized stiffness
end

%--load adjacency matrix
file_adja=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat'];
load(file_adja)

%--average  operator matrix--%
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_mv)
par.MV=MV;

%--fault nodes  coordinates--%
gcfault0=g_coord_fault;
xf0=gcfault0(:,1);
yf0=gcfault0(:,2);

%--strain gages coordinates--%
gcjauge=[x_jauge y_jauge z_jauge];
gcjauge(:,1)=gcjauge(:,1)/cosd(90-ksipar.TH);
gcjauge(:,3)=0;

%--Model resolution matrix
file_resolution=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/Resolution_matrix_r',num2str(ksipar.HMAXF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(file_resolution)

%--load range of model predictions
file_pred=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/predictions_mcmc_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'.mat'];
load(file_pred)

%--compute slip rate from mean reconstructed slip history
for i=1:par.nu
    matv_mean(i,:)=gradient(matu_mean(i,:),tobs);
end

%--calculate maximum slip and max slip rate for colorscale
uaxmax=max(max(matu_mean));
vaxmax=max(max(matv_mean));

%--make a regular grid to interpolate resolution and make a nice contour plot
dx=2*ksipar.R/ksipar.NPTS_RES_MAP;
lx=[dx:dx:2*ksipar.R];vlx=sort([-lx 0 lx]);
ly=lx/2;vly=sort([-ly 0 ly]);
[xq,yq]=meshgrid(vlx,vly);
xq=reshape(xq,length(vlx)*length(vly),1);
yq=reshape(yq,length(vlx)*length(vly),1);

%--initialize figures
figure(1);
tiledlayout(ksipar.NXFIG_UMAP,ksipar.NYFIG_UMAP,"TileSpacing","tight");
figure(2);
tiledlayout(ksipar.NXFIG_UMAP,ksipar.NYFIG_UMAP,"TileSpacing","tight");
figure(3);
tiledlayout(ksipar.NXFIG_UMAP,ksipar.NYFIG_UMAP,"TileSpacing","tight");
figure(4);
tiledlayout(ksipar.NXFIG,ksipar.NYFIG,'TileSpacing','tight');


kfig=0;
matu=zeros(par.nu,par.nt);
for i=1:par.nt
    %--result deterministic inversion
    psi=slip_basis_fn_sc_smooth(data.xdata(i),x1(par.nu+1:2*par.nu),x1(2*par.nu+1:3*par.nu));
    u_d=x1(1:par.nu).*psi;
    u_d=u_d/par.u0;

    %--result mean model mcmc
    u=squeeze(matu_mean(:,i));
    u=u/par.u0;

    %--slip rate
    slip_rate=matv_mean(:,i);

    %--strain
    strain(:,i)=par.gf_ezz*u;
    straind(:,i)=par.gf_ezz*u_d;

    %--meanslip
    meanslip(i)=par.MV*u;
    meanslipd(i)=par.MV*u_d;

    if mod(i,floor(par.nt/ksipar.NFIG_UMAP))==0
        kfig=kfig+1;

        xfig=floor((kfig-1)/ksipar.NYFIG_UMAP)+1;
        yfig=kfig-(xfig-1)*ksipar.NYFIG_UMAP;

        %--plot slip maps
        figure(1);
        nexttile(kfig)

        patch('Faces',g_num_fault,'Vertices',100*[xf0 yf0],'FaceVertexCData', [u]*par.u0*1e3, 'FaceColor', 'interp',...
            'FaceVertexAlphaData',vec_cache_lowres,'AlphaDataMapping','scaled','FaceAlpha','interp');
        alphamap([0.7 1]);

        caxis([0 uaxmax*1e3]);

        hold on
        plot3(100*gcjauge(:,1),100*gcjauge(:,2),max(u)*1e3*ones(par.nj,1),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
        if kfig==1
            txg=text(100*gcjauge(:,1)-1,100*gcjauge(:,2)-0.4,ksipar.STRG);
            set(txg,'Color',[0.7 0 0],'Fontsize',13)
        end



        xlim([-2*ksipar.R 2*ksipar.R]*100);ylim([-ksipar.R ksipar.R]*100);
        if xfig==ksipar.NXFIG_UMAP
            xlabel('$x_1$ (cm)','Interpreter','Latex');
        else
            set(gca,'XTickLabel',[])
        end
        if yfig==1
            ylabel('$x_2$ (cm)','Interpreter','Latex');
        else
            set(gca,'YTickLabel',[])
        end

        if ((xfig==ksipar.NXFIG_UMAP) & (yfig==ksipar.NYFIG_UMAP))
            cb=colorbar;
            ylabel(cb,'slip ($\mu$m)','Interpreter','Latex')
        end

        title([num2str(floor(tobs(i)*100)/100),' s'],'Interpreter','Latex')
        set(gca,'Fontsize',15)


        %--plot slip uncertainty maps
        figure(2);
        nexttile(kfig)

        patch('Faces',g_num_fault,'Vertices',100*[xf0 yf0],'FaceVertexCData', 0.5*[matu_p(:,i)-matu_m(:,i)]*1e3, 'FaceColor', 'interp',...
            'FaceVertexAlphaData',vec_cache_lowres,'AlphaDataMapping','scaled','FaceAlpha','interp');
        alphamap([0.7 1]);

        hold on
        plot3(100*gcjauge(:,1),100*gcjauge(:,2),max(u)*1e3*ones(par.nj,1),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
        if kfig==1
            txg=text(100*gcjauge(:,1)-1,100*gcjauge(:,2)-0.4,ksipar.STRG);
            set(txg,'Color',[0.7 0 0],'Fontsize',13)
        end

        xlim([-2*ksipar.R 2*ksipar.R]*100);ylim([-ksipar.R ksipar.R]*100);
        if xfig==ksipar.NXFIG_UMAP
            xlabel('$x_1$ (cm)','Interpreter','Latex');
        else
            set(gca,'XTickLabel',[])
        end
        if yfig==1
            ylabel('$x_2$ (cm)','Interpreter','Latex');
        else
            set(gca,'YTickLabel',[])
        end

        caxis([0 0.5*uaxmax*1e3]);

        if ((xfig==ksipar.NXFIG_UMAP) & (yfig==ksipar.NYFIG_UMAP))
            cb=colorbar;
            ylabel(cb,'std $\sigma_{\delta}$ ($\mu$m)','Interpreter','Latex')
        end

        title([num2str(floor(tobs(i)*100)/100),' s'],'Interpreter','Latex')
        set(gca,'Fontsize',15)

        %--plot slip rate maps
        figure(3);
        nexttile(kfig)

        patch('Faces',g_num_fault,'Vertices',100*[xf0 yf0],'FaceVertexCData', [slip_rate]*1e3, 'FaceColor', 'interp',...
            'FaceVertexAlphaData',vec_cache_lowres,'AlphaDataMapping','scaled','FaceAlpha','interp');
        alphamap([0.7 1]);

        caxis([0 vaxmax*1e3]);

        hold on
        plot3(100*gcjauge(:,1),100*gcjauge(:,2),max(slip_rate)*1e3*ones(par.nj,1),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
        if kfig==1
            txg=text(100*gcjauge(:,1)-1,100*gcjauge(:,2)-0.4,ksipar.STRG);
            set(txg,'Color',[0.7 0 0],'Fontsize',13)
        end



        xlim([-2*ksipar.R 2*ksipar.R]*100);ylim([-ksipar.R ksipar.R]*100);
        if xfig==ksipar.NXFIG_UMAP
            xlabel('$x_1$ (cm)','Interpreter','Latex');
        else
            set(gca,'XTickLabel',[])
        end
        if yfig==1
            ylabel('$x_2$ (cm)','Interpreter','Latex');
        else
            set(gca,'YTickLabel',[])
        end

        if ((xfig==ksipar.NXFIG_UMAP) & (yfig==ksipar.NYFIG_UMAP))
            cb=colorbar;
            ylabel(cb,'slip rate ($\mu$m.s$^{-1}$)','Interpreter','Latex')
        end

        title([num2str(floor(tobs(i)*100)/100),' s'],'Interpreter','Latex')
        set(gca,'Fontsize',15)


    end


end

%--plot observations and models predictions
%--strain
for i=1:par.nj
    s0min(i)=min([data.ydata(i,:) strain(i,:) strain_p(i,:) strain_m(i,:)]*max(par.s0));
    s0max(i)=max([data.ydata(i,:) strain(i,:) strain_p(i,:) strain_m(i,:)]*max(par.s0));

    figure(4);
    nexttile(i)
    %--data and data uncertainty
    plot(tobs,data.ydata(i,:)*max(par.s0),'-k','Linewidth',2)
    hold on
    plot(tobs,data.ydata(i,:)*max(par.s0)-std_dobs(i),'-','Color',[1 1 1],'Linewidth',2)
    hold on
    plot(tobs,data.ydata(i,:)*max(par.s0)+std_dobs(i),'-','Color',[1 1 1],'Linewidth',2)

    ftest=fill([tobs' fliplr(tobs')], [data.ydata(i,:)*max(par.s0)+std_dobs(i) fliplr(data.ydata(i,:)*max(par.s0)-std_dobs(i))],[0.7 0.7 0.7],'FaceAlpha',0.5,'EdgeColor',[1 1 1]);

    %--models predictions
    ttobs=[0:max(tobs)/19:max(tobs)];
    plot(tobs,strain(i,:)*max(par.s0),'-','Color',[0.7 0 0],'Linewidth',2)
    hold on
    plot(tobs,strain_pp(i,:),'--','Color',[0.7 0 0],'Linewidth',2)
    hold on
    plot(tobs,strain_mm(i,:),':','Color',[0.7 0 0],'Linewidth',2)
    hold on
    plot(tobs,straind(i,:)*max(par.s0),'-b','Linewidth',2)
    txg1(i)=text(80*max(tobs)/100,s0min(i)-std_dobs(i)+0.9*(s0max(i)-s0min(i)+2*std_dobs(i)),['G',num2str(i)]);
    xlim([0 max(tobs)]);ylim([1.1*(s0min(i)-std_dobs(i)) 1.1*(s0max(i)+std_dobs(i))]);
    if i>ksipar.NXFIG*(ksipar.NYFIG-1)
        xlabel('time (s)','Interpreter','Latex');
    else
        set(gca,'XTicklabel',[]);
    end
    if mod(i,ksipar.NXFIG)==1
        ylabel('$\varepsilon_{11}$','Interpreter','Latex');
    end
    set(gca,'Fontsize',15)


end
set(txg1,'Interpreter','latex','Fontsize',15)


u0min=min([data.ydata(par.nj+1,:) meanslip meanslip_p' meanslip_m']*par.u0*1e3);
u0max=max([data.ydata(par.nj+1,:) meanslip meanslip_p' meanslip_m']*par.u0*1e3);

%--average slip
figure(4)
nexttile(par.nj+1)
pltu(1)=plot(tobs,data.ydata(par.nj+1,:)*par.u0*1e3,'-k','Linewidth',2)
hold on
plot(tobs,data.ydata(par.nj+1,:)*par.u0*1e3-std_dobsms*1e3,'-','Color',[0.7 0.7 0.7],'Linewidth',2)
hold on
plot(tobs,data.ydata(par.nj+1,:)*par.u0*1e3+std_dobsms*1e3,'-','Color',[0.7 0.7 0.7],'Linewidth',2)
fill([tobs' fliplr(tobs')], [data.ydata(par.nj+1,:)*par.u0+std_dobsms fliplr(data.ydata(par.nj+1,:)*par.u0-std_dobsms)]*1e3,[0.7 0.7 0.7],'FaceAlpha',0.5,'EdgeColor',[1 1 1]);
hold on
pltu(2)=plot(tobs,meanslip*par.u0*1e3,'-r','Color',[0.7 0 0],'Linewidth',2)
hold on
pltu(3)=plot(tobs,meanslip_p*1e3,'--r','Color',[0.7 0 0],'Linewidth',2)
hold on
pltu(4)=plot(tobs,meanslip_m*1e3,':r','Color',[0.7 0 0],'Linewidth',2)
hold on
pltu(5)=plot(tobs,meanslipd*par.u0*1e3,'-b','Linewidth',2)
xlim([0 max(tobs)]);ylim([u0min-std_dobsms*1e3 u0min-std_dobsms*1e3+1.5*(u0max-u0min+2*std_dobsms*1e3)]);
xlabel('time (s)','Interpreter','Latex');
ylabel('$\delta_m$ ($\mu m$)','Interpreter','Latex');
leg=legend(pltu,'data','mean $\bar{\delta}$','$\bar{\delta}+\sigma_{\delta}$','$\bar{\delta}-\sigma_{\delta}$','deterministic');
set(leg,'Interpreter','latex','location','NorthWest')
set(gca,'Fontsize',15)




%--Calculate rupture speed

xi=xf0(ksipar.FIRSTACTIVATEDNODE);yi=yf0(ksipar.FIRSTACTIVATEDNODE);
df=sqrt((xf0-xi).^2+(yf0-yi).^2);
[df,indf]=sort(df);

matu=matu_mean;



for i=1:par.nu
    if isempty(find(matu(i,:)>ksipar.UTHRESH_RUPSPEED))==0
        tu(i)=tobs(min(find(matu(i,:)>ksipar.UTHRESH_RUPSPEED)));
    else
        tu(i)=0;
    end
    if isempty(find(matu_p(i,:)>ksipar.UTHRESH_RUPSPEED))==0
        tum(i)=tobs(min(find(matu_pp(i,:)>ksipar.UTHRESH_RUPSPEED)));
    else
        tum(i)=max(tobs);
    end
    if isempty(find(matu_m(i,:)>ksipar.UTHRESH_RUPSPEED))==0
        tup(i)=tobs(min(find(matu_mm(i,:)>ksipar.UTHRESH_RUPSPEED)));
    else
        tup(i)=max(tobs);
    end
    dumaxend(i)=max(matu(i,:));
end
vtu=tu;
indz=find(vtu==0);
vtu(indz)=1.1*max(tobs);
tu=tu(indf);tum=tum(indf);tup=tup(indf);
dumaxend=dumaxend(indf);
indpostu=find(tu>0);

mTU = griddata(xf0,yf0,vtu,xq,yq);
mTU=reshape(mTU,length(vlx),length(vly));


figure(5);clf;
subplot(1,2,2,'align')
scatter(df(indpostu)*100,tu(indpostu),50,dumaxend(indpostu)*1e3,'filled')
cb=colorbar;
caxis([0 uaxmax]*1e3);
hold on
errorbar(df(indpostu)*100,tu(indpostu),tu(indpostu)-tum(indpostu),tup(indpostu)-tu(indpostu),'ok');
hold on
plot(df*100,min(tu(indpostu))+df*100./ksipar.VR_REF(1),'--','Color',[.7 0 0],'Linewidth',2)
hold on
plot(df*100,min(tu(indpostu))+df*100./ksipar.VR_REF(2),'--','Color',[.7 0 0],'Linewidth',2)
hold on
plot(df*100,min(tu(indpostu))+df*100./ksipar.VR_REF(3),'--','Color',[.7 0 0],'Linewidth',2)


% tx(1)=text(2.2,36,['100 m.day$^{-1}$']);
% tx(2)=text(4.8,39,['200 m.day$^{-1}$']);
% tx(3)=text(5,27,['500 m.day$^{-1}$']);
% tx(4)=text(0.2,47,'(b)');


% set(tx(1),'Interpreter','latex','Fontsize',15,'Color',[.7 0 0],'Rotation',58)
% set(tx(2),'Interpreter','latex','Fontsize',15,'Color',[.7 0 0],'Rotation',35)
% set(tx(3),'Interpreter','latex','Fontsize',15,'Color',[.7 0 0],'Rotation',15)
% set(tx(4),'Fontsize',20,'Interpreter','latex');

ylabel(cb,'max slip $\delta(t_{max})$ ($\mu$m)','Interpreter','Latex');
ylabel('Time $t_{r}$ (s)','Interpreter','Latex');
xlabel('Distance from maximum of slip (cm)','Interpreter','Latex');
set(gca,'Fontsize',20)


figure(5);
ax(1)=subplot(1,2,1,'align')
patch('Faces',g_num_fault,'Vertices',100*[xf0 yf0],'FaceVertexCData', zeros(par.nu,1), 'FaceColor', [1 1 1]);
cb=colorbar
hold on
contour(vlx*100,vly*100,mTU,[0:max(tobs)/10:max(tobs)],'-','linewidth',2,'ShowText','off')
colormap(ax(1),"jet")
hold on
plot3(100*gcjauge(:,1),100*gcjauge(:,2),max(vtu)*1e3*ones(par.nj,1),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
hold on
plot3(100*xi,100*yi,max(vtu)*1e3,'p','MarkerFaceColor',[0.7 0 0.7],'MarkerEdgeColor',[0.7 0 0.7],'Markersize',15)
tx(5)=text(-1.9*ksipar.R*100,0.9*ksipar.R*100,'(a)');
set(tx(5),'Fontsize',20,'Interpreter','latex');
xlabel('$x_1$ (cm)','Interpreter','Latex');
ylabel('$x_2$ (cm)','Interpreter','Latex');
ylabel(cb,'Time $t_{2.0}$ (s)','Interpreter','Latex');
set(gca,'Fontsize',20)



%--compute and plot standard deviation on slip statistics
[h1,edges1]=histcounts(0.5*1e3*reshape(matu_pp-matu_mm,par.nu*length(tobs),1),50,'Normalization','cdf');
[h2,edges2]=histcounts(0.5*1e3*reshape(matu_pp(find(vec_cache_lowres<0.5),:)-matu_mm(find(vec_cache_lowres<0.5),:),length(find(vec_cache_lowres<0.5))*length(tobs),1),50,'Normalization','cdf');
[h3,edges3]=histcounts(0.5*1e3*reshape(matu_pp(find(vec_cache_lowres>0.5),:)-matu_mm(find(vec_cache_lowres>0.5),:),length(find(vec_cache_lowres>0.5))*length(tobs),1),50,'Normalization','cdf');


figure(6);clf;
plot(0.5*(edges1(1:end-1)+edges1(2:end)),h1,'-k','Linewidth',2)
hold on;plot(0.5*(edges2(1:end-1)+edges2(2:end)),h2,'-b','Linewidth',2)
hold on;plot(0.5*(edges3(1:end-1)+edges3(2:end)),h3,'-r','Linewidth',2)
xlabel('standard deviation of slip $\sigma_{\delta}$ ($\mu$m)','Interpreter','Latex');ylabel('cdf','Interpreter','Latex');
leg=legend('all nodes','resolution $<0.05$','resolution $>0.05$')
set(leg,'Interpreter','latex','Location','SouthEast');
set(gca,'Fontsize',15)

%--rescale and print figures
figure(1)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2*vop(3) 1.5*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/slip_map_meanmcmc_',ksipar.INV_FILE_ID,'.png'];
print(1,'-dpng',filefig);

figure(2)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2*vop(3) 1.5*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/sigslip_map_meanmcmc_',ksipar.INV_FILE_ID,'.png'];
print(2,'-dpng',filefig);

figure(3)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2*vop(3) 1.5*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/slip_rate_map_meanmcmc_',ksipar.INV_FILE_ID,'.png'];
print(3,'-dpng',filefig);

figure(4)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 1.3*vop(3) 1.6*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/fit_jauge_meanmcmc_',ksipar.INV_FILE_ID,'.png'];
print(4,'-dpng',filefig);

figure(5)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2.8*vop(3) 1.3*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/vr_mcmc_',ksipar.INV_FILE_ID,'.png'];
print(5,'-dpng',filefig);

figure(6)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2*vop(3) 1.5*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/cdf_std_slip_',ksipar.INV_FILE_ID,'.png'];
print(6,'-dpng',filefig);



