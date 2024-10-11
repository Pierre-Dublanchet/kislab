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
% KIS_PLT_DET_SLIP: plot deterministic inversion results
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load deterministic inversion result
file_invd=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/inv_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_ls',num2str(ksipar.REGUL),'.mat'];
load(file_invd)

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,nu]=size(gf_ezz);

%--average  operator matrix--%
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_mv)

%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);

%--calculate maximum slip for colorscale
psi=slip_basis_fn_sc_smooth(max(tobs),x1(par.nu+1:2*par.nu),x1(2*par.nu+1:3*par.nu));
uaxmax=max(psi.*x1(1:par.nu));

%--load adjacency matrix
file_adja=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat'];
load(file_adja)

%--Model resolution matrix
file_resolution=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/Resolution_matrix_r',num2str(ksipar.HMAXF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(file_resolution)


%--initial model
tref=max(tobs);
uref=max(abs(dobsms'));
x0=zeros(3*par.nu,1);
x0(1:par.nu)=ksipar.INIT_DET_DU*uref*ones(par.nu,1);
x0(par.nu+1:2*par.nu)=ksipar.INIT_DET_T0*tref*ones(par.nu,1);
x0(2*par.nu+1:3*par.nu)=ksipar.INIT_DET_T*tref*ones(par.nu,1);

%--data error stdev
std_dobs=ksipar.STD_STRAIN*ones(1,par.nj);
std_dobsms=ksipar.STD_SLIP;
std_dobs=std_dobs./sqrt(ksipar.GAUGE_CONF_LEV);

%--fault nodes and jauge coordinates
gcfault0=g_coord_fault;
xf0=gcfault0(:,1);
yf0=gcfault0(:,2);

gcjauge=[x_jauge y_jauge z_jauge];
gcjauge(:,1)=gcjauge(:,1)/cosd(90-ksipar.TH);
gcjauge(:,3)=0;

%--initialize figures
figure(1);
tiledlayout(ksipar.NXFIG_UMAP,ksipar.NYFIG_UMAP,'TileSpacing','tight');
figure(2);
tiledlayout(ksipar.NXFIG,ksipar.NYFIG,'TileSpacing','tight');

%--Plot slip snapshots
kfig=0;
for i=1:par.nt

    %--reconstruct u
    psi=slip_basis_fn_sc_smooth(tobs(i),x1(par.nu+1:2*par.nu),x1(2*par.nu+1:3*par.nu));
    u=psi.*x1(1:par.nu);

    psi0=slip_basis_fn_sc_smooth(tobs(i),x0(par.nu+1:2*par.nu),x0(2*par.nu+1:3*par.nu));
    u0=psi0.*x0(1:par.nu);

    %--compute spatial average of u
    mean_slip(i)=MV*u;
    mean_slip0(i)=MV*u0;

    %--compute strain
    strain(:,i)=gf_ezz*u;
    strain0(:,i)=gf_ezz*u0;


    %--plot u
    if mod(i,floor(par.nt/ksipar.NFIG_UMAP))==0
        kfig=kfig+1;

        xfig=floor((kfig-1)/ksipar.NYFIG_UMAP)+1;
        yfig=kfig-(xfig-1)*ksipar.NYFIG_UMAP;

        figure(1);
        nexttile(kfig)

        patch('Faces',g_num_fault,'Vertices',100*[xf0 yf0],'FaceVertexCData', [u]*1e3, 'FaceColor', 'interp',...
            'FaceVertexAlphaData',vec_cache_lowres,'AlphaDataMapping','scaled','FaceAlpha','interp');
        alphamap([0.7 1]);

        caxis([0 uaxmax*1e3]);

        %--plot gauge position
        hold on
        plot3(100*gcjauge(:,1),100*gcjauge(:,2),max(u)*1e3*ones(par.nj,1),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
        hold on
        if kfig==1
            tx=text(100*gcjauge(:,1)-1,100*gcjauge(:,2)-0.4,ksipar.STRG);
            set(tx,'Color',[0.7 0 0],'Fontsize',13)
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

        if kfig==ksipar.NXFIG_UMAP*ksipar.NYFIG_UMAP
            cb=colorbar;
            ylabel(cb,'slip ($\mu$m)','Interpreter','latex','Fontsize',15)
        end

        set(gca,'Fontsize',15)
        title([num2str(floor(tobs(i)*100)/100),' s'],'Interpreter','Latex')
    end



end

%--plot observations and model predictions
s0max=max([dobs';strain';strain0']);
s0min=min([dobs';strain';strain0']);

%--strain
figure(2);clf;
for j=1:par.nj


    nexttile(j)
    plot(tobs,squeeze(dobs(j,:))-std_dobs(j),'-','Color',[1 1 1],'Linewidth',2)
    hold on
    plot(tobs,squeeze(dobs(j,:))+std_dobs(j),'-','Color',[1 1 1],'Linewidth',2)
    ftest=fill([tobs' fliplr(tobs')], [squeeze(dobs(j,:))+std_dobs(j) fliplr(squeeze(dobs(j,:))-std_dobs(j))],[0.7 0.7 0.7],'FaceAlpha',0.5,'EdgeColor',[1 1 1]);
    hold on
    plot(tobs,squeeze(dobs(j,:)),'-k')
    hold on
    plot(tobs,squeeze(strain0(j,:)),'-b','Linewidth',2)
    hold on
    plot(tobs,squeeze(strain(j,:)),'-r','Linewidth',2);




    txg(j)=text(max(tobs)/100,s0min(j)+0.9*(s0max(j)-s0min(j)),['G',num2str(j)]);
    xlim([0 max(tobs)]);ylim([1.1*s0min(j) 1.1*s0max(j)]);
    if j>(ksipar.NXFIG*(ksipar.NYFIG-1))
        xlabel('time (s)','Interpreter','Latex');
    else
        set(gca,'XTickLabel',[]);
    end
    if mod(j,ksipar.NXFIG)==1
        ylabel('strain $\varepsilon_{11}$','Interpreter','Latex')
    end
    set(gca,'Fontsize',15)

end
set(txg,'Fontsize',20,'Interpreter','Latex')

%--average slip
figure(2);
nexttile(par.nj+1)
plt(1)=plot(tobs,dobsms*1e3,'-k')
hold on
plot(tobs,dobsms-std_dobsms,'-','Color',[1 1 1],'Linewidth',2)
hold on
plot(tobs,dobsms+std_dobsms,'-','Color',[1 1 1],'Linewidth',2)
ftest=fill([tobs' fliplr(tobs')], [dobsms'+std_dobsms fliplr(dobsms'-std_dobsms)],[0.7 0.7 0.7],'FaceAlpha',0.5,'EdgeColor',[1 1 1]);
hold on
plt(2)=plot(tobs,mean_slip0*1e3,'-b','Linewidth',2)
hold on
plt(3)=plot(tobs,mean_slip*1e3,'-r','Linewidth',2)

xlim([0 max(tobs)]);ylim([1.1*min(dobsms*1e3) 1.1*max(dobsms*1e3)]);
leg=legend(plt,'observed','initial model','final model','Location','NorthWest')
set(leg,'Interpreter','latex','Fontsize',15)
xlabel('Time (s)','Interpreter','Latex')
ylabel('$\delta_m$ ($\mu$m)','Interpreter','Latex')
set(gca,'Fontsize',15)

%--rescale and print figures
figure(1)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 1.5*vop(3) 1.5*vop(4)]);


figure(2)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2*vop(3) 2*vop(4)]);

print(1,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/slip_map_',ksipar.INV_FILE_ID,'_ls',num2str(ksipar.REGUL),'.png']);

print(2,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/fit_jauge_',ksipar.INV_FILE_ID,'_ls',num2str(ksipar.REGUL),'.png']);


