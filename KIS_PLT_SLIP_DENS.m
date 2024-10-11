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
% KIS_PLT_SLIP_DENS: plot reconstructed slip history 
%                    for each fault node (slip density)
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load density plot matrices
filedens=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/slip_dens_map_',ksipar.INV_FILE_ID,'.mat']
load(filedens)


%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);
par.nt=length(tobs);
par.s0=max(abs(dobs'));   %--stress scale (Pa)
par.u0=max(abs(dobsms));  %--slip scale (mm)

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,par.nu]=size(gf_ezz);


%--load adjacency matrix
file_adja=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat'];
load(file_adja)

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

%--load mean and min max slips
filemeanslip=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/predictions_mcmc_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'.mat'];
load(filemeanslip)


%--initialize figures
figure(1);tiledlayout(ksipar.NXFIG_NODES,ksipar.NYFIG_NODES,'TileSpacing','none')

for i=1:nnodes

    %--position of the panel in the figure from node number
    fpos(i)=i;
    xf=floor((fpos(i)-1)/ksipar.NYFIG_NODES)+1;
    yf=fpos(i)-(xf-1)*ksipar.NYFIG_NODES;

    indnz=find(isnan(matslipdens(i).mm)==0);
    norm=sum(sum(matslipdens(i).mm(indnz)))*(ttobs(2)-ttobs(1))*(vslipbin(2)-vslipbin(1));
    pdf=matslipdens(i).mm/norm;
    x=ttobs;
    y=0.5*(vslipbin(2:end)+vslipbin(1:end-1))*1e3;

    figure(1);
    ax=nexttile(fpos(i))
    truc=imagesc(x,y,log10(pdf))
    set(truc,'AlphaData',1-isnan(matslipdens(i).mm))
    
    colormap(ax,"jet")
    caxis([-4 2.5])
    hold on
    plot(tobs,matu_mean(i,:)*1e3,'-k')
    plot(tobs,matu_mm(i,:)*1e3,'--k')
    plot(tobs,matu_pp(i,:)*1e3,'--k')
    tx=text(2,8,['node ',num2str(i),'(',num2str(floor(xf0(i)*1000)/10),',',num2str(floor(yf0(i)*1000)/10),')']);
    set(tx,'Fontsize',15,'Interpreter','latex')
    ylim([0 10.5])



    if ((xf==ksipar.NXFIG_NODES))
        xlabel('time (s)','Interpreter','latex')
    else
        set(gca,'XTickLabel',[]);
    end

    if ((yf==1))
        ylabel('slip ($\mu$m)','Interpreter','latex')
        set(gca,'YTick',[0 3 6 9]);
    else
        set(gca,'YTickLabel',[],'YTick',[0 3 6 9]);
    end
    if (fpos(i)==nnodes)
        cb=colorbar;
        ylabel(cb,'log pdf','Interpreter','latex');
    end
    set(gca,'Ydir','Normal','Fontsize',15)

end

%--final slip map
u=squeeze(matu_mean(:,length(tobs)))*1e3;
uaxmax=max(max(matu_pp))*1e3;

figure(2);
patch('Faces',g_num_fault,'Vertices',100*[xf0 yf0],'FaceVertexCData', [u], 'FaceColor', 'interp',...
    'FaceVertexAlphaData',vec_cache_lowres,'AlphaDataMapping','scaled','FaceAlpha','interp');
alphamap([0.7 1]);

caxis([0 uaxmax]);

hold on
plot3(100*gcjauge(:,1),100*gcjauge(:,2),max(u)*1e3*ones(par.nj,1),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
for i=1:nnodes
    txnn(i)=text(xf0(i)*100,yf0(i)*100,num2str(i))
end
set(txnn,'Interpreter','latex','Fontsize',15)

xlim([-2*ksipar.R 2*ksipar.R]*100);ylim([-ksipar.R ksipar.R]*100);
xlabel('$x_1$ (cm)','Interpreter','Latex');
ylabel('$x_2$ (cm)','Interpreter','Latex');
cb=colorbar;
ylabel(cb,'final slip ($\mu$m)','Interpreter','Latex')
set(gca,'Fontsize',15)

%--rescale and print figures
figure(1)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 2*vop(3) 3*vop(4)]);
filefig=['./FIGURES/',ksipar.REP_MANIP,'/slip_density_plots_',ksipar.INV_FILE_ID,'.png'];
print(1,'-dpng',filefig);

figure(2)
filefig=['./FIGURES/',ksipar.REP_MANIP,'/slip_final_map_',ksipar.INV_FILE_ID,'.png'];
print(1,'-dpng',filefig);


