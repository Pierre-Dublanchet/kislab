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
% KIS_PLT_GF: plot green's function, compute and plot resolution 
%             and restitution at specific nodes                
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--make a regular grid to interpolate resolution and make a nice contour plot
dx=2*ksipar.R/ksipar.NPTS_RES_MAP;
lx=[dx:dx:2*ksipar.R];vlx=sort([-lx 0 lx]);
ly=lx/2;vly=sort([-ly 0 ly]);
[xq,yq]=meshgrid(vlx,vly);
xq=reshape(xq,length(vlx)*length(vly),1);
yq=reshape(yq,length(vlx)*length(vly),1);

%--load green's function--%
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[nj,nf]=size(gf_ezz);

%--fault node coordinates
gcfault=g_coord_fault;

%--strain gauges coordinates
gcjauge=[x_jauge y_jauge z_jauge];
gcjauge(:,1)=gcjauge(:,1)/cosd(90-ksipar.TH);

%--choose restitution nodes to be displayed (3)-------%
figure(1);
plot(gcfault(:,1),gcfault(:,2),'ok');
title('click on 8 nodes for restitution')
[xx,yy]=ginput;
for i=1:length(xx)
    rgcfn=sqrt((xx(i)-gcfault(:,1)).^2+(yy(i)-gcfault(:,2)).^2);
    indresti(i)=find(rgcfn==min(rgcfn));
end

%--load matrix to compute spatial averages--%
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
load(file_mv)

%--compute error covariance matrix (diagonal elements)-------%
Cstrain=diag(((1/ksipar.STD_STRAIN).^2)*ones(nj,1));
Cslip=(1./ksipar.STD_SLIP).^2;

%--adjust uncertainty according to gauges confidence level---%
for i=1:nj
    Cstrain(i,i)=Cstrain(i,i)*ksipar.GAUGE_CONF_LEV(i);
end


%--compute resolution matrix
R=gf_ezz'*(Cstrain*gf_ezz)+Cslip*(MV'*MV);
DR=diag(R);
DR=DR/max(DR);
mDR = griddata(gcfault(:,1),gcfault(:,2),DR,xq,yq);
mDR=reshape(mDR,length(vlx),length(vly));

ODRnodeinj=zeros(length(indresti),length(DR));
for i=1:length(indresti)
    ODRnodeinj(i,:)=R(indresti(i),:)/max(R(indresti(i),:));
end

%--load adjacency matrix
file_adja=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat'];
load(file_adja)


%--identify "low" resolution nodes
indlowres=find(DR/max(DR)<=ksipar.RES_CUTOFF);
indhighres=find(DR/max(DR)>ksipar.RES_CUTOFF);
vec_cache_lowres=ones(size(DR));
vec_cache_lowres(indlowres)=0.1;

%--plot resolution and restitution
figure(2);clf;
subplot(ksipar.NXFIG,ksipar.NYFIG,1,'align')
patch('Faces',g_num_fault,'Vertices',100*gcfault(:,1:2),'FaceColor', 'interp','FaceVertexCData',DR/max(DR))
cb=colorbar;
xlim([-2*ksipar.R 2*ksipar.R]*100);ylim([-ksipar.R ksipar.R]*100);
hold on
plot(100*gcjauge(:,1),100*gcjauge(:,2),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
txg=text(100*gcjauge(:,1)-1,100*gcjauge(:,2)-0.4,ksipar.STRG);
set(txg,'Color',[0.7 0 0],'Fontsize',13)
hold on;
contour(vlx*100,vly*100,mDR,[ksipar.RES_CUTOFF ksipar.RES_CUTOFF],'--','Linewidth',2,'Color',[.7 0 0])
hold  on
tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(a)');
set(tx,'Fontsize',15,'Interpreter','latex')
xlabel('$x_1$ (cm)','Interpreter','latex');
ylabel('$x_2$ (cm)','Interpreter','latex');
ylabel(cb,'resolution $r/max(r)$','Interpreter','latex','Fontsize',15);
set(gca,'Fontsize',15);

for i=1:nj

    subplot(ksipar.NXFIG,ksipar.NYFIG,1+i,'align');
    patch('Faces',g_num_fault,'Vertices',100*gcfault(:,1:2),'FaceColor', 'interp','FaceVertexCData',ODRnodeinj(i,:)')
    cb=colorbar;
    xlim([-2*ksipar.R 2*ksipar.R]*100);ylim([-ksipar.R ksipar.R]*100);
    hold on
    plot(100*gcjauge(:,1),100*gcjauge(:,2),'o','MarkerFaceColor',[0.7 0 0],'MarkerEdgeColor',[0.7 0 0])
    hold on
    plot(100*gcfault(indresti(i),1),100*gcfault(indresti(i),2),'ok','MarkerFaceColor',[.7 0 .7],'MarkerSize',10)
    hold on;
    contour(vlx*100,vly*100,mDR,[ksipar.RES_CUTOFF ksipar.RES_CUTOFF],'--','Linewidth',2,'Color',[.7 0 0])
    switch i
        case{1}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(b)');
        case{2}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(c)');
        case{3}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(d)');
        case{4}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(e)');
        case{5}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(f)');
        case{6}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(g)');
        case{7}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(h)');
        case{8}
            tx=text(-2*ksipar.R*100+0.1,ksipar.R*100-0.3,'(i)');
    end
    set(tx,'Fontsize',15,'Interpreter','latex')
    xlabel('$x_1$ (cm)','Interpreter','latex');
    ylabel('$x_2$ (cm)','Interpreter','latex');
    ylabel(cb,'restitution $\rho/max(\rho)$','Interpreter','latex','Fontsize',15);
    set(gca,'Fontsize',15);

end

%--plot node number map and resolution per node
figure(3);clf
subplot(2,1,1,'align')
for i=1:length(DR)
    plot(100*gcfault(i,1),100*gcfault(i,2),'ok');
    hold on
    text(100*gcfault(i,1)+0.03,100*gcfault(i,2)+0.03,num2str(i));
end
xlabel('$x_1$ (cm)','Interpreter','latex');
ylabel('$x_2$ (cm)','Interpreter','latex');
set(gca,'Fontsize',15);

subplot(2,1,2,'align')
plot(DR,'-ok')
ylabel('resolution $r/max(r)$','Interpreter','latex');
xlabel('node number','Interpreter','latex')
set(gca,'Fontsize',15);


%--print figures
figure(2)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) 1.8*vop(3) 1.5*vop(4)]);
print(2,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/resolution_r',num2str(ksipar.HMAXF),'_Pc',num2str(ksipar.PC/1e6),'MPa.png']);

figure(3)
vop=get(gcf,'outerposition');
set(gcf,'outerposition',[0.1*vop(1) 0.1*vop(2) vop(3) 1.5*vop(4)]);
print(3,'-dpng',['./FIGURES/',ksipar.REP_MANIP,'/node_number_map',num2str(ksipar.HMAXF),'_Pc',num2str(ksipar.PC/1e6),'MPa.png'])

%--save resolution data
DR=DR/max(DR);
filesave=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/Resolution_matrix_r',num2str(ksipar.HMAXF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
save(filesave,'mDR','vec_cache_lowres','DR')
