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
% KIS_FAULT_TRI: create fault mesh used for inversion
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--create model object
fmodel=createpde('structural','static-planestress');

%--import geometry from .stl file
meshfile=['SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/',ksipar.MESH_FAULT_FILE]
g=importGeometry(fmodel,meshfile);

%--plot model geometry
figure(1);
subplot(2,1,1,'align')
pdegplot(fmodel,'FaceLabels','on','FaceAlpha',0.5);
xlabel('$x_1$ (m)','Interpreter','latex')
ylabel('$x_2$ (m)','Interpreter','latex')
title('Fault geometry','Interpreter','latex')
set(gca,'Fontsize',10)

%--generate mesh
msh=generateMesh(fmodel,'Hmax',ksipar.HMAXF,'GeometricOrder','linear');

%--compuet elements area
[A,AE] = area(msh);

%-plot mesh
figure(1);
subplot(2,1,2,'align')
pdeplot(fmodel)
xlabel('$x_1$ (m)','Interpreter','latex')
ylabel('$x_2$ (m)','Interpreter','latex')
title('Fault mesh for inversion','Interpreter','latex')
set(gca,'Fontsize',10)

%--extract nodes from fault
xfninv=msh.Nodes(1,:);
yfninv=msh.Nodes(2,:);

%--save fault nodes coordinates for inversion
filesave=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/Nodes_xy_inversion_linear_r',num2str(ksipar.HMAXF),'.mat'];
save(filesave,'xfninv','yfninv','A','AE');

%--save adjacency matrix for the fault
[p,e,t]=meshToPet(msh);
g_num_fault=t(1:3,:)';
filesave=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat'];
save(filesave,'g_num_fault');