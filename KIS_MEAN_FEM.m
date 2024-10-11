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
% KIS_MEAN_FEM: compute spatial average operator                 %
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load adjacency matrix for fault nodes
file_adjacency=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat']
load(file_adjacency);

%--load fault node coordinates and elements areas
file_fnodes=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/Nodes_xy_inversion_linear_r',num2str(ksipar.HMAXF),'.mat']
load(file_fnodes);

x=xfninv;y=yfninv;
g_coord_fault=[x' y'];
S=A;Se=AE;
n=length(x);
g_num_fault=g_num_fault';
[nc,ntr]=size(g_num_fault);

%--initialize operator
MV=zeros(1,n);

%--loop on elements
for i=1:ntr
    MV(g_num_fault(:,i))=MV(g_num_fault(:,i))+(Se(i)/3)*ones(1,3);
end
MV=MV/S;

%--save average operator
filesave=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat'];
save(filesave,'MV');