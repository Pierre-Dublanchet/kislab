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
% KIS_GRADIENT_FEM: compute gradients operator                   %
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load adjacency matrix for fault nodes
file_adjacency=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/g_num_fault_sawcut_r',num2str(ksipar.HMAXF),'.mat']
load(file_adjacency);

%--load fault node coordinates and element areas
file_fnodes=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/Nodes_xy_inversion_linear_r',num2str(ksipar.HMAXF),'.mat']
load(file_fnodes);

x=xfninv;y=yfninv;
g_coord_fault=[x' y'];
S=A;Se=AE;
n=length(x);
g_num_fault=g_num_fault';

%--initialize gradient matrices
G=zeros(n,n);  %--d/dx operator
H=zeros(n,n);  %--d/dy operator

%--gradient of basis functions
A=[-1 1 0;-1 0 1];

%--loop on the nodes
for i=1:n 

    %--identify triangles attached to the node
    [iel,jel]=find(g_num_fault==i);
    Stot=sum(Se(jel));
    for j=1:length(jel)  %--loop on the triangles attached to a node
        B=g_coord_fault(g_num_fault(:,jel(j)),:);
        JAC=A*B;
        GN=inv(JAC)*A;
        G(i,g_num_fault(:,jel(j)))=G(i,g_num_fault(:,jel(j)))+squeeze(GN(1,:))*Se(jel(j))/Stot;
        H(i,g_num_fault(:,jel(j)))=H(i,g_num_fault(:,jel(j)))+squeeze(GN(2,:))*Se(jel(j))/Stot;
    end

end

%--save gradient operators
filesave=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/GHmatrix_r',num2str(ksipar.HMAXF),'.mat'];
save(filesave,'G','H');