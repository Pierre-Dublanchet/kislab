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
% KIS_GF_SAWCUT: compute Green's function relating strain at     %
%                gauges to unit slip on fault                    %
%----------------------------------------------------------------%

clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--gauge coordinates
file_gauge=['./GAUGE_LOC/',ksipar.REP_MANIP,'/',ksipar.GAUGE_FILE]
load(file_gauge);

nj=length(ksipar.GAUGES_OK);
x0_jauge=squeeze(jauges(ksipar.GAUGES_OK,1));
y0_jauge=squeeze(jauges(ksipar.GAUGES_OK,2));
z0_jauge=squeeze(jauges(ksipar.GAUGES_OK,3));

%--coordinates of the center of gauge plane (used to correct gauge position)
xcj=-ksipar.EJ*cosd(ksipar.TH);
ycj=0;
zcj=-ksipar.EJ*sind(ksipar.TH);

%--Correct gauge position to fit in the mesh
x_jauge=x0_jauge+ksipar.ALPHAJ*(xcj-x0_jauge);
y_jauge=y0_jauge+ksipar.ALPHAJ*(ycj-y0_jauge);
z_jauge=z0_jauge+ksipar.ALPHAJ*(zcj-z0_jauge);


%--create model object
smodel=createpde('structural','static-solid');

%--import geometry from .stl file
meshfile=['SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/',ksipar.MESH_SAMPLE_FILE]
g=importGeometry(smodel,meshfile);

%--plot model geometry and strain gauges positions
figure(1);
subplot(1,2,1,'align')
pdegplot(smodel,'FaceLabels','on','FaceAlpha',0.5);
hold on
plot3(x_jauge,y_jauge,z_jauge,'or')
for i=1:nj
    hold on
    text(x_jauge(i),y_jauge(i),z_jauge(i),num2str(i))
end
xlabel('x (m)','Interpreter','latex'); 
ylabel('y (m)','Interpreter','latex'); 
zlabel('z (m)','Interpreter','latex');
title('Sample geometry','Interpreter','latex')

%--generate mesh
msh=generateMesh(smodel,'Hmax',ksipar.HMAXGF);     %--default quadratic mesh
eltsize=[msh.MinElementSize msh.MaxElementSize];

%-plot mesh
figure(1);
subplot(1,2,2,'align')
pdeplot3D(smodel)
xlabel('x (m)','Interpreter','latex'); 
ylabel('y (m)','Interpreter','latex'); 
zlabel('z (m)','Interpreter','latex');
title('Mesh for Greens function computation','Interpreter','latex')


%--find fault nodes
nodes_fault = findNodes(msh,"region","Face",3);

xn=msh.Nodes(1,:);
yn=msh.Nodes(2,:);
zn=msh.Nodes(3,:);
xfn=msh.Nodes(1,nodes_fault);
yfn=msh.Nodes(2,nodes_fault);
zfn=msh.Nodes(3,nodes_fault);

nnf=length(xfn);

%--load fault inversion nodes coordinates
file_fnodes=['./SAMPLE_MESH_FILES/',ksipar.REP_MANIP,'/Nodes_xy_inversion_linear_r',num2str(ksipar.HMAXF),'.mat']
load(file_fnodes);

nnfinv=length(xfninv);

%--set structural properties
structuralProperties(smodel,'YoungsModulus',ksipar.E,'PoissonsRatio',ksipar.NU);

%--set boundary conditions
structuralBC(smodel,'Face',2,'constraint','fixed');
structuralBoundaryLoad(smodel,'Face',1,'Pressure',ksipar.PC);

%--initialize Green's function
gf_ezz=zeros(nj,nnfinv);        % axial strain (e_zz: z is principal stress direction)
gf_exx=zeros(nj,nnfinv);      % e_xx
gf_eyy=zeros(nj,nnfinv);      % e_yy
gf_exy=zeros(nj,nnfinv);      % e_xy
gf_exz=zeros(nj,nnfinv);      % e_xz
gf_eyz=zeros(nj,nnfinv);      % e_yz

%--loop on impulse location
for k=1:nnfinv

    disp(['impulse on node ',num2str(k),' over ',num2str(nnfinv)]);

    ugf=zeros(3,nnfinv);
    ugf(1,k)=-ksipar.U0*sind(ksipar.TH);
    ugf(2,k)=0;
    ugf(3,k)=ksipar.U0*cosd(ksipar.TH);

    th=ksipar.TH;

    save GF_fault_bc.mat ugf xfninv yfninv th

    structuralBC(smodel,'Face',3,'XDisplacement',@movingxPulseFcn,'YDisplacement',0,'ZDisplacement',@movingzPulseFcn);

    R=solve(smodel);

    intrpStrain = interpolateStrain(R,x_jauge,y_jauge,z_jauge);
    %--for e_zz
    if isempty(find(isnan(intrpStrain.ezz)==1))==1
        gf_ezz(:,k)=intrpStrain.ezz;
    else
        disp('NaN value of strain found')
        pause
    end

    %--for e_xx
    if isempty(find(isnan(intrpStrain.exx)==1))==1
        gf_exx(:,k)=intrpStrain.exx;
    else
        disp('NaN value of strain found')
        pause
    end

    %--for e_yy
    if isempty(find(isnan(intrpStrain.eyy)==1))==1
        gf_eyy(:,k)=intrpStrain.eyy;
    else
        disp('NaN value of strain found')
        pause
    end

    %--for e_xy
    if isempty(find(isnan(intrpStrain.exy)==1))==1
        gf_exy(:,k)=intrpStrain.exy;
    else
        disp('NaN value of strain found')
        pause
    end

    %--for e_xz
    if isempty(find(isnan(intrpStrain.exz)==1))==1
        gf_exz(:,k)=intrpStrain.exz;
    else
        disp('NaN value of strain found')
        pause
    end

    %--for e_yz
    if isempty(find(isnan(intrpStrain.eyz)==1))==1
        gf_eyz(:,k)=intrpStrain.eyz;
    else
        disp('NaN value of strain found')
        pause
    end
end

%--fault nodes coordinates
g_coord_fault=[xfninv' yfninv'];

%--save Green's functions
filesave=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
save(filesave,'x_jauge','y_jauge','z_jauge','g_coord_fault','gf_ezz','gf_exx','gf_eyy','gf_exy',...
    'gf_exz','gf_eyz','eltsize');

function ux=movingxPulseFcn(region,state)

load GF_fault_bc.mat

xx=region.x*sind(th)-region.z*cosd(th);
yy=region.y;
zz=region.x*cosd(th)+region.z*sind(th);

F=scatteredInterpolant(xfninv',yfninv',ugf(1,:)');
ux=F(xx,yy);

end

function uz=movingzPulseFcn(region,state)

load GF_fault_bc.mat

xx=region.x*sind(th)-region.z*cosd(th);
yy=region.y;
zz=region.x*cosd(th)+region.z*sind(th);

F=scatteredInterpolant(xfninv',yfninv',ugf(3,:)');
uz=F(xx,yy);

end
