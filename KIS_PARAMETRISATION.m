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
% KIS_PARAMETRISATION: define parameters
%----------------------------------------------------------------%
clear all
close all

ksipar.REP_MANIP='MANIP_TEST';                          %--Name of the experiment directory
ksipar.MESH_FAULT_FILE='Fault2DRS.stl';                 %--Fault geometry file (generated by GMSH)
ksipar.MESH_SAMPLE_FILE='RockSampleMesh3D.stl';         %--Sample geometry file (generated by GMSH)
ksipar.GAUGE_FILE='jauges.mat';                         %--Name of strain gauges coordinates file
ksipar.FILEDATA='data_ellcrack.mat';      %--Name of data file
ksipar.INV_FILE_ID='90MPaellcrack';                         %--Event ID
ksipar.HMAXF=0.015;                                     %--Fault mesh size control parameter (small value = fine mesh) for inversion (sparse mesh)
ksipar.HMAXGF=0.002;                                    %--Sample mesh size control parameter (small value = fine mesh) for Green's function computation (fine mesh)
ksipar.E=65e9;                                          %--Young's Modulus (Pa) 
ksipar.NU=0.25;                                         %--Poisson ratio 
ksipar.H=0.0856;                                        %--Sample height (m) 
ksipar.R=0.019765;                                          %--Sample radius (m) 
ksipar.TH=30;                                           %--Fault angle w.r.t. principal stress direction sigma_1
ksipar.PC=9e7;                                          %--Confining Pressure (Pa)
ksipar.EJ=0.0024;                                       %--Mean distance between fault and gauges (useful if not mentionned in gauges coordinates)
ksipar.ALPHAJ=0.07;                                     %--Correction factor for gauge position (to include gauges in the mesh)
ksipar.U0=0.001;                                        %--Imposed displacement magnitude (m) to compute GF
ksipar.GAUGES_OK=[1:8];                                 %--Gauges to consider
ksipar.STRG={'G1','G2','G3','G4','G5','G6','G7','G8'};  %--Name of gauges
ksipar.NGOK=length(ksipar.GAUGES_OK);                   %--Number of strain gauges to consider
ksipar.REGUL=0.1;                                       %--Value of the regularization parameter (Deterministic inversion)


ksipar.NITERINV_D=2000;                                 %--Max number of iterations for deterministic inversion
ksipar.TOLINV_D=1e-8;                                   %--Tolerance for deterministic inversion (LBFGS algotithm)
ksipar.STD_STRAIN=1.78e-6;                              %--Standard deviation on strain measurements (strain units) to compute error variance on strain 
ksipar.STD_SLIP=0.3e-3;                                 %--Standard deviation on slip measurements (mm) to compute error variance on mean slip
ksipar.GAUGE_CONF_LEV=[0.920737558509644;               %--Between 0 and 1 (confidence level on gauges: 0 bad, 1 perfect) just for operating gauges (should be same size as ksipar.GAUGES_OK)
         0.755012307207327;
         0.890726295731803;
         1;
         0.778857789308194;
         0.355270958916782;
         0.836644286765981;
         0.958266239801565];
ksipar.INIT_DET_DU=0.5;                                 %--Initial model of deterministic inversion: initial Delta u tested (for all nodes). Provide a fraction (>0 can be >1) of the maximum average slip recorded
ksipar.INIT_DET_T0=0.5;                                 %--Initial model of deterministic inversion: initial t0 tested (for all nodes). Provide a fraction (between 0 and 1) of the maximum time
ksipar.INIT_DET_T=0.3;                                  %--Initial model of det inversion: initial T tested (for all nodes). Provide a fraction (>0 can be >1) of the maximum time
ksipar.METH_INIT_MCMC='90MPaellcrack';                  %--Initial mcmc model defined by hand (empty) or from the result of the deterministic inversion (give the INV_FILE_ID). 
%--If ksipar.METH_INIT is empty:--%
ksipar.INIT_MCMC_DU=0.5;                                %--Initial model of mcmc inversion: initial Delta u tested (for all nodes). Provide a fraction (>0 can be >1) of the maximum average slip recorded
ksipar.INIT_MCMC_T0=0.5;                                %--Initial model of mcmc inversion: initial t0 tested (for all nodes). Provide a fraction (between 0 and 1) of the maximum time
ksipar.INIT_MCMC_T=0.3;                                 %--Initial model of mcmc inversion: initial T tested (for all nodes). Provide a fraction (>0 can be >1) of the maximum time
%-----------------------%
ksipar.SWITCH_STRAIN_DATA=ones(ksipar.NGOK,1);          %--1 if strain data considered in det inversion, 0 if not
ksipar.SWITCH_SLIP_DATA=1;                              %--1 if slip data considered in det inversion, 0 if not

ksipar.NITER_MCMC=2e5;                                  %--Number of mcmc iterations to build a result file
ksipar.NITER_MCMC_TOTAL=20e5;                           %--Total number of mcmc iterations
ksipar.MCMC_TIME_SAMP=50;                               %--Number of data samples used for mcmc exploration

ksipar.MCMC_STD_STEP_DU=1.5e-3;                         %--Standard deviation in Delta u (normalized by u0) for mcmc exploration step
ksipar.MCMC_STD_STEP_T0=1.5e-3;                         %--Standard deviation in t0 (normalized by tmax) for mcmc exploration step
ksipar.MCMC_STD_STEP_T=1.5e-3;                          %--Standard deviation in T (normalized by tmax) for mcmc exploration step

ksipar.MCMC_DU_MINMAX=[0 4];                            %--Min and max values of Delta u (normalized by u0) to be explored by mcmc
ksipar.MCMC_T0_MINMAX=[0 1];                            %--Min and max values of t0 (normalized by tmax) to be explored by mcmc
ksipar.MCMC_T_MINMAX=[0 4];                             %--Min and max values of T (normalized by tmax) to be explored by mcmc

ksipar.NITER_MCMC_BURNIN=4e5;                           %--Number of initial iterations to remove in mcmc chain (burn in phase)
ksipar.MCMC_CHAIN_SKIP=1;                               %--Keep every *models of the mcmc iteration in the chain
ksipar.MCMC_FRECU=2;                                    %--Reconstruct slip for every *models in the chain
ksipar.NSBIN=200;

ksipar.NXFIG=3;                                         %--Number of columns in strain gauge plots
ksipar.NYFIG=3;                                         %--Number of rows in strain gauge plots   (nxfig*nyfig=number of gauges+1)
ksipar.NXFIG_NODES=4;                                   %--Number of rows in pdf plots (or in figures with 1 panel/node)
ksipar.NYFIG_NODES=6;                                   %--Number of columns in pdf plots (or in figures with 1 panel/node) (nxfignodes*nyfignodes=number of fault nodes)

ksipar.NPTS_PDF=40;                                     %--Number of hor/vert points in the grid used to compute 2D joint PDF

ksipar.RES_CUTOFF=0.05;                                 %--Resolution limit used to make shaded areas in slip map plots
ksipar.NPTS_RES_MAP=20;                                 %--Number of grid points used to discretize resolution map (just used to make contour of resolution limit)


ksipar.NFIG_UMAP=12;                                    %--Number of panels in slip history figures
ksipar.NXFIG_UMAP=4;                                    %--Number of columns in slip history figures 
ksipar.NYFIG_UMAP=3;                                    %--Number of rows in slip history figures  (nxfigumap*nyfigumap=nfigumap)

ksipar.UTHRESH_RUPSPEED=2.5e-3;                         %--Slip threshold used to estimate rupture speed (in mm) 
ksipar.FIRSTACTIVATEDNODE=18;                           %--Reference node number to estimate the rupture speed 19 Nice90MPaEvt4, 22 Nice30MPaEvt3
ksipar.VR_REF=[100*100/(24*3600)...                     %--Reference rupture speeds (in cm/s)  (for plots)
    200*100/(24*3600) 500*100/(24*3600)];           


save KIS_PARAM.mat ksipar


