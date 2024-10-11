## KISLAB: Kinematic Inversion of Slip in the LABoratory 

### Introduction

Package to infer quasi-static slip in space and time along a laboratory rock interface, from strain gauges and shortening measurements.
Details are provided in:

***Dublanchet, P., Passelegue, F.X., Chauris, H., Gesret, A., Twardzik, C., Noël C. (2024) Kinematic inversion of aseimic slip during the nucleation of
laboratory earthquakes, Journal of Geophysical Research: Solid Earth***

### Third party source code

The source code of KISLAB includes external source codes (see licenses
notices in licenses):

<table>
<thead>
<tr class="header">
<th style="text-align: left;"><strong>Name</strong></th>
<th style="text-align: left;"><strong>License</strong></th>
<th style="text-align: left;"><strong>URL</strong></th>
<th style="text-align: left;"><strong>Copyright</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">fminlbfgs.m</td>
<td style="text-align: left;">see licenses</td>
<td style="text-align: left;"><a href="https://fr.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer"
class="uri">https://fr.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer</a></td>
<td style="text-align: left;">Copyright (c) 2009, Dirk-Jan Kroon</td>
</tr>
<tr class="even">
<td style="text-align: left;">mcmc.m</td>
<td style="text-align: left;">see licenses</td>
<td style="text-align: left;"><a
href="https://fr.mathworks.com/matlabcentral/fileexchange/47912-markov-chain-monte-carlo-sampling-of-posterior-distribution"
class="uri">https://fr.mathworks.com/matlabcentral/fileexchange/47912-markov-chain-monte-carlo-sampling-of-posterior-distribution</a></td>
<td style="text-align: left;">Copyright (c) 2015, Aslak Grinsted</td>
</tbody>
</table>


### Running instructions

Create a main directory, with the following subdirectories:

-DATA
-FIGURES
-GAUGE_LOC
-GRADIENT_MEAN_OPERATORS
-GREEN_FUNCTIONS
-PARAM_INVERSION
-SAMPLE_MESH_FILES

Each of these subdirectories should contain a directory labeled by experiment (REPMANIP in the following).

Follow the steps (1) to (14) in the main directory. Here is a definition of some variables used in this file:

REPMANIP=name of experiment directory
HMAXF=measure of mesh resolution for the fault (inversion, small value = fine mesh)
HMAXGF=measure of mesh resolution for the computation of Green's function (inversion, small value = fine mesh)
FILEID=reference for inversion result file
REGUL=value of the regularization parameter used

#### (1) Define all parameters
KIS_PARAMETRISATION.m

#### (2) Construct fault triangulation used in inversion (sparse)

KIS_FAULT_TRI.m: 

    - input:
    
        - ./KIS_PARAM.mat (generated by KIS_PARAMETRISATION.m)
        - ./SAMPLE_MESH_FILES/REPMANIP/Fault2DRS.stl  fault geometry file obtained from GMSH for instance
        
    - output:
    
        - ./SAMPLE_MESH_FILES/REPMANIP/Nodes_xy_inversion_linear_rHMAXF.mat (nodes coordinates, surface of elements)
        - ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (adjacency matrix for fault nodes)

#### (3) Compute static Green's function:

KIS_GF_SAWCUT.m
    -input:
        -./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        -./GAUGE_LOC/REPMANIP/jauges.mat (strain gauges coordinates in the sample reference frame)
        -./SAMPLE_MESH_FILES/REPMANIP/RockSampleMesh3D.stl  rock sample geometry file obtained from GMSH for instance
        -./SAMPLE_MESH_FILES/REPMANIP/Nodes_xy_inversion_linear_rHMAXF.mat (nodes coordinates, obtained from KIS_FAULT_TRI.m )
    -output:
        -./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF_PcXXMPa.mat Green function, fault nodes coordinates, gauges coordinates, elt size

#### (4) Build gradient (d/dx, d/dy) and Spatial average operators for f(x,y) where (x,y) are fault nodes coordinates (in the fault reference frame):

KIS_GRADIENT_FEM.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (KIS_FAULT_TRI.m )
        ./SAMPLE_MESH_FILES/REPMANIP/Nodes_xy_inversion_linear_rHMAXF.mat (KIS_FAULT_TRI.m )
    output:
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/GHmatrix_rHMAXF.mat

KIS_MEAN_FEM.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (KIS_FAULT_TRI.m )
        ./SAMPLE_MESH_FILES/REPMANIP/Nodes_xy_inversion_linear_rHMAXF.mat (KIS_FAULT_TRI.m)
    output:
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat

#--------------------------------------------------------------------------------------------------------------------#
#--(5) Green's function representation and resolution analysis-------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
KIS_PLT_GF.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat (from KIS_GF_SAWCUT.m)
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat (from KIS_MEAN_FEM.m)
        ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (KIS_FAULT_TRI.m )
    output:
        ./FIGURES/REPMANIP/resolution_rHMAXF.eps(png) plot of resolution and restitution for a set of nodes
        ./FIGURES/REPMANIP/node_number_map_rHMAXF.eps(png) map plotof node numbers
        ./GREEN_FUNCTIONS/REPMANIP/Resolution_matrix_rHMAXF.mat value of resolution at each node


#--------------------------------------------------------------------#
#--(6) Deterministic inversion---------------------------------------#
#--------------------------------------------------------------------#
KIS_D.m
    includes:
        ./fminlbfgs_version2c/fminlbfgs.m
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/GHmatrix_rHMAXF.mat Gradient operator
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat Mean value operator
    output:
        ./PARAM_INVERSION/REPMANIP/inv_param_strain_rHMAXF_FILEID_lsREGUL.mat'



#----------------------------------------------------------------------#
#--(7) Compute epistemic uncertainty from the output of a KIS_D.m run--#
#----------------------------------------------------------------------#
KIS_EPISTEMIC.m
    includes:
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/REPMANIP/inv_param_strain_rHMAXF_FILEID_lsREGUL.mat (output of KIS_D.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./GRADIENT_MEAN_OPERATORS/MANIP_LONDRES/MVmatrix_rHMAXF.mat Mean value operator
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
    output:
        .print value of epistemic uncertainty for strain and slip in title of figure 2

REMARK: this should be done after a first run of KIS_D.m, then the output should be used to adjust std strain and slip values (modify KIS_PARAMETRSATION), then run again
step (1) and step (6)

#---------------------------------------------------------------------------#
#--(8) Bayesian inversion (MCMC)--------------------------------------------#
#---------------------------------------------------------------------------#
KIS_MCMC.m
    includes:
        ./slip_basis_fn_sc_smooth.m
        ./mcmc.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/REPMANIP/inv_param_strain_rHMAXF_FILEID_lsREGUL.mat (output of KIS_D.m) or specify initial model by hand
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat Mean value operator
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
    output:
        ./PARAM_INVERSION/MANIP_NICE/mcmc_param_strain_rHMAXF_FILEID_nit*.mat (mcmc models chains (chain), rejection rate (alpha, should be between 0.7 and 0.8), -log likelihood function (logP))

#------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#--(9) Plot rms evolution and acceptance rate for MCMC chain------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------#
KIS_RMS_DET_MCMC.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/MANIP_NICE/mcmc_param_strain_rHMAXF_FILEID_nit*.mat (output of KIS_MCMC.m)
    output:
        ./FIGURES/REPMANIP/rmsidetmcmc_FILEID.png 

#----------------------------------------------------------------------------------------------------------------------------------------------------#
#--(10) Post-process MCMC chain 1 (reconstruct mean slip history and slip history range from the MCMC chain)--------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------#
KIS_SLIP_MCMC.m
    includes:
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/MANIP_NICE/mcmc_param_strain_rHMAXF_FILEID_nit*.mat (output of KIS_MCMC.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat Mean value operator
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
    output:
        ./PARAM_INVERSION/REPMANIP/predictions_mcmc_rHMAXF_FILEID.mat: mean, max, and min slip, strain and average slip from mcmc chain


KIS_DENS_SLIP.m
    includes:
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/MANIP_NICE/mcmc_param_strain_rHMAXF_FILEID_nit*.mat (output of KIS_MCMC.m)
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat Mean value operator
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
    output:
        ./PARAM_INVERSION/REPMANIP/slip_dens_map_FILEID.mat: pdf of slip history for each fault node

#----------------------------------------------------------------------------------------------------------------------------------------------------#
#--(11) Plot results of deterministic inversion-------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------#
KIS_PLT_DET_SLIP.m
    includes:
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/REPMANIP/inv_param_strain_rHMAXF_FILEID_lsREGUL.mat (from KIS_D.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat Mean value operator
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
        ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (KIS_FAULT_TRI.m )
        ./GREEN_FUNCTIONS/REPMANIP/Resolution_matrix_rHMAXF.mat value of resolution at each node (from KIS_PLT_GF.m)
    output:
        ./FIGURES/REPMANIP/slip_map_FILEID_lsREGUL.png(eps): plot of slip history
        ./FIGURES/REPMANIP/fit_jauge_FILEID_lsREGUL.png(eps): plot of data vs. best model predictions (strain and slip)

#----------------------------------------------------------------------------------------------------------------------------------------------------#
#--(12) Plot results of Bayesian inversion-----------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------#
KIS_PLT_MCMC_SLIP.m
    includes:
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/REPMANIP/predictions_mcmc_rHMAXF_FILEID.mat (from KIS_SLIP_MCMC.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./GRADIENT_MEAN_OPERATORS/REPMANIP/MVmatrix_rHMAXF.mat Mean value operator
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
        ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (KIS_FAULT_TRI.m )
        ./GREEN_FUNCTIONS/REPMANIP/Resolution_matrix_rHMAXF.mat value of resolution at each node (from KIS_PLT_GF.m)
    output:
        ./FIGURES/REPMANIP/slip_map_meanmcmc_FILEID.eps(png): plot of slip map from mean mcmc model 
        ./FIGURES/REPMANIP/sliprate_map_meanmcmc_FILEID.eps(png): plot of slip rate map from mean mcmc model 
        ./FIGURES/REPMANIP/sigslip_map_meanmcmc_FILEID.eps(png) : plot of slip uncertainty (1sigma) from mcmc 
        ./FIGURES/REPMANIP/fit_jauge_meanmcmc_FILEID.eps(png): plot of data vs. best model predictions (strain and slip) 
        ./FIGURES/REPMANIP/vr_mcmc_FILEID.eps(png): rupture speed plots
        ./FIGURES/REPMANIP/stdu_mcmc_FILEID.eps(png): distribution of slip uncertainty plot

#----------------------------------------------------------------------------------------------------------------------------------------------------#
#--(13) Plot joint PDF from MCMC exploration---------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------#
KIS_PLT_PDF.m
    includes:
        ./slip_basis_fn_sc_smooth.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/MANIP_NICE/mcmc_param_strain_rHMAXF_FILEID_nit*.mat (output of KIS_MCMC.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
    output:
        ./FIGURES/REPMANIP/pdf_Tdu_mcmc_FILEID.eps(png): pdf plot delta u vs T (one panel per node)
        ./FIGURES/REPMANIP/pdf_Tto_mcmc_FILEID.eps(png): pdf plot t0 vs T (one panel per node)
        ./FIGURES/REPMANIP/pdf_todu_mcmc_FILEID.eps(png): pdf plot delta u vs t0 (one panel per node)
    
#----------------------------------------------------------------------------------------------------------------------------------------------------#
#--(14) Plot slip density (pdf of slip history)------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------#
KIS_PLT_SLIP_DENS.m
    input:
        ./KIS_PARAM.mat (from KIS_PARAMETRISATION.m)
        ./PARAM_INVERSION/REPMANIP/predictions_mcmc_rHMAXF_FILEID.mat (from KIS_SLIP_MCMC.m)
        ./PARAM_INVERSION/REPMANIP/slip_dens_map_FILEID.mat: pdf of slip history for each fault node (from KIS_DENS_SLIP.m)
        ./GREEN_FUNCTIONS/REPMANIP/strain_gf_saw_cut_sample_tr_quadratic_rHMAXF_rHMAXGF.mat Green's function
        ./DATA(OR SYNTHETIC_DATA)/REPMANIP/datafile.mat strain and slip data
        ./SAMPLE_MESH_FILES/REPMANIP/g_num_fault_sawcut_rHMAXF.mat (KIS_FAULT_TRI.m )
        ./GREEN_FUNCTIONS/REPMANIP/Resolution_matrix_rHMAXF.mat value of resolution at each node (from KIS_PLT_GF.m)
    output:
        ./FIGURES/REPMANIP/slip_density_plots_FILEID.eps(png): pdf of slip history (one panel per node)



Data file (datafile.mat) should contain three variables:

tobs :time of strain and slip measurement (1xnt vector, nt being the number of observations)
dobs :strain values (ngxnt matrix, ng being the number of strain gauges)
dobsms :mean slip values (1xnt vector)  

Strain gauge coordinates file (jauges.mat) should contain a matrix jauges (ng x 3) with
column 1: x coordinates
column 2: y coordinates
column 3: z coordinates
