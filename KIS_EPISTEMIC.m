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
% KIS_EPISTEMIC: compute epistemic uncertainty from the result   
%               of a preliminary deterministic step
%----------------------------------------------------------------%
clear all
close all

%--load parametrisation
load KIS_PARAM.mat

%--load results from a preliminary deterministic inversion
fileinv=['./PARAM_INVERSION/',ksipar.REP_MANIP,'/inv_param_strain_r',num2str(ksipar.HMAXF),'_',ksipar.INV_FILE_ID,'_ls',num2str(ksipar.REGUL),'.mat'];
load(fileinv);

%--load Green's function
filegf=['./GREEN_FUNCTIONS/',ksipar.REP_MANIP,'/strain_gf_saw_cut_sample_tr_quadratic_r',num2str(ksipar.HMAXF),'_r',num2str(ksipar.HMAXGF),'_Pc',num2str(ksipar.PC/1e6),'MPa.mat'];
load(filegf);
[par.nj,nu]=size(gf_ezz);

%--load spatial average operator
file_mv=['./GRADIENT_MEAN_OPERATORS/',ksipar.REP_MANIP,'/MVmatrix_r',num2str(ksipar.HMAXF),'.mat']
load(file_mv);

%--load data
filedata=['./DATA/',ksipar.REP_MANIP,'/',ksipar.FILEDATA];
load(filedata);

%--initialize slip matrix
matu=zeros(par.nu,par.nt);

for i=1:par.nt
    
    %--reconstruct u
    psi=slip_basis_fn_sc_smooth(tobs(i),x1(par.nu+1:2*par.nu),x1(2*par.nu+1:3*par.nu));
    u=psi.*x1(1:par.nu);
     
    %--keep u
    matu(:,i)=u;
    
    %--compute spatial average of u
    mean_slip(i)=MV*u;

    %--compute strain
    strain(:,i)=gf_ezz*u;

end


s0max=max([dobs';strain']);
s0min=min([dobs';strain']);
numjauge=[1:1:par.nj];

%--compute epistemic uncertainty from rms of preliminary inversion
ds=[];
for i=1:par.nj
    ds=[ds (strain(i,:)-dobs(i,:))];
end
du=(mean_slip'-dobsms);
disp(['epistemic uncertainty on strain (rms) = ',num2str(std(ds))]);
disp(['epistemic uncertainty on average slip (rms) = ',num2str(std(du*1e3)),' microns']);


%--plot data and model predictions
figure(1);clf;
for j=1:par.nj

    subplot(ksipar.NYFIG,ksipar.NXFIG,j,'align')
    plot(tobs,squeeze(strain(j,:)),'-r','Linewidth',2);
    hold on
    plot(tobs,squeeze(dobs(j,:)),'-k')
     txg(j)=text(max(tobs)/100,s0min(j)+0.9*(s0max(j)-s0min(j)),['G',num2str(j)]);
     xlim([0 max(tobs)]);ylim([1.1*s0min(j) 1.1*s0max(j)]);
     xlabel('time (s)','Interpreter','Latex')
     ylabel('axial strain $\varepsilon_{11}$','Interpreter','Latex')
     set(gca,'Fontsize',15)

end
set(txg,'Fontsize',20,'Interpreter','Latex')

figure(1);
subplot(ksipar.NYFIG,ksipar.NXFIG,par.nj+1,'align')
plot(tobs,dobsms*1e3,'-k')
hold on
plot(tobs,mean_slip*1e3,'-r','Linewidth',2)
xlim([0 max(tobs)]);ylim([1.1*min(dobsms*1e3) 1.1*max(dobsms*1e3)]);
leg=legend('observed','model','Location','NorthWest')
set(leg,'Interpreter','latex','Fontsize',15)
xlabel('Time (s)','Interpreter','Latex')
ylabel('$\delta_m$ ($\mu$m)','Interpreter','Latex')
set(gca,'Fontsize',15)

%--plot error distribution
figure(2);clf;
subplot(2,1,1,'align')
histogram(ds);
xlabel('axial strain $\varepsilon_{11}-\varepsilon_{11}^{obs}$','Interpreter','latex');
ylabel('count','Interpreter','latex');
set(gca,'Fontsize',15)

subplot(2,1,2,'align')
histogram(du*1e3);
xlabel('average slip $\delta_m-\delta_m^{obs} (\mu m)$','Interpreter','latex');
ylabel('count','Interpreter','latex');
set(gca,'Fontsize',15)




