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
% slip_basis_fn_sc: function to reconstruct slip                 %
%----------------------------------------------------------------%

function psi=slip_basis_fn_sc_smooth(t,t0,T)

psi=zeros(length(t0),1);
%-------------------------%
%--Parametrisation num 1--%
%-------------------------%
% indp=find(t0<=t);
% if isempty(indp)==0
%     psi(indp)=1-exp((t0(indp)-t)./T(indp));
% end
%-------------------------%
%--Parametrisation num 2--%
%-------------------------%
% psi=0.5*(tanh((t-t0-T)./T)+1);
%-------------------------%
%--Parametrisation num 3--%
%-------------------------%
indp=find(t0<=t);
if isempty(indp)==0
    psi(indp)=0.5*(1-cos(2*pi*(t-t0(indp))./T(indp)));
end
indpp=find(t0+0.5*T<=t);
psi(indpp)=ones(length(indpp),1);