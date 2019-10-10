

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subONMF: this function solves the following problem 
%
%   min_{z}    <nabla_i f(x),z-xi>+1/gamma_iDh(x-U_i(z-xi),x)+g_i(z), 
%
% which appears in BPALM and A-BPALM.
%
% INPUT:
%
% xk                   % point xk in the above minimization problem
% block                % the block index
% gammaki              % stepsize parameter
% gfk                  % gradient of f at xk
% ghk                  % gradient of h at xk
%
% OUTPUT:
%
% x                    % the best approximation of the minimizer          
%            
% WRITTEN BY: 
%
% Masoud Ahookhosh
% Department of Electrical Engineering(ESAT-STADIUS), KU Leuven, Belgium
%
%
% LAST UPDATE: 
%
% August 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = subONMF( beta1,alpha2,beta2,xk,optsub )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Error messages for input and output %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 5
    error('The number of input arguments is more than what is needed');
elseif nargin < 5
    error('The number of input arguments is not enough');
end

if isempty(optsub)
    error('subONMF needs the structure array optsub to be defined');
elseif ~isa(optsub,'struct')
    error('optsub should be a structure array');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Main body of subONMF.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blockFlag = optsub.block;
gammaki   = optsub.gammaki;
gfki      = optsub.gfki;
ghki      = optsub.ghki;

if blockFlag == 1
    Uk      = xk{1};
    Vk      = xk{2};
    normVk2 = norm(Vk,'fro')^2;
    eta1    = alpha2/4*normVk2^2+beta2/2*normVk2+1;
    mu1     = gammaki/(beta1*eta1);
    Uk1     = max(Uk-mu1*gfki,0);
    x       = Uk1;
elseif blockFlag == 2
    Uk        = xk{1};
    normUk2   = norm(Uk,'fro')^2;
    eta2      = 0.5*beta1*normUk2+1;
    Gk        = ghki-gammaki*gfki;
    proj_Gk   = max(Gk,0);
    Nproj_Gk2 = norm(proj_Gk,'fro')^2;
    tau1      = -beta2^2/3;
    tau2      = (-2*eta2^2*beta2^3-27*alpha2*Nproj_Gk2)/(27*eta2^2);
    tau22     = tau2/2;
    sk        = sqrt((tau22)^2+(tau1/3)^3);
    tk        = beta2/3+nthroot(-tau22+sk,3)+nthroot(-tau22-sk,3);
    Vk1       = 1/(eta2*tk)*proj_Gk;
    x         = Vk1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of subONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%