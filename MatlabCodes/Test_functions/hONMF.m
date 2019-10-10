

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% hONMF is a kernel generating both function values and gradient  
% evaluations of
% 
% h(U,V)=(beta1/2||U||_F^2+1)*(alpha2/4 ||U||_F^4+beta2/2 ||U||_F^2+1).
%
% INPUT:
%
% xk         % current point;
% opth       % structure includes required parameters;
%
% OUTPUT:
%
% hxk        % function value of f at xk
% ghki       % gradient of f at xk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [hxk, ghki] = hONMF(opth,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of ONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 3
    error('The number of input arguments is not valid');
end   

if nargout >= 3 
    error('The number of output arguments is not valid');
end

xk = varargin{1};
i  = varargin{2};
   
% ======================== Function evaluation =========================
beta1  = opth.beta1;
beta2  = opth.beta2;
alpha2 = opth.alpha2;
U      = xk{1};
V      = xk{2};

normV2 = norm(V,'fro')^2;
normU2 = norm(U,'fro')^2;
eta2   = 0.5*beta1*normU2+1;
eta1   = 0.25*alpha2*normV2^2+0.5*beta2*normV2+1;
hxk    = eta2*eta1;
    
% ======================== gradient evaluation =========================
if nargout > 1
    if i == 1
        ghki  = (beta1*eta1)*U;
    elseif i == 2
        ghki  = eta2*(alpha2*normV2+beta2)*V;
    end
end 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of ONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%