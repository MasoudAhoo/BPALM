

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ONMF is a function generating both function values and gradient  
% evaluations of the convex test function:
% 
%  f(U,V) = 1/2||X-UV||_F^2+1/2||I-VV?T||_F^2+delat_{U>=0}+delat_{V>=0}.
%
% INPUT:
%
% xk         % current point;
% opt        % structure includes required parameters;
%
% OUTPUT:
%
% fxk % function value of f at xk
% gfk % gradient of f at xk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fxk, gfki] = ONMF(opt,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of ONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 3
    error('The number of input arguments is not valid');
end   

if nargout >= 3 
    error('The number of output arguments is not valid');
end

xk     = varargin{1};
i      = varargin{2};
   
% ======================== Function evaluation =========================
X      = opt.X;
lambda = opt.lambda;
U     = xk{1};
V     = xk{2};

X_UV  = X-U*V;
[m,r] = size(U);
VVt   = V*V';
I_VVt = eye(r)-VVt;
fxk   = 0.5*norm(X_UV,'fro')^2+0.5*lambda*norm(I_VVt,'fro')^2;
    
% ======================== gradient evaluation =========================
if nargout > 1
    if i == 1
        gfki = U*VVt-X*V';
    else
        gfki = (U'*U)*V-U'*X+2*lambda*(VVt*V-V);
    end
end 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of ONMF.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%