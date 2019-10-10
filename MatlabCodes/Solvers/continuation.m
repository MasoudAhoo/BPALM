

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% continuation.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,output] = continuation(contOpts) 
% continuation is a continuation procedure for penalty method in
%                    min f(x) + sum_{i=1}^N g_i(x_i) 
%      where x=(x_1,...,x_N) and
%            f is relatively smooth;
%            g_i (i=1,...,N) are proper and lsc. 
%
% INPUT:
%
% contOpts             % structure including the parameteres of scheme
%
%   .func              % function handle for the objective function
%   .kernel            % kernel of Bregman distance
%   .subprob           % function handle for associated subproblems
%   .x0                % initial point
%   .options           % structure including the parameteres
%   .method            % one of BPALM, A-BPALM1, and A-BPALM2
%   .contIter          % maximum number of iteration for inner method
%   .contTimeLimit     % maximum running time for inner method
%   .c                 % continuation constant
%   .lambda            % penalty parameter
%   .X                 % data matrix X
%   .stop_crit         % stopping criterion for continuation
%
%                      % 1 : stop if MaxNumIter is reached (default)
%                      % 2 : stop if MaxNumFunEval is reached
%                      % 3 : stop if MaxNumGradEval is reached
%                      % 4 : stop if TimeLimit is reached
%
% OUTPUT:
%
% x                    % the best approximation of the optimizer
% f                    % the best approximation of the optimum
% output               % structure including more output information
%
%   .T                 % running time
%   .Niter             % total number of iterations
%   .Nfunc             % total number of function evaluations
%   .Ngrad             % total number of gradient evaluations
%   .F                 % array including all function values             
%   .Status            % reason of termination
%
% REFERENCE: 
%
% [1] M. Ahookhosh, L.T.K. Hien, N. Gillis, and P. Patrinos, 
% Multi-block Bregman proximal alternationg liearized minimization and 
% its application to orthogonal nonnegative matrix factorization, 
% Submitted,(2019)
%           
% WRITTEN BY: 
%
% Masoud Ahookhosh
% Department of Electrical Engineering(ESAT-STADIUS), KU Leuven, Belgium
%
% Le Thi Khaneh Hien
% Department of Mathematics and Operational Research, University of Mons
% Mons, Belgium
%
% LAST UPDATE: 
%
% October 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,f,output] = continuation(contOpts)

% =================== initializing the parameters ======================
% ===== user has requested viewing the default values of "options" =====
x0            = contOpts.x0;
X             = contOpts.X;
func          = contOpts.func;
kernel        = contOpts.kernel;
subprob       = contOpts.subprob;
options       = contOpts.options;
method        = contOpts.method;
contIter      = contOpts.contIter; 
contTimeLimit = contOpts.contTimeLimit;
c             = contOpts.c;
lambda        = contOpts.lambda;
stop_crit     = contOpts.stop_crit;

[MaxNumIter,MaxNumFunEval,MaxNumGradEval,TimeLimit,flag_time, ...
                               Stopping_Crit] = Initialization(options);
                           
options.MaxNumIter = contIter;
options.TimeLimit  = contTimeLimit;
Stopping_Crit      = stop_crit;
Niter              = 0;
Nfunc              = 0;
Ngrad              = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Main body of continuation.m %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
StopFlag = 0;
F        = [];
T0       = tic;

% ======================= start of the main loop =======================
while ~StopFlag

    if strcmp(method,'BPALM')  
        [x,f,out]  = BPALM( x0,func,kernel,subprob,options );
        if stop_crit==1
            Niter = Niter+out.Niter;
        elseif stop_crit==2
            Nfunc = Nfunc+out.Nfunc;
        elseif stop_crit==3
            Ngrad = Ngrad+out.Ngrad;
        end
        F_new      = out.F;
        F          = [F,F_new'];
        x0         = x;
        lambda     = c*lambda;
        opt.X      = X;
        opt.lambda = lambda;
        func       = @ (varargin) ONMF(opt,varargin{:});
    elseif strcmp(method,'A-BPALM1')
        [x,f,out] = A_BPALM1( x0,func,kernel,subprob,options ); 
        if stop_crit==1
            Niter = Niter+out.Niter;
        elseif stop_crit==2
            Nfunc = Nfunc+out.Nfunc;
        elseif stop_crit==3
            Ngrad = Ngrad+out.Ngrad;
        end
        F_new      = out.F;
        F          = [F,F_new'];
        x0         = x;
        lambda     = c*lambda;
        opt.X      = X;
        opt.lambda = lambda;
        func       = @ (varargin) ONMF(opt,varargin{:});
    elseif strcmp(method,'A-BPALM2')
        [x,f,out] = A_BPALM2( x0,func,kernel,subprob,options ); 
        if stop_crit==1
            Niter = Niter+out.Niter;
        elseif stop_crit==2
            Nfunc = Nfunc+out.Nfunc;
        elseif stop_crit==3
            Ngrad = Ngrad+out.Ngrad;
        end
        F_new      = out.F;
        F          = [F,F_new'];
        x0         = x;
        lambda     = c*lambda;
        opt.X      = X;
        opt.lambda = lambda;
        func       = @ (varargin) ONMF(opt,varargin{:});
    end
        
    % ================== checking stopping criteria ====================
    Time = toc(T0); 
    [StopFlag, Status] = StopCriterion(Niter,Nfunc,Ngrad, ...
    Time,MaxNumIter,MaxNumFunEval,MaxNumGradEval,TimeLimit, ...
    Stopping_Crit);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Status
T             = Time;
output.T      = T;
output.F      = F';
output.Niter  = Niter;
output.Nfunc  = Nfunc;
output.Ngrad  = Ngrad;
output.Status = Status;
        
if flag_time == 1
    output.Time = Time;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of continuation.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
