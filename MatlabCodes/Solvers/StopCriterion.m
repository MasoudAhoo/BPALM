

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% StopCriterion.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StopCriterion is a function checking that one of the stopping criteria 
% holds to terminate AGA1, AGA2, UGA1, UGA2, and NESUN. It also perepare 
% the status determining why the algorithm is stopped.
%
% INPUT:
% 
% Niter                 % number of iterations 
% Nfunc                 % number of function evaluations
% Ngrad                 % number of subgradient evaluations
% Time                  % running time
% MaxNumIter            % maximum number of iterations
% MaxNumFunEval         % maximum number of function evaluations
% MaxNumGradEval        % maximum number of gradient evaluations
% TimeLimit             % maximum running time
% Stopping_Crit         % stopping criterion
%
%                       % 1 : stop if MaxNumIter is reached
%                       % 2 : stop if MaxNumFunEval is reached
%                       % 3 : stop if MaxNumGradEval is reached
%                       % 4 : stop if TimeLimit is reached
  
%
% OUTPUT:
%
% StopFlag              % 1: if one of the stopping criteria holds
%                       % 0: if none of the stopping criteria holds
% Status                % the reason of termination
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function [StopFlag, Status] = StopCriterion(Niter,Nfunc,Ngrad, ...
    Time,MaxNumIter,MaxNumFunEval,MaxNumGradEval,TimeLimit, ...
    Stopping_Crit)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Main body of StopCriterion.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Stopping_Crit
  case 1
    if Niter >= MaxNumIter
      StopFlag = 1;
      Status   = 'Maximum number of iterations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 2
    if Nfunc >= MaxNumFunEval
      StopFlag = 1;
      Status   = 'Maximum number of function evaluations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end 
  case 3
    if Ngrad >= MaxNumGradEval
      StopFlag = 1;
      Status   = 'Maximum number of gradient evaluations is reached';
    else
      StopFlag = 0;
      Status   = [];
    end
  case 4
    if Time >= TimeLimit
      StopFlag = 1;
      Status   = 'Time limit is reached';
    else
      StopFlag = 0;
      Status   = [];
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of StopCriterion.m %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
