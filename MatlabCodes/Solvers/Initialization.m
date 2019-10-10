

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization is a function for initializing the parameters of AGA1,
% AGA2, UGA1, UGA2, and NESUN. If some parameters specified by the user 
% Initialization uses these parameters. Otherwise, the default values 
% will be employed.
%
% INPUT: 
%
% x0                   % initial point
% options              % structure including the parameteres of schemes
%
%   .MaxNumIter        % maximum number of iterations
%   .MaxNumFunEval     % maximum number of function evaluations
%   .MaxNumGradEval    % maximum number of gradient evaluations
%   .TimeLimit         % maximum running time 
%   .Stopping_Crit     % stopping criterion
%
%                      % 1 : stop if MaxNumIter is reached (default)
%                      % 2 : stop if MaxNumFunEval is reached
%                      % 3 : stop if MaxNumSubGradEval is reached
%                      % 4 : stop if TimeLimit is reached 
%
% OUTPUT:
%
% MaxNumIter           % maximum number of iterations
% MaxNumFunEval        % maximum number of function evaluations
% MaxNumGradEval       % maximum number of gradient evaluations
% TimeLimit            % maximum running time 
% Stopping_Crit        % stopping criterion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [MaxNumIter,MaxNumFunEval,MaxNumGradEval,TimeLimit, ...
    flag_time,Stopping_Crit] = Initialization(options)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Main body of Initialization.m %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(options,'TimeLimit') 
    TimeLimit = options.TimeLimit;
else
    TimeLimit = inf;
end

if isfield(options,'MaxNumIter') 
    MaxNumIter = options.MaxNumIter;
else
    MaxNumIter = 5000;
end

if isfield(options,'MaxNumFunEval') 
    MaxNumFunEval = options.MaxNumFunEval;
else
    MaxNumFunEval = 10000;
end

if isfield(options,'MaxNumGradEval') 
    MaxNumGradEval = options.MaxNumGradEval;
else
    MaxNumGradEval = 10000;
end

if isfield(options,'Stopping_Crit') 
    Stopping_Crit = options.Stopping_Crit;
else
    Stopping_Crit = 1;
end

if isfield(options,'flag_time') 
    flag_time = options.flag_time;
else
    flag_time = 0;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% End of Initialization.m %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%