%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ONMF_driver.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact 
clear
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Slist = {'BPALM','A-BPALM1','A-BPALM2'};
flag  = 'dense';
%flag  = 'sparse';

lambda        = 10;
Ttol          = 15;
maxit         = 10000;
Stopping_Crit = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1
    % ==================================================================
    % Dimensions of input matrix and rank
    m = 200;
    n = 2000;
    r = 10;
    
    if strcmp(flag,'dense') 
        % Generation of synthetic data set
        U = rand(m,r);
        V = zeros(r,n);
        for k = 1 : n
            V(randi(r),k) = rand(1);
        end
        for j = 1 : r
            V(j,:) = V(j,:)/norm(V(j,:));
        end
        X = U*V; % no noise
        R = rand(m,n);
        X = X + 0.05*R/norm(R,'fro')*norm(X,'fro');
   
    elseif strcmp(flag,'sparse') 
        % Generation of synthetic data set
        V = sprand(r,n,0.5);
        U = sprand(m,r,0.5);
        X = U*V; % no noise
        R = sprand(m,n,0.5);
        X = X+0.05*R; 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% Start of impelementations %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialization
    [U0,V0] = NNDSVD(X',r,0); % SVD-based initialization
    x0(1)   = {V0'};
    x0(2)   = {U0'};
    
    nu     = 2;
    beta1  = 1;
    alpha2 = 1;
    beta2  = 1;
    L0(1)  = 2/(beta1*beta2);
    gamma1 = 1/L0(1)-eps; 
    max3t  = max([lambda/alpha2; 2*lambda/(beta1*beta2);lambda/beta2]);
    L0(2)  = 6*max3t;
    gamma2 = 1/L0(2)-eps;
    
    opt.X      = X;
    opt.lambda = lambda;
    func       = @ (varargin) ONMF(opt,varargin{:});
    
    opth.beta1  = beta1;
    opth.alpha2 = alpha2;
    opth.beta2  = beta2;
    kernel      = @ (varargin) hONMF(opth,varargin{:});
    
    subprob = @ (varargin) subONMF(beta1,alpha2,beta2,varargin{:});
    
    options.MaxNumIter     = maxit;
    options.MaxNumFunEval  = 1000;
    options.MaxNumGradEval = 5000;
    options.Stopping_Crit  = Stopping_Crit;
    options.TimeLimit      = Ttol;
    options.epsilon        = 1e-5;
    options.gamma0         = [gamma1,gamma2];
    options.lambda         = lambda;
    options.L0             = L0;
    options.nu             = nu;
  
    % ======================== applying solvers ========================
    for j = 1:length(Slist) 

        switch Slist{j}
                           
            case 'BPALM'
                fprintf('Running BPALM ...\n')
                tic
                [x,f,out] = BPALM( x0,func,kernel,subprob,options );         
                toc
                F           = out.F;
                F_best      = F(1);
                fun_BPALM   = F(end);
                t_BPALM     = out.T;
                nfunc_BPALM = out.Nfunc;
                f_eval{j}   = F';
                   
            case 'A-BPALM1'
                fprintf('Running A-BPALM1 ...\n')
                Lm             = (1e-2)*L0;
                options.L0     = Lm;
                options.gamma0 = 1./Lm-eps;
                tic
                [ x,f,out ] = A_BPALM1( x0,func,kernel,subprob,options );         
                toc
                F             = out.F;
                F_best        = F(1);
                fun_abpalm1   = F(end);
                t_abpalm1     = out.T;
                nfunc_abpalm1 = out.Nfunc;
                f_eval{j}     = F';
            
            case 'A-BPALM2'
                fprintf('Running A-BPALM2 ...\n')
                Lm             = (1e-1)*L0;
                options.L0     = Lm;
                options.gamma0 = 1./Lm-eps;
                tic
                [ x,f,out ] = A_BPALM2( x0,func,kernel,subprob,options );         
                toc
                F             = out.F;
                F_best        = F(1);
                fun_abpalm2   = F(end);
                t_abpalm2     = out.T;
                nfunc_abpalm2 = out.Nfunc;
                f_eval{j}     = F';

        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
for k = 1: length(Slist)
    switch Slist{k}                
        case 'BPALM'
            xk = floor(exp((1:20)*(log(length(f_eval{k}))/20)));
            Fk = f_eval{k}(xk);
            semilogy(xk,Fk,'Color',[0.53,0.33,0.65], ...
                        'Marker','v','LineWidth',2)
        case 'A-BPALM1'
            xk = floor(exp((1:15)*(log(length(f_eval{k}))/15)));
            Fk = f_eval{k}(xk);
            semilogy(xk,Fk,'Color',[0.83,0.36,0], ...
                     'Marker','d','LineWidth',2) 
        case 'A-BPALM2'
            xk = floor(exp((1:10)*(log(length(f_eval{k}))/10)));
            Fk = f_eval{k}(xk);
            semilogy(xk,Fk,'Color',[0.26,0.63,0.79],...
                        'Marker','p','LineWidth',2) 
    end  
    hold on
    xlabel('iterations','FontSize',20);
    ylabel('function values','FontSize',20);
    legend('BPALM','A-BPALM1','A-BPALM2','FontSize',14);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of the script %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
