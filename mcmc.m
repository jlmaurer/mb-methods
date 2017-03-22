function [xhats, all_likes, all_dpreds, accept_rat] = mcmc(Niter, stepsize, ...
    fun, x0, k, bounds,A, b, bcut, varargin)
%MCMC: This function does MCMC sampling on fun, rejecting all samples that
%fall outside bounds, given an initial guess x0. It uses a simple uniform
%distribution as the test distribution. 
%
% Inputs: 
%       Niter:      Number of simulations to perform
%       stepsize:   characteristic step distance: Can be a scalar or column
%                   vector the same size as x0
%       fun:        function to be evaluated
%       x0:         initial guess for x (x is mx1)
%       bounds:     a mx2 matrix specifying lower (col.1) and upper (col.2)
%                   bounds on x.
%       A:          matrix of conditions: enforced as Ax =<b
%       b:          vector of conditions: enforced as Ax =<b
%       varargin:   input variables for fun
%
% Outputs: 
%       results:    the results of sampling fun 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jeremy Maurer, Stanford University, 2017
% License: MIT (see LICENSE.txt in this repository)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% specify values
N= length(x0); 

%% specify default conditions
if isempty(x0)
    x0 = (bounds(:,2) - bounds(:,1))./2;
end

if isempty(A)
    A = eye(N); 
    b = inf([N, 1]); 
end

if isempty(bcut)
    bcut = 10; 
end

%% intialize loop
xprev = x0; 
[loglikeprev, dpred] = feval(fun, x0, varargin{:});
lb = bounds(:,1); 
ub = bounds(:,2);

xhats = zeros(N, Niter/k); 
all_dpreds = zeros(length(dpred), Niter/k); 
[accept_rats, all_likes] = deal(zeros(Niter/k, 1));

accept_count = 0; 

%% Run MCMC loop
for loop=1:Niter

        x = xprev + stepsize.*(rand(size(xprev))-.5);

        mins = x<lb;
        maxs = x>ub;
        bounds = A*x; 
        test = bounds> b;
        
        % test if bounds are violated
        if any(mins) || any(maxs)
            accept =0; 
        
        % test secondary criteria
        elseif any(test)
                accept = 0; 
        else
            [loglike, dpred] = feval(fun, x, varargin{:});
            lograt = exp(loglike - loglikeprev); 
            
            if lograt>1
                accept=1;
            else
                r=rand;
                if r<lograt
                    accept=1;
                else
                    accept=0;
                end
            end
        end

        if accept==1
            xprev = x;
            loglikeprev = loglike; 
            accept_count = accept_count+1;
        else
            x = xprev;
            loglike = loglikeprev; 
        end

        %save every kth sample
        if mod(loop,k)==0
            xhats(:,loop/k) = x;
            all_likes(loop/k) = exp(loglike); 
            all_dpreds(:,loop/k) = dpred; 
            accept_rats(loop/k) = accept_count/loop; 
        end   
end

%% Trim initial burn-in period
xhats(:,1:bcut) = []; 
all_likes(1:bcut) = []; 
all_dpreds(:,1:bcut) = []; 
accept_rats(1:bcut) = []; 
accept_rat = accept_rats(end); 
end
