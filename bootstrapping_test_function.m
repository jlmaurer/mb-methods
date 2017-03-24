function [results] = bootstrapping_test_function (Nboot, numModels,fault,  data,solver_opts, flag)
% This function runs the bootstrap method

% Patch vector, units such that multipled by mm/yr gives N m/yr
Aeq = fault.mu*fault.patch_areas'./1000;

% setup bootstrap loop 
all_potencies = zeros(Nboot, numModels); 
[m,n] = size(data.G);
nums = linspace(1,lbd, lbd)';

% iterate over slip models 
for outer_loop = 1:numModels
    % run bootstrap 
    for inner_loop =1:Nboot
        randomnums = datasample(nums, m);
        
        % choose whether to use data or residual bootstrap - 1 = data
        if flag==1
           dboot =  data.dd(randomnums,outer_loop);
           Gboot = data.GG(randomnums,:);
        else
            dboot = data.dd(:,outer_loop) + residuals(randomnums);
            Gboot = data.GG;
        end
        
        mhat = lsqlin(Gboot, dboot, [],[], [],[],faults.bounds(:,1),...
            faults.bounds(:,2), [], solver_opts);
        all_potencies(inner_loop,outer_loop) = Aeq*mhat;
    end
end
results.Potencies = all_potencies;
end