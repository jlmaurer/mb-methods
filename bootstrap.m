function [results] = bootstrap (Nboot, numModels,fault,  data,solver_opts, flag)
% This function runs the bootstrap method

% setup bootstrap loop 
all_potencies = zeros(Nboot, numModels); 
[m,n] = size(data.G);
nums = linspace(1,m, m)';

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
            dboot = (data.data_obs(:,outer_loop) + residuals(randomnums))./data.sig;
            Gboot = data.GG;
        end
        
        mhat = lsqlin(Gboot, dboot, [],[], [],[],fault.bounds(:,1),...
            fault.bounds(:,2), [], solver_opts);
        all_potencies(inner_loop,outer_loop) = fault.avec*mhat;
    end
end
results.Potencies = all_potencies;
end