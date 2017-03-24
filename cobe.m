function [results] = cobe (data, fault, nPs, numModels, solver_opts)
% This function runs COBE for a set of data realizations

% Test values, units of moment (N m/yr)
Mtest = linspace(0, fault.maxM, nPs);

% Patch vector, units such that multipled by mm/yr gives N m/yr
Aeq = fault.mu*fault.patch_areas'./1000;

% Run ucurve coverage
[all_chi2s,out]= deal(zeros(nPs, numModels)); 
for loop = 1:numModels
    for p = 1:nPs
        [~, all_chi2s(p,loop)] = lsqlin(data.GG,data.dd(:,loop),[], [], ...
            Aeq, Mtest(p), fault.bounds(:,1), fault.bounds(:,2),...
            [], solver_opts); 
    end
end
y = all_chi2s - repmat(min(all_chi2s),nPs,1); 
f = exp(-.5*y); 

for loop = 1:numModels
    out(:,loop) = f(:,loop)./trapz(Mtest, f(:,loop));     
end
results.fCOBE = out; 
results.Mtest = Mtest; 

end
