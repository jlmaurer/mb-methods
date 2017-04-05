function [emp_perc] = coverage (results, conf_level, numModels, test_type)
% this function returns the percentage of distributions that includes the
% true MDR at the specified confidence level

    [all_lbs, all_ubs] = deal(zeros(length(conf_level),numModels));
    if test_type == 1
        YY = prctile(results.Potencies, [50-conf_level/2, 50+conf_level/2]);
        if numModels>1
            all_lbs = YY(1:end/2, :); 
            all_ubs = YY(end/2+1:end, :); 
        else
            all_lbs = YY(1:end/2)'; 
            all_ubs = YY(end/2+1:end)'; 
        end
    elseif test_type == 2
        for loop = 1:numModels
            [all_lbs(:,loop), all_ubs(:,loop)] = get_CI_bounds(results.Mtest,results.fCOBE(:,loop), conf_level./100 );
        end
    else
        error('Not a valid test_type')
    end
    emp_perc = test_ci( all_ubs, all_lbs, results.truM );
end
