function [emp_perc] = coverage_cobe (results, conf_level, numModels)
% this function returns the percentage of distributions that includes the
% true MDR at the specified confidence level

    [all_lbs, all_ubs] = deal(zeros(length(conf_level),numModels));
    for loop = 1:numModels
        f = results.fCOBE(:,loop); 
        [all_lbs(:,loop), all_ubs(:,loop)] = get_CI_bounds(results.Mtest, f, conf_level );
    end
    emp_perc = test_ci( all_ubs, all_lbs, results.truM );
end
