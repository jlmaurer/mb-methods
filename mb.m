%% 1D Anti-plane strike-slip fault geodetic synthetic test script for bounding MDR
% This set of functions run different synnthetic tests for estimating MDR
% on an anti-plane fault. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%   snr        -should be a struct with type and value: type can be snr or
%                 sd, value is a scalar. sd is the standard deviation of 
%                 Gaussian noise added. 
%   slipmod    -Distribution on slip to use in the synthetic trials. Can be
%                 1 - random  locked patches partially smeared out, 2 - 
%                 single dislocation (i.e. fully locked fault), 3 - a fixed
%                 smooth slip model, 5 - random locked patches, 7 - fixed 
%                 model, 9 - random Gaussian slip (use with COBLE). Some 
%                 of these values can result in slip outside the bounds,
%                 particularly if Mtest is specified, so care is needed when 
%                 using the COBE or bootstrap methods, as their solutions
%                 assume slip is bounded. 
%   N          -Number of fault patches. Should be small enough such that
%                 the closest station to the fault is at least one patch 
%                 width away. 
%
%
% Example parameters:  
% snr = struct('type', 'snr', 'val', 40); 
% N = 20; 
% slipmod = 1;
% Mtest = 0.35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following calls will produce the same type of figures as are shown in
% the paper (Bounding the moment deficit rate on crustal faults using
% geodetic data, Maurer et al. 2017). Check the function definitions for
% the parameters called. The parameters used here, along with those defined
% above, should produce close replicas (with some random variation) of the
% figures in the paper denoted below. 
%
% Fig. 2:   [out] = mb('simple_bootstrap_test');
% Fig. 4a:  [out] = mb('nloop_mcmc',0.55,snr, 0.85, slipmod);
% Fig. 6a:  [out] = mb('xloop_cobe',N, snr, 0.35, slipmod);
% Fig. 6b:  [out] = mb('nloop_cobe',0.55, snr, 0.75, slipmod);
% Fig. 7:   [out] = mb('coble_cov',0.3, N);
% Fig. 7:   [out] = mb('coble_cov',[0.1, 0.4, 0.7, 1, 1.5, 2, 2.5], N);
% Fig. 8:   [out] = mb('compare_cobe_coble');
% Fig. 11a: [out] = mb('singleStationResolution'); 
% Fig. 11b: [out] = mb('singleStationVarReduce');  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies: this function requires the mcmc.m function in this
% repository.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jeremy Maurer, Stanford University, 2017
% License: MIT (see LICENSE.txt in this repository)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = mb (varargin)
[varargout{1:nargout}]=feval(varargin{:}); 
end

%% Plots from Maurer et al., 2017

% Bootstrap test using the residual bootstrap method (Fig. 2)
function [out] = simple_bootstrap_test ()

    % parameters to test
    snrs = [2,200];
    snr.type = 'snr';
    xobs1 = [.04:.005:.1]; 
    xobs2 = xobs1+1; 

    % set other parameters
    N = 30; 
    M = .35; 

    % toggle warnings off or on
    warning('off', 'all')

    out = {};
    for k = 1:4
        if mod(k,2)==0, snr.val = snrs(2); else, snr.val = snrs(1); end
        if k<3, x = xobs2; else, x=xobs1; end
        out{k} = mb('test_mb','boot_resid', x, N, M, snr);
    end
    
    % plot results
    figure; 
    try
        subplot(2,2,1), histogram(out{1}.Potency, 'normalization','pdf'), hold on
        title('High noise, far station'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
        subplot(2,2,2), histogram(out{2}.Potency, 'normalization','pdf'), hold on
        title('Low noise, far stations'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
        subplot(2,2,3), histogram(out{3}.Potency, 'normalization','pdf'), hold on
        title('HIgh noise, close stations'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
        subplot(2,2,4), histogram(out{4}.Potency, 'normalization','pdf'), hold on
        title('Low noise, Close stations'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
    catch
        nbins = 15; 
        [h,x] = hist(out{1}.Potency,nbins);
        subplot(2,2,1), bar(x,h./trapz(x,h),'style','histc'),hold on
        title('High noise, far station'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
        [h,x] = hist(out{2}.Potency, nbins);
        subplot(2,2,2), bar(x,h./trapz(x,h),'style','histc'),hold on
        title('Low noise, far stations'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
        [h,x] = hist(out{3}.Potency, nbins);
        subplot(2,2,3),bar(x,h./trapz(x,h),'style','histc'),hold on
        title('HIgh noise, close stations'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
        
        [h,x] = hist(out{4}.Potencynbins);
        subplot(2,2,4), bar(x,h./trapz(x,h),'style','histc'),hold on
        title('Low noise, Close stations'), ax = axis; xlim([0,1])
        plot([out{1}.True_Potency, out{1}.True_Potency], [ax(3) ax(4)], 'k');
       
    end

end

% changing N, Naive MCMC (Fig. 4a)
% NOTE: this is a very simple MCMC sampler, with no adaptive step size,
% etc. Be careful when interpreting the results. If you have your own
% sampler, feel free to substitute it for mine. Just look at the inputs 
% to see how to call your function. 
function [out] = nloop_mcmc (xobs, snr, moment, slipmod)
    Ns = [2, 10, 20, 40, 80]; 
    lb = 0;
    ub = 1;
    Niter = 1e6; 
    params = struct();
    params.k = 100; 
    params.bcut = 100; 
    params.step = 0.1; 
    
    if moment > ub || moment < lb
        error('Error: Not a valid value for moment. Choose values between the bounds')
    end
    if nargin < 4, slipmod = 1; end
    if nargin < 5, loopdat = 0; end
    
    data = mb('make1Ddata', xobs, max(Ns), lb,ub, snr, moment, slipmod);

    [aps, all_likes] = deal(zeros(Niter./params.k - params.bcut,length(Ns))); 
    for loop = 1:length(Ns)
        N = Ns(loop); 
        params.init = (ub-lb)*ones(N,1)./2; 
        G = mb('make_G', N, 0, 1, xobs, data.slip); 
        [outdat] = mb('run_mcmc', data, N, lb, ub, Niter, params, G);
        aps(:,loop) = outdat.Potencies(:);
        all_likes(:,loop) = outdat.Likelihoods(:); 
    end
    
    out = struct('Samples', aps, 'Likelihoods', all_likes, 'Data', data, ...
        'N_vals', Ns);
    
    % plot results
    nbins = round(sqrt(Niter./params.k - params.bcut)); 
    figure; 
    hold on
    for loop = 1:length(Ns)
        try
            histogram(aps(:,loop), 'normalization', 'pdf')
        catch
            [h,x] = hist(aps(:,loop), nbins); 
            bar(x,h,'style', 'histc')
        end
    end
    legend(num2str(Ns(:)))
    title('Changing N, Naive MCMC')
    xlabel('MDR')
    ylabel('Probability')
    xlim([0,1])
    ax = axis; 
    plot([data.trueP, data.trueP], [ax(3) ax(4)], 'k'); 

end

% Changing N, COBE (Fig. 6b)
function [out]=nloop_cobe (xobs, snr, Mtest, slipmod)
    Ns = [4,10,20,40,80]; 
    lb = 0;
    ub = 1;

    data = mb('make1Ddata', xobs, Ns(end), lb, ub, snr, Mtest, slipmod); 
    for loop = 1:length(Ns)
        N = Ns(loop); 
        G = mb('make_G', N, 0, 1, xobs, data.slip); 
        all_likes(:,loop) = mb('run_cobe', data, N, lb, ub, G);
    end
    
    plotloopdata(data, Ns, 'N', all_likes)
    out = struct('Data', data, 'Likelihoods', all_likes, 'N_vals', Ns); 
end

% Changing x/D (Fig. 6a)
function [out] = xloop_cobe (N, snr, Mtest, slipmod)
    xs = [0.2, 0.5, 1.1, 1.8, 2.5]; 
    lb= 0; ub = 1;   
    
    % call functions
    for i = 1:length(xs)
        xobs = xs(i); 
        data = mb('make1Ddata', xobs, N, lb,ub, snr, Mtest, slipmod); 
        disp(['mean snr = ',num2str(abs(data.d./data.sd))])
        [likes_fmb, Ps] = mb('run_cobe', data, N, lb, ub); 
        all_likes(:,i) = likes_fmb; 
    end
      
    plotloopdata(data, xs, 'x_{obs}', all_likes)
    out = struct('Data', data, 'X_vals', xs, 'Likelihoods', all_likes, ...
        'Potencies', Ps); 
end

% Comparing COBE, COBLE for different constraints on MDR (Fig. 8)
function [out] = compare_cobe_coble (snr, N, Mtest)

    % Use three different measurement configurations that represent
    % well-constrained, moderately constrained, and poorly constrained MDR. 
    x{1} = [.1:.1:.6, .9, 1.4, 1.7,2,2.5,3,3.5];    % many measurements including far-field 
    x{2} = [0.1:.05:.6];                            % several measurements but not as far-field
    x{3} = 0.3;                                     % only a single relatively close station

    % parameters
    if nargin <1, snr.val = 20; snr.type = 'snr'; end
    if nargin < 2, N = 20; end
    if nargin < 3, Mtest = 0.3; end
    slipmod =1;
    Nstep = 200; 

    figure;  
    [liks_cobe, liks_coble,Ps] = deal(zeros(Nstep+1,3)); 
    [mu_coble, var_coble] = deal(zeros(3,1)); 
    for loop = 1:3
        xobs = x{loop};
        data{loop}= mb('make1Ddata', xobs, N, 0,1, snr, Mtest, slipmod, 1/Nstep);
        [liks_cobe(:,loop), Ps(:,loop)] = mb('run_cobe', data{loop}, N, 0, 1); 
        [liks_coble(:,loop),~, mu_coble(loop), var_coble(loop)] = run_coble (data{loop}, N); 

        subplot(3,1,loop)
        plot(Ps(:,loop), liks_cobe(:,loop), '-b', 'linewidth',2)
        hold on
        plot(Ps(:,loop), liks_coble(:,loop), '-r', 'linewidth',2)
        ax = axis; 
        plot([data{loop}.trueP, data{loop}.trueP], [ax(3) ax(4)], 'k');
        legend('COBE', 'COBLE')
    end
    out = struct('Data', data, 'Ps', Ps, 'COBE', liks_cobe, 'COBLE', liks_coble, ...
        'COBLE_mean', mu_coble, 'COBLE_var', var_coble); 
end

% MDR resolution with changing x/D and single station for 1D fault (Fig. 11)
function [out] = singleStationResolution ()
    xs = [.05, .1:.1:1, 1.5:.5:2.5]; 
    Mtest = .3; 
    N = 50;
    snr.type = 'snr';
    snr.val = 10;
    slipmod = 1; 

    figure;
    hold on
    ColorSet = mb('varycolor', length(xs));
    set(gca, 'ColorOrder', ColorSet);
    set(gcf, 'Colormap', ColorSet);

    for loop = 1:length(xs)
        x = xs(loop); 
        data = mb('make1Ddata', x, N, 0,1,snr, Mtest, slipmod);
        G = data.G; 
        [~,~,V] = svd(G); 
        r(:,loop) = diag(V(:,1)*V(:,1)'); 
        plot(r(:,loop), [1:N])    
    end
    % colorbar colors by distance but the scale does not reflect this
%     colorbar      
    title('Resolution for changing x')
    set(gca, 'ydir', 'reverse', 'xscale', 'log')
    ylabel('Patch # (increases with depth)')
    
    out = struct('X_vals', xs, 'Data', data, 'ResMats', r); 
end

% Variance reduction for the single station as a function of distance using
% different noise levels: Fig. 11b. We use a simple dislocation model for
% slip corresponding to Mtest = 1. 
function [out] = singleStationVarReduce ()

    N = 30; 
    if nargin < 2; Mtest = 0.4; end
    if nargin < 1; slipmod = 2; end
    
    % Prior variance on MDR assuming MDR in [0,1]
    priorVar = 1/12; 
    
    sds = [3.2e-3, 1.6e-2, 3.2e-2, 8e-2, 1.5e-1]; 
    snr.type='sd'; 
    xs = linspace(0.1, 5); 
    [snrs, varRatio,d] = deal(zeros(length(sds), length(xs))); 
    for loop = 1:length(sds)
        snr.val = sds(loop); 
        for loop2 = 1:length(xs)
            x = xs(loop2); 
            data = mb('make1Ddata', x, N, 0,1,snr, Mtest, slipmod);
            d(loop,loop2) = (data.G*data.slip);
            snrs(loop,loop2) = abs(d(loop,loop2)/data.sd); 
            [~,~, ~, varM] = mb('run_coble',data, N); 
            varRatio(loop,loop2) = varM./priorVar; 
        end
    end
    meansnr = mean(snrs,2)'; 
    
    out = struct('X_vals', xs, 'SIG_vals', sds, 'VarianceRatios', varRatio,...
        'SNR_vals', meansnr);
    
    % plot results
    figure; 
    plot(xs, varRatio)
    legend(num2str(round(meansnr(:))))
    xlabel('x/D')
    ylabel('\sigma^2_{post}/\sigma^2_{prior}')
    
end

%% Coverage plot functions

% Generates a coverage plot for COBLE/BEGS (Fig. 7)
function [out] = coble_cov (xobs, N)
    snr.val = 10; 
    snr.type = 'snr'; 
    lb = 0; 
    ub = 1;
    warning('off', 'all')
    numModels = 500; 
    ci = [0:.05:1];
    zscores = sqrt(2)*erfinv(ci)';
    hit = zeros(length(ci), 1); 
    
    % slipmod = 9 samples slip from the appropriate multi-variate Gaussian,
    % where the variance on each patch is scaled as in Eq. 9. Mtest not needed
    slipmod = 9;
    
    Mtru = zeros(numModels,1); 
    for mloop = 1:numModels
        data = make1Ddata(xobs,N, lb, ub, snr, [], slipmod);
        Mtru(mloop) = data.trueP;
        [lik,Ps, meanM, varM] = run_coble (data, N);
        steps = zscores*sqrt(varM); 
        test = Mtru(mloop) >= meanM - steps & Mtru(mloop) <= meanM + steps; 
        hit(test) = hit(test)+1; 
        
        % uncomment the following lines to draw
%         if (mod(mloop,20) == 0 || mloop == numModels)
%           clf; plot([0 1], [0 1], 'k--', ci, hit/mloop, 'b-');
%           drawnow;
%         end
    end
    emp_perc = hit/numModels*100; 
    conf_levels = ci*100; 
    out = struct('ep', emp_perc,'ci', conf_levels);
    
    % plot coverage plot
    figure; 
    plot(conf_levels, conf_levels, '--k')
    hold on
    plot(conf_levels, emp_perc, '-*b')
    legend('Ideal', 'COBLE')
    title(['COBLE Coverage plot for ',num2str(numModels), ' Distributions'])
    
end

%% Test each kind of inversion scheme individually
function [out] = test_mb (testtype, xobs, N,M, snr, slipmod)
% Inputs: 
%   testtype:   specifies which method to use, can be 'cobe', 'coble',
%                 'boot', 'boot_resid', or 'mcmc. 
%   xobs:       vector of measurement locations
%   N:          Number of fault patches
%   snr:        struct containing the type of noise (either SNR or noise
%                 level sd) and the value (scalar)
%   slipmod:    Synthetic slip mode to use. Can be 1, 3, 5, 7, 9. See 
%                 make1Ddata for more details.
%

    if nargin < 6, slipmod = 1; end
    
    % bounds on slip
    lb= 0; ub = 1;   
    
    % data
    data = make1Ddata(xobs,N, lb, ub, snr, M, slipmod); 

    switch testtype
        case 'cobe'
            [likes, Ps] = mb('run_cobe', data, N, lb, ub); 
            plot_results(data, likes)
            out = struct('Likelihood', likes, 'Potency', Ps, 'True_Potency', data.trueP); 
            
        case 'coble'
            [likes,Ps, meanM, sigM] = mb('run_coble',data, N);
            out = struct('Potency', Ps, 'True_Potency', data.trueP, ...
                'data', data, 'Likelihood', likes, 'Mhat', meanM, 'Msig', sigM); 
            
        case 'boot'
            % this is a data bootstrap, for residual bootstrap call
            % 'boot_resid'
            Nboot = 500; 
            [allPs] = mb('run_boot', 'data', data, N, lb, ub, Nboot); 
            histdat(allPs);
            out = struct('Potency', allPs, 'True_Potency', data.trueP, 'data', data); 
            
        case 'boot_resid'
            % residual bootstrap method
            Nboot = 500; 
            [allPs] = mb('run_boot', 'resid', data, N, lb, ub, Nboot); 
            histdat(allPs);
            out = struct('Potency', allPs, 'True_Potency', data.trueP, 'data', data);
            
        case 'mcmc'
            [out] = mb('run_mcmc', data, N, lb, ub); 
    end
end

%% Method functions

% COBE to estimate MDR
function [likes_COBE,Ps, chi2, slip, all_flags] = run_cobe (data, N, lb, ub, G)
    options= optimset('MaxIter', 3000); 
    r = length(data.Ps);
    Ps = data.Ps; 
    ubs = ub*ones(N,1); 
    lbs = lb*ones(N,1); 
    sd = data.sd; 
    if nargin < 5
        G = data.G; 
    end
    Gw = G./repmat(sd, 1, size(G, 2)); 
    dw = data.d./sd; 
    avec = double(ones(1,N))./N;
    [all_flags, chi2] = deal(zeros(r,1)); 
    slip = zeros(N,r); 
     for k = 1:r
        testP = Ps(k);
        Beq = double(testP);

        [slip(:,k), chi2(k), ~, all_flags(k)] = lsqlin(Gw, dw, [], [], avec, Beq, lbs, ubs, [], options); % use bounds on slip
     end
    chi2s = chi2 - min(chi2); 
    y = exp(-.5*chi2s); 
    area = trapz(Ps, y(:)); 
    likes_COBE = y./area;
end

% COBLE/BEGS to estimate MDR. Can export the moment/likelihoods or simply
% the mean and variance of the posterior. Note that the moment/likelihoods
% export does not extend beyond the bounds [0,1].
function [likes,Ps, meanM, varM] = run_coble (data, N, G)

    if nargin < 3, G = data.G; end

    r = length(data.Ps);
    Ps = data.Ps; 
    [f_coble] = deal(zeros(r, 1)); 
    d = data.d; 
    evec = ones(N,1);
    avec = (1/N)*evec; 

    % data error
    sd = data.sd;
    if numel(sd) == numel(data.d)
        Di = diag(1./sd.^2);
    elseif numel(sd)==1
        Di = eye(length(data.d))/(sd^2);
    else
        Di = inv(sd.^2); 
    end
        
    % Moment and slip prior (assumes moment is bounded between 0 and 1)
    slipmu = 0.5; 
    varMprior = 1/12; 
    sigslip = varMprior*N;           % scale slip variance for mesh independence
    Si = (1/sigslip)*eye(N); 
    
    % c, A
    c = G'*Di*d + slipmu*Si*evec; 
    A = Si + G'*Di*G; 
    
    % COBLE
    for i = 1:r
        testP = double(Ps(i));
        s_star = A\(c - avec*(avec'*(A\c)-testP)/(avec'*(A\avec))); 
        r1 = G*s_star - d; 
        r2 = s_star - slipmu*evec; 
        f_coble(i) = exp(-.5*(r1'*Di*r1 + r2'*Si*r2));
     end
    area = trapz(Ps, f_coble(:)); 
    likes = f_coble./area;
    
    % BEGS
    meanM = avec'*(A\c); 
    varM = avec'*(A\avec); 
end

% Bootstrap to estimate MDR, 'type' can be 'resid' or 'data' for residual
% or data bootstrap methods.
function [all_potencies] = run_boot (type, data, N, lb, ub, Nboot, G)
    if nargin < 7
        G = data.G; 
    end
    
   ubs = ub*ones(N,1); 
   lbs = lb*ones(N,1); 
   sd = mean(data.sd);
   Gw = G./sd
   dw = data.d./sd;
   nd = length(dw);
   nums = 1:nd;
   options = optimset('MaxIter', 2500); 
   all_slips=zeros(N,Nboot); 
   
   switch type
       case 'data'
           for i = 1:Nboot
               rnums = datasample(nums, nd); 
               dboot = dw(rnums); 
               Gboot = Gw(rnums,:); 
               sb = lsqlin(Gboot, dboot, [], [], [], [], lbs, ubs, [], options); 

               all_slips(:,i) = sb; 
           end
       case 'resid'
           sbest = lsqlin(Gw,dw,[],[],[],[],lbs,ubs,[],options);
           rs =(G*sbest - data.d); 
           for i = 1:Nboot
               dp = (data.d+ datasample(rs,nd))./data.sd; 
               sb = lsqlin(Gw, dp, [], [], [], [], lbs, ubs, [], options); 
               all_slips(:,i) = sb; 
           end
   end
   all_potencies = mean(all_slips,1);  
end

function [outdat] = run_mcmc (data, N, lb, ub, Niter, params, G)
    % get data parameters
    L = 1/N; 
    xi = [0:L:1-L]; 
    ubs = ub*ones(N,1); 
    lbs = lb*ones(N,1); 
    bounds = [lbs, ubs]; 
    sd = data.sd; 
        
    if nargin < 6
        G = data.G; 
    end
    % set sampling parameters
    if nargin < 5, Niter = 1e6; end
    k = params.k;
    bcut = params.bcut; 
    stepsize = params.step; 
    if isempty(params.init), s_init = (ub-lb)*ones(N,1)./2; else, s_init=params.init; end
    avec = (1/N)*[ones(N,1)]; 
    
    % run mcmc for slip
    [shats, slikes, dpreds, ~] = mcmc(Niter, stepsize, @compute_likelihood,...
        s_init, k, bounds, [], [], bcut, data.d, G, sd); 

    % convert slips to moment
    Ps = shats'*avec;
    
    % output results
    outdat= struct('Mean', mean(Ps), 'Standard_Error', std(Ps), ...
        'Likelihoods', slikes, 'Potencies', Ps, 'Predicted_data', dpreds); 
end

% best-fitting moment
function [M] = Msolve (data, N, lb, ub)
   ubs = ub*ones(N,1); 
   lbs = lb*ones(N,1); 
   sd = mean(data.sd);
   Gw = data.G./repmat(sd, size(data.G,1), size(data.G, 2)); 
   dw = data.d./sd;

   options = optimset('MaxIter', 2000); 
   s = lsqlin(Gw, dw, [], [], [], [], lbs, ubs, [], options); 
   M = mean(s,1);  
end

function [loglike, dpred, likelihood] = compute_likelihood (x, d, G, sd)
%COMPUTE_LIKELIHOOD: This function computes the data likelihood used in
%MCMC simulations. 
% 
% Inputs: 
%       x:          vector of x values (can also be a matrix where the
%                       columns are individual vectors)
%       bounds:     vectors of bounds on x
%       d:          vector of data values
%       G:          matrix of Green's functions
%       sd:         data standard deviation
%
% Outputs: 
%       loglike:    logarithm of the likehood
%       dpred:      predicted data
%       likelihood: actual likelihood (not log)
    
    likelihood = zeros(size(x,2),1); 
    for k = 1:size(x, 2)
        xtest = x(:,k); 
        
        % compute predicted data
        dpred = G*xtest; 
        
        % residuals
        resids = d - dpred;
        m= length(resids); 
        
        % compute log-likelihood
        if numel(sd) == 1
            tau = 1/(sd)^2; 
            loglike = (-m/2)*log(tau) + (-.5*tau).*(norm(resids)^2); 
        elseif isvector(sd)     
            loglike = (-.5)*(sum((resids.^2)./(m*sd.^2))); 
        else
            loglike = (-.5)*(resids'*inv(sd.^2)*resids); 
        end      
        
        % and likelihood
        likelihood(k) = exp(loglike);
    end
    
end

%% Data functions
function [data] = make1Ddata (xobs, N, lb,ub, snr, const, slipmod, stepsize)
    if nargin < 8, stepsize = []; end
    Ps = genP(lb, ub, const, 'uniform', stepsize);
    
    % Slip model to use. 1:a fixed smooth slip model. 5:random fully locked 
    % patches. 3: smoothed random locked patches. 7: fixed set of locked patches.
    % 9: slip is sampled from a multivariate Gaussian (for COBLE/BEGS). 
    % BE CAREFUL when changing modes, as some models can give slip outsides
    % the bounds for some values of M, especially M close to the upper
    % bound (1). 
    if slipmod ==1
        slip = zeros(N,1);
        slip(randperm(N,round(N*const))) = ub;
        slip = smooth(slip); 
        slip = [((const*N)/sum(slip))*slip]';
        slip = slip(:);     
    elseif slipmod == 2
        slip = ub*ones(N,1); 
    elseif slipmod == 3
        %NOTE: For some values of M, this mode can give slip outside the
        %bounds
        Ns = [1:N];
        s1 = 20*normpdf(Ns, 0, max(Ns)/10);
        s2 = 200*normpdf(Ns, max(Ns), max(Ns)/4);
        slip = s1+s2; 
        slip = [((const*N)/sum(slip))*slip]'; 
    elseif slipmod == 5
        slip = zeros(N,1);
        slip(randperm(N,round(N*const))) = ub;
        slip = [((const*N)/sum(slip))*slip]';
        slip = slip(:); 
    elseif slipmod ==7        
        frac = floor(N/4); 
        slip = zeros(N,1);
        slip(end/2:end/2+frac) = ub; 
        slip = slip(:);
    elseif slipmod == 9
        slip = normrnd(0.5, sqrt(N/12), [N,1]); 
        slip = slip(:); 
    else
        error('Not a valid slip model')
    end
    
    % make data
    [dw, Gw, trueP, ~, sd] = make_data(N, 0, 1, xobs, snr, slip);

    % data structure
    data = struct('d', dw, 'G', Gw, 'trueP', trueP, 'sd', sd, 'Ps', Ps, ...
        'slip', slip, 'x', xobs, 'ub', ub, 'lb', lb, 'SNR', snr); 
end

function [ data, G, true_potency, L, sd] = make_data (N, ztop, zbot, xobs, snr, slip)
%MAKE_DATA This function computes synthetic data and green's function
%matrix for a 1-D fault. 
%   Detailed explanation goes here

% create synthetic data - 1 station
if nargin <6
    slip = ones(N,1);
end

L = (zbot-ztop)/N;
ztops = ztop:L:zbot-L;

G = [];
for loop = 1:N
    [~,g1] = screw_disloc(90, 0, ztops(loop), L, xobs, slip(loop));
    G(:,loop) = g1; 
end

% create synthetic data and add noise
u3_syn = G*slip;

% can set sd equal to 1/snr or u3_syn/snr = constant sd or constant snr,
% respectively

if strcmp(snr.type,'snr')==1
    sd = abs(mean(u3_syn/snr.val));
elseif strcmp(snr.type,'sd')==1
    sd = snr.val; 
end
u3 = u3_syn + sd.*randn(size(u3_syn)); 

data = u3; 

true_potency = sum(slip)*L; 
end

function G = make_G (N, ztop, zbot, xobs, slip)

L = (zbot-ztop)/N;
ztops = ztop:L:zbot-L;

G = [];
for loop = 1:N
    [~, g1] = screw_disloc(90, 0, ztops(loop), L, xobs, slip(loop));
    G(:,loop) = g1; 
end

end

function [ u , G] = screw_disloc ( dip,xtop, ztop, L, x, slip)
%compute bottom of fault x and z values
xbot = xtop + L*cosd(dip);
zbot = ztop + L*sind(dip);
s = slip;
u = (-s/pi)*(atan2(zbot,(x-xbot)) - atan2(ztop,x));
G = ((-1/pi)*(atan2(zbot,(x-xbot)) - atan2(ztop,x)));
end

function Ps = genP (lb, ub, const, samptype, step, specPs)
    
    maxP = (ub-lb);
    if nargin < 6,specPs = []; end
    if nargin < 5, step = []; end
    if nargin < 4,samptype = 'uniform'; end
    if isempty(step), step = maxP/200; end
   
    %%%% Uniform sampling from all Potencies
    if strcmp(samptype, 'uniform') == 1 || strcmp(samptype, 'u') == 1
        Ps = 0:step:maxP;
    
    %%%% Non-uniform sampling from all Potencies
    elseif strcmp(samptype, 'nonu') == 1 || strcmp(samptype, 'nonuniform') == 1
        lp = max(const - .1*ub, lb);
        mlp = max(const - .2*ub, lb); 
        up = min(const + .1*ub, ub); 
        mup = min(const + .2*ub, ub); 
        stepfine = (up-lp)/20;
        stepmid = (up - lb)/15; 
        stepcoarse = (maxP/10); 
        
        Ps = [lb:stepcoarse:mlp, mlp:stepmid:lp, lp:stepfine:up, up:stepmid:mup,...
            mup:stepcoarse:ub]; 
        Ps = unique(Ps); 
        
        indzero = find(Ps == 0); 
        if length(indzero) > 1
            Ps(1:length(indzero)-1) = []; 
        end
        
    elseif strcmp(samptype, 'spec')== 1 || strcmp(samptype, 'specified') == 1
        Ps = specPs;
    else
        error('The type of Potency stepping specified is not valid');
    end
end

%% Plotting functions
function plot_results (data, likes)
    Ps = data.Ps; 
    trueP = data.trueP; 
    figure; hold on
    plot(Ps, likes)
    ax = axis; 
    plot([trueP, trueP], [ax(3),ax(4)], 'k')
    xlabel('Potency')
    ylabel('Likelihood')
end

function [h, x]= histdat (input_data, nbins)
    if nargin < 2
        nbins = floor(length(input_data)/10); 
    end
    [h1, x] = hist(input_data, nbins); 
    y = trapz(x, h1);
    h = h1/y; 
    figure; hold on
    bar(x, h, 'histc')
end

function plotloopdata (data, loop_param, lname, all_likes)
    Ps = data.Ps;
    C= {'-b', '-r', '--k', '-g', '-c', '-.k', ':xb'};
    figure; hold on
    for i = 1:length(loop_param)
        plot(Ps, all_likes(:,i), C{i})
    end
    trueP = data.trueP; 
    plot([trueP, trueP], [0, max(all_likes(:))], 'k')
    legend(num2str(loop_param'))
    xlabel('Potency')
    ylabel('Likelihood')
    str = ['Loop over ', lname,', snr =', num2str(data.SNR.val), ...
        ', x_{obs} = ', num2str(data.x)]; 
    title(str)
end

%% Functions for computing credible intervals

function [ lbs, ubs] = get_CI_bounds (x, d, ps )
%GET_CI_BOUNDS: This function computes the confidence interval given
%desired confidence level ps and given density d evaluated at every
%observation point x. 
%   Detailed explanation goes here
[pd_lo, pd_hi] = deal(zeros(size(ps)));
pdf_smooth = spline(x, d); 
XX = linspace(min(x), max(x), 1000); 
Vpdf = ppval(pdf_smooth, XX); 

% compute smoothed and interpolated ECDF
F = mb('calc_cdf_from_pdf', x, d); 
cdf_smooth = spline(x, F);
V = ppval(cdf_smooth, XX); 

% now loop thru confidence levels and estimate bounds
    for i = 1:numel(ps)
        delta = 0.5*(1 - ps(i));

        [~, loind] = min(abs(delta-V));
        [~, hiind] = min(abs((1-delta)-V)); 

        pd_lo(i) = XX(loind);
        pd_hi(i) = XX(hiind);

    end

% finally, assign bounds to output
lbs=pd_lo;
ubs= pd_hi;
end

function c = calc_cdf_from_pdf (x, p)
  x = x(:)';
  p = p(:)';
  c = [0, cumsum(0.5*(p(1:end-1) + p(2:end)).*diff(x))];
  c = c/c(end);
end

function [ empirical_conf_test ] = test_ci ( CI_ubs, CI_lbs, true_mean )
%TEST_CI: this function computes the number of times confidence interval
%excludes the true mean. 
%   Inputs: 
%       CI_ubs - Matrix of upper confidence bounds
%       CI_lbs - Matrix of lower confidence bounds
%       true_mean - true mean for the population
%
%   This function outputs the actual number of observations that fall
%   within the upper and lower confidence limits given by CI_ubs and
%   CI_lbs. These may be in scalar, vector, or matrix format. The algorithm
%   computes the number of times the true mean falls within the bounds
%   given and returns the answer as a percentage of the total trials. 

if size(CI_ubs, 1) < size(CI_ubs, 2) && size(CI_ubs, 1) == 1
    CI_ubs = CI_ubs'; 
    CI_lbs = CI_lbs'; 
end
test_stat = zeros(size(CI_ubs)); 

for num_ds_loop = 1:size(CI_ubs, 2)
    for test_loop = 1:size(CI_ubs,1)
        if true_mean>CI_ubs(test_loop,num_ds_loop) || true_mean<CI_lbs(test_loop,num_ds_loop)
            test_stat(test_loop, num_ds_loop) = 0; 
        else
            test_stat(test_loop,num_ds_loop) = 1; 
        end 
    end
end

empirical_conf_test = 100*sum(test_stat, 2)./size(test_stat, 2); 

end

function ColorSet=varycolor(NumberOfPlots)
% VARYCOLOR Produces colors with maximum variation on plots with multiple
% lines.
%
%     VARYCOLOR(X) returns a matrix of dimension X by 3.  The matrix may be
%     used in conjunction with the plot command option 'color' to vary the
%     color of lines.  
%
%     Yellow and White colors were not used because of their poor
%     translation to presentations.
% 
%     Example Usage:
%         NumberOfPlots=50;
%
%         ColorSet=varycolor(NumberOfPlots);
% 
%         figure
%         hold on;
% 
%         for m=1:NumberOfPlots
%             plot(ones(20,1)*m,'Color',ColorSet(m,:))
%         end

%Created by Daniel Helmick 8/12/2008

%Accessed by Jeremy Maurer, last accessed from 
% https://www.mathworks.com/matlabcentral/fileexchange/21050-varycolor 
% on March 20, 2017

error(nargchk(1,1,nargin))%correct number of input arguements??
error(nargoutchk(0, 1, nargout))%correct number of output arguements??

%Take care of the anomolies
if NumberOfPlots<1
    ColorSet=[];
elseif NumberOfPlots==1
    ColorSet=[0 1 0];
elseif NumberOfPlots==2
    ColorSet=[0 1 0; 0 1 1];
elseif NumberOfPlots==3
    ColorSet=[0 1 0; 0 1 1; 0 0 1];
elseif NumberOfPlots==4
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1];
elseif NumberOfPlots==5
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0];
elseif NumberOfPlots==6
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0; 0 0 0];

else %default and where this function has an actual advantage

    %we have 5 segments to distribute the plots
    EachSec=floor(NumberOfPlots/5); 
    
    %how many extra lines are there? 
    ExtraPlots=mod(NumberOfPlots,5); 
    
    %initialize our vector
    ColorSet=zeros(NumberOfPlots,3);
    
    %This is to deal with the extra plots that don't fit nicely into the
    %segments
    Adjust=zeros(1,5);
    for m=1:ExtraPlots
        Adjust(m)=1;
    end
    
    SecOne   =EachSec+Adjust(1);
    SecTwo   =EachSec+Adjust(2);
    SecThree =EachSec+Adjust(3);
    SecFour  =EachSec+Adjust(4);
    SecFive  =EachSec;

    for m=1:SecOne
        ColorSet(m,:)=[0 1 (m-1)/(SecOne-1)];
    end

    for m=1:SecTwo
        ColorSet(m+SecOne,:)=[0 (SecTwo-m)/(SecTwo) 1];
    end
    
    for m=1:SecThree
        ColorSet(m+SecOne+SecTwo,:)=[(m)/(SecThree) 0 1];
    end
    
    for m=1:SecFour
        ColorSet(m+SecOne+SecTwo+SecThree,:)=[1 0 (SecFour-m)/(SecFour)];
    end

    for m=1:SecFive
        ColorSet(m+SecOne+SecTwo+SecThree+SecFour,:)=[(SecFive-m)/(SecFive) 0 0];
    end
end
end
