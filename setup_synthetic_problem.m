function [fault, data, sig, SIG] = setup_synthetic_problem(fault, data, numModels)
% This function creates the synthetic fault, slip models, data, and Green's
% functions to use to solve the 3-D synthetic tests shown in the paper. 
%
% Inputs: 
%   faults  - a struct containing the following parameters:
%       H       - transition depth (maximum depth of locking) in km
%       L       - fault length in km
%       dip     - fault dip in positive degrees
%       nve     - number of elements in the vertical direction on the fault.
%       N       - approximate number of fault patches desired
%       coup    - coupling ratio (fixed in this version)
%       max_ss  - maximum strike-slip slip rate on the fault
%       max_ds  - maximum dip-slip slip rate on the fault. 
%
%   data    - a struct containing the data-related parameters:
%       xlocmod - one of several geodetic measurement location configurations.
%                   This can also be supplied as an Mx2 matrix of positions. 
%       sd  - scalar standard deviation of Gaussian noise to be added to data 
%
%   numModels - the number of slip models to use for the coverage test
%
% Outputs: 
%   fault   - the fault struct with the following added parameters: 
%       pm      - patch model for the fault
%       maxM    - maximum possible MDR on the fault
%       truM    - actual MDR
%       nPatch  - actual number of fault patches
%       nLock   - number of patches to lock
%       slipTru - slip models used to construct synthetic data
%
%   data    - the data struct with the following added: 
%       G           - Green's function matrix
%       xysites     - geodetic measurement locations
%       data_obs    - Observed data (true + noise)
%       dd          - weighted data vector (data_obs/sd)
%       GG          - weighted Green's function matrix (G/sd)


[fault] = make_fault_geometry(fault);

% station locations
data.xysites = make_xy(data.xlocmod);  

% G matrix
data.G = make_G3D(fault.pm, data.xysites, fault.alphas); 

[data_true, fault.slipTru] = make_syn_data(numModels, fault, data);
[data.dd, data.GG, data.data_obs, sig, SIG] = data_uncertainties(data.sd, data.G, data_true);

end

% create the fault geometry
function [fault] = make_fault_geometry(fault, faults)
% This function creates a fault model and Greens' function matrix given 
% a set of parameters. 

% Fault patch parameters
fault.nhe = ceil(fault.N/fault.nve);

%The actual number of patches in the model
fault.nPatch = sum(fault.nve*fault.nhe);

% set number of locked patches
fault.nLock = floor(fault.nPatch*fault.coupling); 

%long-term slip rates in mm/yr
longterm_slip_ss = fault.max_ss*ones(fault.nPatch,1);
longterm_slip_ds = fault.max_ds*ones(fault.nPatch,1);

% rake direction for each patch
fault.alphas = atan2(longterm_slip_ds, longterm_slip_ss);

% upper bound on slip in each patch
rs = sqrt(longterm_slip_ds.^2 + longterm_slip_ss.^2);
fault.bounds = [zeros(fault.nPatch, 1), rs(:)]; 

% fault model
if nargin < 2, faults = [fault.L, fault.H, 0, -fault.dip, 0, 100,100]; end
[fault.pm, fault.patch_areas] = make_pm(fault, faults);

% max MDR and true MDR
% units here are Pa*m^3/yr = N m/yr = moment
fault.maxM = fault.mu*fault.patch_areas'*rs/1000; 
fault.truM = (fault.nLock/fault.nPatch)*fault.maxM; 
end

%% make G
function [G] = make_G3D(pm, xysites, alphas)
% make the G matrix for a specified patchmodel pm and measurement locations
% xy

%displacement-slip matrix for GPS data
dispG=make_dispG(pm,xysites);
dispG(3:3:end,:)=[]; %no vertical

% rearrange rows of G to be all eastings then all northings
newdispG=dispG;
newdispG(1:end/2,:)=dispG(1:2:end,:);
newdispG(end/2+1:end,:)=dispG(2:2:end,:);
dispG=newdispG;

% get the strike-slip and dip-slip parts and combine, weighted by the slip rake
Gss = dispG(:,1:end/2); 
Gds = dispG(:,end/2+1:end); 
Gds_hat = Gds.*repmat(sin(alphas)', size(Gds,1), 1); 
Gss_hat = Gss.*repmat(cos(alphas)', size(Gss,1), 1); 
G = Gss_hat + Gds_hat;
end

%% Make the fault model as discretized patches. 
function [pm, patch_areas] = make_pm(fault, fault_loc)
% make the patch model
%fault_loc=[length, width, *depth, dip, strike(degrees), *north offset, 
% *east offset], *depth, north and east offsets refer to location of 
% center of top edge. 
%
% NOTE: this code can only handle the simple geometries used in the tests
% in the paper; more general geometries will not necessarily be handled
% correclty here. 

fault_loc2 = fault_loc; 
fault_loc2(4) = 180+fault_loc(4);
fault_loc2(3) = fault_loc(3) - fault_loc(2)*sind(fault_loc(4)); 
fault_loc2(6) = fault_loc(6) + fault_loc(2)*cosd(fault_loc(4)); 
fault_loc = [fault_loc2, 1, 1, 0]; 

% Create slip patches
pm = patchfault(fault_loc,fault.nhe,fault.nve);
pm(:,3)=pm(:,3)+10^-4;   %sometimes a strange error if patch breaks surface

% area of each fault patch in m^2
patch_areas = pm(:,1).*pm(:,2)*1e6;     
end

%% geodetic measurement sites
% If xlocmod is a matrix of stations locations, return them. Otherwise
% xlocmod should be a scalar that specifies the desired station
% configuration. Configurations are noted in the main script that were used
% in the paper. 
function [xy] = make_xy(xlocmod)
% specify which station location configuration to use, or pass specified
% site locations. Location format is [x(:), y(:)]. 
if numel(xlocmod)>1
    xy = xlocmod; 
else
    offset = 3; % this is to move GPS stations from off directly on the fault
    switch xlocmod
        % grid of stations 
        case 1
            xs = [40:10:160]+offset;
            ys = 40:10:160;
        
        % single station, close case
        case 2     
            xs= 70; 
            ys = 100;
        
        % set of four stations
        case 3
            xs = [110, 90];  
            ys = [100,110]; 
            
        % line of stations
        case 4
            xs = [50:2:150];  
            ys = [100]; 
            
        % close case
        case 5
            xs =  [57:1.5:72]+22;
            ys =92:1:96;
        
        % middle case
        case 6
            xs = [57:1.5:72];
            ys =92:1:96;
            
        % far case
        case 7
            xs = [57:1.5:72]-40;
            ys =92:1:96;
            
        % a few more cases
        case 8
            xs = [25:10:95];
            ys =92:1:96;
        case 9
            xs = [127:3:142]-50;
            ys =92:2:96;
        case 10
            xy = [40+randi(120, [35,1]), 70 + randi(60, [35,1])]; 
    end

    if xlocmod ~= 10
        [x, y] = meshgrid(xs,ys);
        x = x(:); y=y(:);
        xy = [x y];
    end
    
end

end

%% Make synthetic slip models and data
function [synthetic_data, synthetic_slip] = make_syn_data(numModels, fault, data)

[nData,N] = size(data.G); 

%create data matrices to fill
synthetic_slip = zeros(N, numModels);
synthetic_data = zeros(nData, numModels);

% generate synthetic slip models and data for each distribution
for loop=1:numModels
    
    switch fault.slip_model
        % random locked patches
        case 1
            index = randperm(fault.nPatch,fault.nLock);      %randomly generate index of locked patches
        % fixed set of locked patches
        case 2
            end_factor = 5; 
            index = [(fault.nPatch-fault.nPatch-end_factor):fault.nPatch-end_factor]; 
        % top two rows locked
        case 3
            index = [1:fault.nve:fault.nPatch, 2:fault.nve:fault.nPatch/2];
        % all patches locked (dislocation)
        case 4
            index = 1:fault.nPatch;
    end
    
    %make slip models
    % Note that right now, I can only handle either dip-slip or strike-slip
    % consistently, not both. 
    synthetic_slip(index,loop) = fault.bounds(1,2); 

    % Run the forward model
    synthetic_data(:,loop) = data.G*synthetic_slip(:,loop);
end

end

%% Add noise to data with standard devation sd
function [dd, GG, data, sig, SIG] = data_uncertainties(sd, G, dtru)
% This function computes the weighted d and G matrices to be used in the
% inversion process, using the vector of uncertainties sig and the matrix
% SIG. 

    % for the synthetic tests, assume uncorrelated iid variables
    [M,N] = size(G); 
    sig = sd*ones(M,1); 
    SIG = repmat(sig,1, N);
    GG = G./SIG;

    % Can handle multiple data vectors at once. Columns are different
    % realizations (vectors)
    data = dtru + sd*randn(size(dtru));
    dd = data./repmat(sig, 1, size(data,2)); 
end

function dispG=make_dispG(pm,xystats)

xloc=[xystats';zeros(1,size(xystats',2))];

npatches=size(pm,1);

for k=1:npatches
   m1=[pm(k,:) 1 0 0]';
   m2=[pm(k,:) 0 1 0]';
   
   [U1,D,S]=disloc3d(m1,xloc,1,.25);
   [U2,D,S]=disloc3d(m2,xloc,1,.25);
   
   G1(:,k)=U1(:);
   G2(:,k)=U2(:);
   
 
end

dispG=[G1 G2];
end
