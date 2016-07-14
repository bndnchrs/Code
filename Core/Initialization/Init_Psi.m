%% Parameters of the dual distribution

% Using the current thicknesses, update gridded size/thickness categories
% This produces things like FSTD.meshH, FSTD.meshHmid, etc
update_grid; 

% Initial Maximum Thicknesses
FSTD.H_max_i = FSTD.H_max;
FSTD.H_mid_max_i = FSTD.Hmid(end); 

FSTD.Hmid_i = FSTD.Hmid; 
FSTD.dA_i = FSTD.dA; 

% If we haven't defined psi, define it as a matrix of zeros.
% This is the FSTD!
if ~isfield(FSTD,'psi')
    
    % FSTD.psi(i,j) is the total fractional area that belongs to floes
    % that have a size between R(i) and R(i+1) and H(j) and H(j+1).
    % This means that FSTD.psi has units of area/m^2. It is not a spectrum!
    % Floe Distribution
    FSTD.psi = zeros(length(FSTD.Rint),length(FSTD.H));
    
    
end

% A matrix of ones
FSTD.one = 0*FSTD.psi + 1;

% Updating variables
% The change in a single timestep
FSTD.diff = 0*FSTD.psi;
% The change in open water fraction

FSTD.opening = 0;
FSTD.eflag = 0; % Flag for when we run into error

% The Number distribution
FSTD.NumberDist = FSTD.psi./(pi*FSTD.meshRmid.^2);

% Open Water and Ice Concentration
FSTD.conc = integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0);

FSTD.phi = 1 - FSTD.conc;

FSTD.time_now = 0; % The current model time
FSTD.i = 0; % The time step index.

% The initial mean floe thickness
OPTS.H_0 = integrate_FSTD(FSTD.psi,FSTD.Hmid,FSTD.dA,1);

if OPTS.H_0 == 0
    OPTS.H_0 = OPTS.h_p;
end

FSTD.resid_adjust = 0*FSTD.psi; 