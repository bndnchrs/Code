function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 1; % Whether or not to use the ocean model. 
MECH.DO = 1;
THERMO.DO = 1;
WAVES.DO = 1; % Wave Fracture
ADVECT.DO = 1; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 1; % Plot Diagnostics?

%% Get Mechanics External Forcing
OPTS.Domainwidth = 1e5; 

% Pull in the forcing fields
EXFORC = load_forcing_fields(EXFORC,OPTS,FSTD.time);


%% Set Waves Options and External Forcing
EXFORC.wavespec = zeros(OPTS.nt,length(FSTD.Rmid));

% Pick a swell wave with a peak wavelength near 100 meters
ind = find(FSTD.Rmid > 100,1); 

spec = 0*EXFORC.wavespec(1,:); 
spec = normpdf(FSTD.Rmid,100,10); 
spec = spec / sum(spec); 

EXFORC.wavespec = repmat(spec,[OPTS.nt 1]);  

%% Set Ocean External Forcing Options

OPTS.rho_water = 1000; 

EXFORC = load_seasonal_cycle(EXFORC,FSTD.time);
OCEAN.presc_evap = 0;

OCEAN.do_LH = 1; 

OCEAN.kappa_turb = 25^2 / (7*86400); 

%% Set Deep Ocean Properties and Initial Temp
OCEAN.T = -1.8;
OCEAN.S = 33; 

OCEAN.T_b = @(z) -1.8;
OCEAN.S_b = @(z) 33;

%% Set Advection Options

var = [5^2 .125^2];

psi = mvnpdf([0*FSTD.meshR(:) FSTD.meshH(:)],[0 1],var);
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));

[~,b] = find(FSTD.Rint > 5,1);

RR = FSTD.meshR; 
RR(1:b-1,:) = Inf; 

psi = psi .* (RR.^(-OPTS.plaw));
psi(isnan(psi)) = 0; 
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));
ADVECT.FSTD_in = 1*psi/ sum(psi(:).*FSTD.dA(:)); 

ADVECT.prescribe_ice_vels = 1; 

%% Initial Conditions

% Initial concentration is 75%
FSTD.psi = .75*ADVECT.FSTD_in; 