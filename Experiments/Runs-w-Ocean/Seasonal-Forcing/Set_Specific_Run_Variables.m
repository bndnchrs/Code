function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 1; % Whether or not to use the ocean model. 
MECH.DO = 0;
THERMO.DO = 1;
WAVES.DO = 0; % Wave Fracture
ADVECT.DO = 0; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 1; % Plot Diagnostics?

%% Set Waves Options and External Forcing
EXFORC.wavespec = zeros(OPTS.nt,length(FSTD.Rmid));

% Pick a swell wave with a peak wavelength near 100 meters
ind = find(FSTD.Rmid > 100,1); 

EXFORC.wavespec(:,ind) = 1; 

%% Set Ocean External Forcing Options

OPTS.rho_water = 1000; 

EXFORC = load_seasonal_cycle(EXFORC,FSTD.time);
OCEAN.presc_evap = 1; 

OCEAN.do_LH = 1; 

OCEAN.kappa_turb = 25^2 / (30*86400); 

%% Set Deep Ocean Properties and Initial Temp
OCEAN.T = 0;
OCEAN.S = 32; 

OCEAN.T_b = @(z) -1.8;
OCEAN.S_b = @(z) 33;




%% Set Advection Options
OPTS.Domainwidth = 1e5; 
var = [5^2 .125^2];

psi = mvnpdf([0*FSTD.meshR(:) FSTD.meshH(:)],[0 1],var);
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H))./FSTD.dA;

ADVECT.FSTD_in = 1*psi/ sum(psi(:).*FSTD.dA(:)); 

ADVECT.prescribe_ice_vels = 1; 

OCEAN.UVEL = zeros(2,OPTS.nt); 

% UVEL(2) needs to be equal to UVEL(1) unless mechanics is turned on. 
OCEAN.UVEL(1,:) = 1;
OCEAN.UVEL(2,:) = 1; 

%% Set Mechanics Options and External Forcing
EXFORC.nu = .1 * OCEAN.UVEL' / OPTS.Domainwidth; 
EXFORC.nu(:,1) = 0; 

%% Initial Conditions
var = [5^2 .125^2];
% Make a Gaussian at thickness 1.5 m and size 25 m with variance var.

psi = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[750 1.5],var);
psi = psi/sum(psi(:));
psi = FSTD.meshR.^(-2) ./ FSTD.dA;
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));


% Initial concentration is 50%
FSTD.psi = .5*psi/ sum(psi(:).*FSTD.dA(:)); 
% FSTD.psi = 0*ADVECT.FSTD_in; 