function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 1; % Whether or not to use the ocean model. 
MECH.DO = 0;
THERMO.DO = 0;
WAVES.DO = 1; % Wave Fracture
ADVECT.DO = 1; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 0; % Plot Diagnostics?

%% Set Thermo Options and External Forcing

%% Set Wave Fracture Options and External Forcing
WAVES.bandwidth = 25; 
WAVES.epscrit = 3e-5; 
WAVES.dobennetts = 1; 

%% Set Advection Options and External Forcing

var = [2.5^2 .125^2];

psi = FSTD.meshRmid .^(0); 
psi(FSTD.meshRmid < 400) = 0;
psi = psi / sum(psi(:)); 

ADVECT.FSTD_in = psi; 
ADVECT.prescribe_ice_vels = 1;
ADVECT.stressreducer = .1; 


OCEAN.UVEL = zeros(2,OPTS.nt); 
% Velocity at left edge of domain
OCEAN.UVEL(1,:) = .1;
% Velocity at right edge of domain
OCEAN.UVEL(2,:) = .1; 

%% Initial Conditions
% Initial Distribution has all ice at one floe size. 
var = [2.5^2 .125^2];
% Make a Gaussian at thickness 1.5 m and size 25 m with variance var.

psi = FSTD.meshRmid .^(-1); 
psi = psi / sum(psi(:)); 

% Initial concentration is 50%
FSTD.psi = psi/sum(psi(:));