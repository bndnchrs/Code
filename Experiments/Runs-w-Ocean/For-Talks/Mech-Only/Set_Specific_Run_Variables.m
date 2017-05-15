function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 0; % Whether or not to use the ocean model. 
MECH.DO = 1;
THERMO.DO = 0;
WAVES.DO = 0; % Wave Fracture
ADVECT.DO = 1; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 1; % Plot Diagnostics?


%% Set Advection Options
OPTS.Domainwidth = 1e5; 

% Pull in the forcing fields
EXFORC = load_forcing_fields(EXFORC,OPTS,FSTD.time);

var = [5^2 .125^2];

psi = mvnpdf([0*FSTD.meshR(:) FSTD.meshH(:)],[0 1],var);
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));

[~,b] = find(FSTD.Rint > 5,1);

RR = FSTD.meshR; 
RR(1:b-1,:) = Inf; 

psi = psi .* (RR.^(-2));
psi(isnan(psi)) = 0; 
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));
ADVECT.FSTD_in = 1*psi/ sum(psi(:).*FSTD.dA(:)); 

ADVECT.prescribe_ice_vels = 1; 

%% Initial Conditions

% Initial concentration is 75%
FSTD.psi = .75*ADVECT.FSTD_in; 