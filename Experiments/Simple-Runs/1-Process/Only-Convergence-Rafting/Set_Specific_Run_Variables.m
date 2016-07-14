function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 0; % Whether or not to use the ocean model. 
MECH.DO = 1;
THERMO.DO = 0;
WAVES.DO = 0; % Wave Fracture
ADVECT.DO = 0; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 1; % Plot Diagnostics?

%% Set Mechanics Options and External Forcing
EXFORC.nu = -(1/(30*86400)) * ones(OPTS.nt,2); 
EXFORC.nu(:,2) = 0; 

%% Initial Conditions
% Initial Distribution has all ice at one floe size. 
var = [15^2 .5^2];
% Make a Gaussian at thickness 1.5 m and size 50 m with variance var.

psi = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[50 1],var);
psi = psi/sum(psi(:));
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));
psi = 0*psi; 
psi(end-1,end-1) = 1; 
psi(end,end) = 0; 
FSTD.psi = psi/sum(psi(:).*FSTD.dA(:));