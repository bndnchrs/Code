function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 0; % Whether or not to use the ocean model. 
MECH.DO = 1;
THERMO.DO = 1;
WAVES.DO = 0; % Wave Fracture
ADVECT.DO = 0; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 0; % Plot Diagnostics?

%% Set Thermo Options and External Forcing
THERMO.fixQ = 1; % Fix the heat flux

Qin = -300; 

THERMO.fixed_Q = Qin + zeros(1,OPTS.nt); % To be Q_fixed

%% Set Mechanics Options and External Forcing
EXFORC.nu = -(1/(60*86400)) * ones(OPTS.nt,2); 
EXFORC.nu(:,2) = 0; 

%% Initial Conditions
% Initial Distribution has all ice at one floe size. 
var = [5^2 .125^2];
% Make a Gaussian at thickness 1.5 m and size 25 m with variance var.

psi = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[100 1.5],var);
psi = psi/sum(psi(:));
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));


% Initial concentration is 50%
FSTD.psi = .75*psi/ sum(psi(:).*FSTD.dA(:)); 