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

%% Set Thermo Options and External Forcing

SW = 0;% * sin(FSTD.time/43200);
SW(SW < 0) = 0; 
LW = 225; 

EXFORC.QSW = SW + zeros(1,OPTS.nt); % To be Q_fixed
EXFORC.QLW = LW + zeros(1,OPTS.nt); 
EXFORC.UATM = 0*EXFORC.QLW + 6; 

EXFORC.PATM = 0*EXFORC.QLW + 100; 
EXFORC.QATM = 0*EXFORC.QLW + 1.5e-3; 
EXFORC.PRECIP = 0*EXFORC.QLW + 0;
EXFORC.TATM = 0*EXFORC.QLW - 15; 

OCEAN.do_LH = 0; 

%% Initial Conditions
% Initial Distribution has all ice at one floe size. 
var = [5^2 .125^2];
% Make a Gaussian at thickness 1.5 m and size 25 m with variance var.
psi = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[100 1.5],var);
psi = psi/sum(psi(:));
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));

% Initial concentration is 50%
FSTD.psi = 0*psi/ sum(psi(:).*FSTD.dA(:)); 