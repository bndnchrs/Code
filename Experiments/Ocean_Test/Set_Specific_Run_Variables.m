function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 1; % Whether or not to use the ocean model. 
MECH.DO = 1;
THERMO.DO = 1;
WAVES.DO = 0; % Wave Fracture
ADVECT.DO = 1; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 1; % Plot Diagnostics?

%% Set Thermodynamics Options
THERMO.mergefloes = 1; 

%% Set Ocean External Forcing Options

SW = 180*sin(2*pi*FSTD.time/86400);
SW(SW < 0) = 0; 

SW = 100; 
LW = 200; 

EXFORC.QSW = SW + zeros(1,OPTS.nt); % To be Q_fixed
EXFORC.QLW = LW + zeros(1,OPTS.nt); 
EXFORC.UATM = 0*EXFORC.QLW + 6; 


EXFORC.TATM = 0*EXFORC.QLW - 15; 

Hml = cos(2*pi*[0 FSTD.time] / (365*86400)); 
Hml = 0*max(Hml,0); 

EXFORC.Hml = 20 * (1 + Hml); 

% Since we aren't using liquid exchange, these can be zero
OCEAN.do_LH = 0; 
EXFORC.QATM = 0*EXFORC.QLW; 
EXFORC.PRECIP = 0*EXFORC.QLW;
EXFORC.PATM = 0*EXFORC.QLW; 


%% Initial Conditions
% Initial Distribution has all ice at one floe size. 
var = [5^2 .125^2];
% Make a Gaussian at thickness 1.5 m and size 25 m with variance var.
psi = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[100 1.5],var);
psi = psi/sum(psi(:));
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));
psi = 0*psi; 
psi(25,5) = 1; 


% Initial concentration is 50%
FSTD.psi = .9*psi/ sum(psi(:).*FSTD.dA(:)); 