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
DIAG.DOPLOT = 0; % Plot Diagnostics?

%% Set Thermo Options and External Forcing
THERMO.fixQ = 0; % Fix the heat flux

if THERMO.fixQ
    
    THERMO.fixed_Q = -100 + zeros(1,OPTS.nt); % The fixed heat flux as a function of index
    
else
    % 200 heating from solar
    % 200 cooling from longwave.
    EXFORC.QSW = 350*cos(2*pi*OPTS.time/OPTS.year); 
    EXFORC.QSW(EXFORC.QSW < 0) = 0; 
    EXFORC.QLW = 200 + zeros(1,OPTS.nt);
    
end



%% Initial Conditions
% Initial Distribution has all ice at one floe size. 
var = [2.5^2 .125^2];

ps1 = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[15 1.5],var);
ps2 = mvnpdf([FSTD.meshR(:) FSTD.meshH(:)],[25 1.5],var);
psi = ps2/sum(ps2(:));
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));
FSTD.psi = .5*psi/sum(psi(:));

OPTS.H_0 = sum_FSTD(FSTD.psi,FSTD.Hmid,1);