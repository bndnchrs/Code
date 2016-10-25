function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) 
%% Which Processes Will Be Used?
% In order to validate a process, it must be added to this call here

FSTD.DO = 1; % Do the main model stuff
OCEAN.DO = 0; % Whether or not to use the ocean model. 
MECH.DO = 0;
THERMO.DO = 0;
WAVES.DO = 1; % Wave Fracture
ADVECT.DO = 1; % Advection Package
DIAG.DO = 1; % Diagnostics Package

%% Diagnostics Options
DIAG.DOPLOT = 1; % Plot Diagnostics?


%% Set Waves Options and External Forcing
EXFORC.wavespec = zeros(OPTS.nt,length(FSTD.Rmid));

% Pick a swell wave with a peak wavelength near 100 meters
ind = find(FSTD.Rmid > 100,1); 

spec = 0*EXFORC.wavespec(1,:); 
spec = normpdf(FSTD.Rmid,100,10); 
spec = spec / sum(spec); 

EXFORC.wavespec = repmat(spec,[OPTS.nt 1]);  


%% Set Advection Options
OPTS.Domainwidth = 1e5; 

var = [5^2 .125^2];
% Make a Gaussian at thickness 1.5 m and size 25 m with variance var.

psi = mvnpdf([0*FSTD.meshR(:) FSTD.meshH(:)],[0 1.5],var);
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));

% psi = FSTD.meshR.^(-2) ./ FSTD.dA);
[~,b] = find(FSTD.Rint > 5,1);

RR = FSTD.meshR; 
RR(1:b-1,:) = Inf; 

psi = psi .* (RR.^(-2));
psi(isnan(psi)) = 0; 
psi = reshape(psi,length(FSTD.Rint),length(FSTD.H));

ADVECT.FSTD_in = psi/ sum(psi(:).*FSTD.dA(:));  

ADVECT.prescribe_ice_vels = 1; 

%% Initial Conditions

% Initial concentration is 75%
FSTD.psi = .75*ADVECT.FSTD_in; 