function DIAG = Drive_FSTD(OPTS)
%% Drive_FSTD
% Updated 12/8/2015 - Chris Horvat

% This routine initializes and then executes multiple runs as the main
% driver of the FD code. It returns diagnostic output requested in the
% structure DIAG. The diagnostics may be turned on and off in the file
% Diagnostics_Init.m which appears in the Packages/Diagnostics directory,
% or by setting DIAG.DO = 0 in Set_General_Run_Variables.m in the run
% directory

% This is called from a Run_Wrapper.m file, which adds it to the path and
% executes it.

% It is written in general format for future runs

% There are five major structure files which govern the development of the
% FSTD. They are

% FSTD: which contains PSI as well as other related variables
% THERMO: containing the thermodynamic options
% MECH: similar, for mechanics
% SWELL: similar, for swell fracture
% OPTS: containing global options

% There are several structures that need to be passed.
% struct FSTD % FSTD op tions
% struct THERMO % Thermodynamics options
% struct MECH % Mechanics options
% struct WAVES % Wave fracture options
% struct OPTS % General options
% struct OCEAN % hehe . Contains information about the ocean model
% struct DIAG % Contains diagnostics
% struct EXFORC % Contains External Forcing

% Add path to Initialization Files
addpath([OPTS.path_of_code 'Core/Initialization'])

%% Start The Initialization Process
[FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_Preamble(OPTS);

%% Now this actually runs the damn thing
[FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_Run(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT) ;

% This saves it. Sent as a function to be compatible with parfor
Save_Run_Output(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT);

end

function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_Preamble(OPTS)

% This sets the general variables
% Set_General_Run_Variables is defined in the run directory
[FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = Set_General_Run_Variables(OPTS);

%% This initializes the main model components
% Set_General_Run_Variables is defined in the run directory
[FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = Initialize_Model(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT); 

%% This sets the specific variables for each model component
% Set_Specific_Run_Variables is defined in the run directory
[FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = Set_Specific_Run_Variables(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT);

%% Initialize all variables for all packages
[FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = Initialize_Packages(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT);

end
