function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT]  = Set_General_Run_Variables(OPTS)
% This creates general run variables to use.
% Updated 12/8/2015 - Chris Horvat


% Create the structures
FSTD = struct(); 
THERMO = struct(); 
MECH = struct(); 
WAVES = struct(); 
OCEAN = struct(); 
DIAG = struct(); 
EXFORC = struct(); 
ADVECT = struct();

%% General Options
OPTS.saveplots = 0;


length = 20* ... % years
    12* ... % months
    30* ... % days
    86400; % seconds

OPTS.nt = 20 * ... % years
    12 * ... % months
    30; % nt / month
    
OPTS.dt = length / OPTS.nt;

OPTS.nr = 90; 
OPTS.nh = 12; 


%% Set General Mechanics Options
MECH.simple_oc_sr = 0; 
MECH.prescribenu = 1; 
MECH.rafting = 1; 
MECH.ridging = 1; 
MECH.try_to_load = 1; 

%%
% Variables that support the semtner thermodynamic package
THERMO.dosemtner = 1; 
THERMO.mergefloes = 0; 
OCEAN.compute_turb_deep = 1; 

%% 
WAVES.maxcounts = 1; 

%% Plotting Options
DIAG.plot_realtime = 0; 
DIAG.PLOT_FSTD = 1; 
DIAG.PLOT_OCEAN = 1; 
OPTS.saveplots = 1; 