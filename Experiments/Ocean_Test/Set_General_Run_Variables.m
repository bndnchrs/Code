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

% Variables that support the semtner thermodynamic package
THERMO.dosemtner = 1; 

OPTS.nt = 120*24; 
OPTS.dt = 3600;

DIAG.DO_PLOT_FSTD = 1; 
DIAG.DO_PLOT_OCEAN = 1; 

OCEAN.compute_turb_deep = 1; 