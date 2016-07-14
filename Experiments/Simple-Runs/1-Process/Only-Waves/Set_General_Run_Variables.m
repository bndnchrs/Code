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

OPTS.nt = 24*2; 
OPTS.dt = 1800; 

%% Set General Waves Options
WAVES.maxcounts = 2; 


