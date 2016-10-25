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
OPTS.dt = 86400/24;
OPTS.nr = 90; 
OPTS.nh = 13; 

OPTS.nt = 14*24;

%% Plotting Options
DIAG.plot_realtime = 0; 
DIAG.PLOT_FSTD = 1; 
DIAG.PLOT_OCEAN = 0; 
OPTS.saveplots = 1; 