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

OPTS.nt = 7*4*12 ;
OPTS.dt = 3600*6; 
OPTS.nr = 90; 

OPTS.Domainwidth = 1e5; 

OPTS.plot_inds = [1:4:28 OPTS.nt];

%% Set General Mechanics Options
MECH.simple_oc_sr = 0; 
MECH.prescribenu = 1; 
MECH.rafting = 1; 
MECH.ridging = 1; 
MECH.try_to_load = 1; 

%% 
WAVES.maxcounts = 1; 