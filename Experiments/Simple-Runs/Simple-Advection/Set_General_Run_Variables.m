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


% Set General options
OPTS.nt = 84; % Number of timesteps
OPTS.dt = 21600; % Timestep duration
OPTS.nh = 13; % No. of thickness categories 

DH = .2; % Thickness increment (m)

% Linearly spaced time vector
OPTS.time = linspace(OPTS.dt,OPTS.nt*OPTS.dt,OPTS.nt); 

% Initial discretization. Spaced at spacing to guarantee conservation of
% volume using mechanics. 

FSTD.Rint(1) = .5;

for i = 2:65
    FSTD.Rint(i) = sqrt(2*FSTD.Rint(i-1)^2 - (4/5) * FSTD.Rint(i-1)^2);
end

FSTD.Rint = linspace(FSTD.Rint(1),FSTD.Rint(end),10); 

OPTS.nr = length(FSTD.Rint); % Number of size categories
FSTD.H = .5:.5:3; % Thickness Vector

OPTS.r_p = .5; % Minimum floe size category
OPTS.h_p = .1; % Minimum thickness category

%% First Run Initialization to get all the default fields we will use 

% Initialize the FSTD Main Parts

%% Set General Advection Options

