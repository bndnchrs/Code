function [FSTD,OPTS,THERMO,MECH,SWELL,DIAG,EXFORC,OCEAN,ADVECT] = Initialize_Model(FSTD,OPTS,THERMO,MECH,SWELL,DIAG,EXFORC,OCEAN,ADVECT)
%% Initialize_Model
% This routine initializes the main floe distribution variables, including
% the floe distribution itself and many of the diagnostics
fprintf('INITIALIZING MAIN MODEL \n');

addpath([OPTS.path_of_code '/Utilities/'])
addpath([OPTS.path_of_code '/Core/Timestepping/'])

%% Timestepping
% Times - hour/day/year/etc
OPTS.hour = 3600;
OPTS.day = 24*OPTS.hour;
OPTS.month = 30*OPTS.day;
OPTS.year = 365*OPTS.day;

OPTS.totnum = 0; % Total number of timesteps that have been run
OPTS.numSC = 0; % The number of sub-cycles in a given timestep

% Debugging
if ~isfield(OPTS,'debug')
    OPTS.debug = 0; % Flag for debugging
end

% Time stepping (s)
if ~isfield(OPTS,'dt')
    OPTS.dt = (1/8) * OPTS.day; % s
end

% Number of time steps
if ~isfield(OPTS,'nt')
    OPTS.nt = round(15*OPTS.year / OPTS.dt);
end

% First timestep
if ~isfield(OPTS,'start_it')
    OPTS.start_it = 1;
end

% Last timestep
if ~isfield(OPTS,'end_it')
    OPTS.end_it = OPTS.nt;
end

% The frequency with which we save output
if ~isfield(OPTS,'save_index')
    OPTS.save_index = Inf; % Default is not to set it
end

% Width of a grid cell (m).
if ~isfield(OPTS,'domain_width')
    OPTS.domain_width = 1e1*1e3; % 1 km
end

% Time Vector
FSTD.time = OPTS.time; %s

% Time in years
FSTD.time_years = FSTD.time/OPTS.year;

% Final Time
OPTS.tend = FSTD.time(end); %s

% iteration number
FSTD.i = 1;

%% Details of the Floe Size Distribution

% Number of size classes
if ~isfield(OPTS,'nr')
    OPTS.nr = 10;
end

% Minimum floe size (m)
if ~isfield(OPTS,'r_p')
    OPTS.r_p = .5; % m
end

% Floe size vector: linearly spaced if not already specified
if ~isfield(FSTD,'Rint')
    % Vector of Sizes
    FSTD.Rint = linspace(OPTS.r_p,OPTS.r_p + (OPTS.nr-1)*OPTS.dr,OPTS.nr);
end

% Increment between floe sizes
if ~isfield(OPTS,'dr')
    % Size Increment
    OPTS.dr = 10*OPTS.r_p; % m
end

%% Parameters of the Ice Thickness Distribution

% Number of Thickness Classes
if ~isfield(OPTS,'nh')
    OPTS.nh = 10;
end

% Smallest Thicknes (m)
if ~isfield(OPTS,'h_p')
    OPTS.h_p = .1; % m
end

% Thickness Increment (m)
if ~isfield(OPTS,'dh')
    OPTS.dh = 5*OPTS.h_p; % m
end

% Vector of Thickness
if ~isfield(FSTD,'H')
    FSTD.H = linspace(OPTS.h_p,OPTS.h_p + (OPTS.nh-1)*OPTS.dh,OPTS.nh); % m
end

% Maximum Floe Size (this does not change in time)
if ~isfield(FSTD,'R_max')
    FSTD.R_max = max(FSTD.Rint) + FSTD.Rint(end) - FSTD.Rint(end-1); %m
end

% Maximum Ice Thickness (this changes in time)
if ~isfield(FSTD,'H_max')
    % Maximum Thickness
    FSTD.H_max = max(FSTD.H) + OPTS.dh; %m
    % If there are no thickness categories, set it to the minimum thickness
    if OPTS.nh == 0
        FSTD.H_max = OPTS.h_p;
    end
end

% Initialize stuff related to FSTD.psi; 
Init_Psi; 