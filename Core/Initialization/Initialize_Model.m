function [FSTD,OPTS,THERMO,MECH,SWELL,DIAG,EXFORC,OCEAN,ADVECT] = Initialize_Model(FSTD,OPTS,THERMO,MECH,SWELL,DIAG,EXFORC,OCEAN,ADVECT)
%% Initialize_Model
% This routine initializes the main floe distribution variables, including
% the floe distribution itself and many of the diagnostics
fprintf('INITIALIZING MAIN MODEL \n');

addpath([OPTS.path_of_code '/Utilities/'])
addpath([OPTS.path_of_code '/Core/Timestepping/'])
addpath([OPTS.path_of_code '/Core/Plotting/'])

% Add output folders

% This gets the folder in which figpath is contained.
[a,~,~] = fileparts(OPTS.figpath);

[~,~,~] = mkdir(OPTS.savepath);

[~,~,~] = mkdir(a);



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
    OPTS.dt = 3600; % s
end


% Number of time steps
if ~isfield(OPTS,'nt')
    OPTS.nt = 24*30;
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

if ~isfield(OPTS,'time')
    % Linearly spaced time vector
    OPTS.time = linspace(OPTS.dt,OPTS.nt*OPTS.dt,OPTS.nt);
end

if ~isfield(OPTS,'saveplots')
    
    OPTS.saveplots = 1; 
    
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
    OPTS.nr = 65;
end

% Minimum floe size (m)
if ~isfield(OPTS,'r_p')
    OPTS.r_p = .5; % m
end

% Floe size vector: linearly spaced if not already specified
if ~isfield(FSTD,'Rint')
    
    FSTD.Rint(1) = .5;
    for i = 2:OPTS.nr
        FSTD.Rint(i) = sqrt(2*FSTD.Rint(i-1)^2 - (4/5) * FSTD.Rint(i-1)^2);
    end
    
end

% Increment between floe sizes
if ~isfield(OPTS,'dr')
    % Size Increment
    OPTS.dr = 10*OPTS.r_p; % m
end

%% Parameters of the Ice Thickness Distribution

% Number of Thickness Classes
if ~isfield(OPTS,'nh')
    OPTS.nh = 13; % No. of thickness categories
end

% Smallest Thicknes (m)
if ~isfield(OPTS,'h_p')
    OPTS.h_p = .1; % m
end

% Thickness Increment (m)
if ~isfield(OPTS,'dh')
    OPTS.dh = .2; % m
end

% Vector of Thickness
if ~isfield(FSTD,'H')
    FSTD.H = linspace(OPTS.h_p,OPTS.h_p + (OPTS.nh-1)*OPTS.dh,OPTS.nh); % m
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

% Maximum Floe Size (this does not change in time)
if ~isfield(FSTD,'R_max')
    FSTD.R_max = max(FSTD.Rint) + FSTD.Rint(end) - FSTD.Rint(end-1); %m
end


% Initialize stuff related to FSTD.psi;
Init_Psi;
