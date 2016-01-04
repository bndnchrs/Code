%% Parameters of the Ice Thickness Distribution
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

