%% Parameters of the Ice Thickness Distribution
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
    % Vector of Sizes
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

