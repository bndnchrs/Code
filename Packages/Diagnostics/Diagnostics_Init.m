%% FD_initialize_diagnostics
% This code initializes all the diagnostics arrays that will be used in the
% execution of the FSD code. Whenever a diagnostic is added to this file,
% it must also be added to the script FSTD_Diagnostics.m; 

% Concerns whether we plot out output. Default is no. 
if ~isfield(DIAG,'PLOT_DO')
    DIAG.PLOT_DO = 0; 
end

if FSTD.DO
    
    % Create a structure for the FSTD Diagnostics. 
    DIAG.FSTD = struct(); 
    
        % Diagnostics for general FSTD. We may add as we go. 

        
    %% Multi-dimensional diagnostics
    
    dum_diag = zeros([size(FSTD.psi) OPTS.nt + 1]); % Something as large as psi
    DIAG.FSTD.psi = dum_diag; % The full FSTD
    DIAG.FSTD.diff = dum_diag; % The total change per timestep for all components
    
    %% One-dimensional diagnostics
    dum_diag = zeros(1,OPTS.nt+1);
    
    DIAG.FSTD.time = dum_diag; % The model time
    DIAG.FSTD.numSC = dum_diag; % The total number of sub-intervals in each timestep
    DIAG.FSTD.conc = dum_diag; % The ice concentration at each timestep
    DIAG.FSTD.Rmeanarea = dum_diag; % The mean floe size (area-weighted)
    DIAG.FSTD.Rmeannum = dum_diag; % The mean floe size (number-weighted)
    DIAG.FSTD.Hmean = dum_diag; % The mean ice thickness
    DIAG.FSTD.Vtot = dum_diag; % Total ice volume
    DIAG.FSTD.Hmax = dum_diag; % Maximum Ice Thickness Category
    DIAG.FSTD.Amax = dum_diag; % Maximum Ice Thickness Category Area
    
    
end

if MECH.DO == 1
    % We may add these as we go
    DIAG.MECH = struct(); 
    
end

if THERMO.DO == 1
    % Diagnostics for thermodynamic mode
    % We may add as we go
    
    dum_diag = zeros(length(FSTD.H),OPTS.nt+1);
    
    DIAG.THERMO.dhdt = dum_diag;
    DIAG.THERMO.Tice = dum_diag;
    DIAG.THERMO.Q_cond = dum_diag;
    
    
    dum_diag = zeros(1,OPTS.nt+1);
    
    DIAG.THERMO = struct(); 
    DIAG.THERMO.Q_lead = dum_diag; 
    DIAG.THERMO.Q_lat = dum_diag;
    DIAG.THERMO.Q_o = dum_diag;
    DIAG.THERMO.Q_bas = dum_diag;
    DIAG.THERMO.dV = dum_diag; 
    DIAG.THERMO.drdt = dum_diag; 
    DIAG.THERMO.dc_adv = dum_diag; 
    DIAG.THERMO.dc_pan = dum_diag; 
    DIAG.THERMO.dc_edge = dum_diag; 
    DIAG.THERMO.dc_tot = dum_diag;
end

if WAVES.DO
    % Diagnostics for wave mode
    % We may add as we go
    DIAG.WAVES = struct(); 
    
end

if OCEAN.DO
    % Diagnostics for ocean mode
    % Add as we go
    DIAG.OCEAN = struct(); 
    
end

% The diagnostics at i = 1 are the diagnostics of initialization. To allow
% this, call FSTD_Timestep_Diagnostics here

FSTD_Diagnostics;
