%% FD_initialize_diagnostics
% This code initializes all the diagnostics arrays that will be used in the
% execution of the FSD code. Whenever a diagnostic is added to this file,
% it must also be added to the script FSTD_Diagnostics.m;

% Concerns whether we plot out output. Default is no.
if ~isfield(DIAG,'DO_PLOT_FSTD')
    DIAG.DO_PLOT_FSTD = 0;
end

if ~isfield(DIAG,'DO_PLOT_OCEAN')
    DIAG.DO_PLOT_OCEAN = 0;
end

if ~isfield('DIAG','plot_realtime')
    DIAG.plot_realtime = 1;
end

if FSTD.DO
    
    % Create a structure for the FSTD Diagnostics.
    DIAG.FSTD = struct();
    
    % Diagnostics for general FSTD. We may add as we go.
    DIAG.FSTD.R = FSTD.Rint;
    DIAG.FSTD.H = FSTD.H;
    
    %% Multi-dimensional diagnostics
    
    dum_diag = zeros([size(FSTD.psi) OPTS.nt + 1]); % Something as large as psi
    DIAG.FSTD.psi = dum_diag; % The full FSTD
    DIAG.FSTD.dA = dum_diag; % The DA matrix
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
    
    DIAG.FSTD.resid_adjust = dum_diag; 
    
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
    DIAG.THERMO.Q_open = dum_diag;
    DIAG.THERMO.Q_bas = dum_diag;
    DIAG.THERMO.dV = dum_diag;
    DIAG.THERMO.drdt = dum_diag;
    DIAG.THERMO.dc_adv = dum_diag;
    DIAG.THERMO.dc_pan = dum_diag;
    DIAG.THERMO.dc_edge = dum_diag;
    DIAG.THERMO.dc_tot = dum_diag;
    DIAG.THERMO.diffnet = dum_diag;
    
end

if THERMO.DO && THERMO.mergefloes
    
    DIAG.THERMO.alphamerge = dum_diag;
    DIAG.THERMO.mn1 = dum_diag;
    DIAG.THERMO.diffmerge = dum_diag; 
    
end

if THERMO.DO && THERMO.dosemtner
    DIAG.THERMO.QSW = dum_diag; 
    dum_diag = zeros(length(FSTD.H),OPTS.nt+1);
    DIAG.THERMO.QLW = dum_diag; 
            
end

if WAVES.DO
    % Diagnostics for wave mode
    % We may add as we go
    DIAG.WAVES = struct();
    
    dum_diag = zeros(1,OPTS.nt+1);
    DIAG.WAVES.tau = dum_diag;
    DIAG.WAVES.diffnet = dum_diag;
    
    dum_diag = zeros([size(FSTD.psi) OPTS.nt+1]);
    DIAG.WAVES.In = dum_diag;
    DIAG.WAVES.Out = dum_diag;
    
end

if MECH.DO
    DIAG.MECH = struct();
    
    dum_diag = zeros(1,OPTS.nt+1);
    DIAG.MECH.mag = dum_diag;
    DIAG.MECH.epsI = dum_diag;
    DIAG.MECH.epsII = dum_diag;
    DIAG.MECH.diffnet = dum_diag;
    
end

if ADVECT.DO
    DIAG.ADVECT = struct();
    
    dum_diag = zeros(1,OPTS.nt+1);
    DIAG.ADVECT.diffnet = dum_diag;
    DIAG.ADVECT.FSTD_in = ADVECT.FSTD_in;
    DIAG.ADVECT.dc_adv = dum_diag;
    
end

if OCEAN.DO
    % Diagnostics for ocean mode
    % Add as we go
    DIAG.OCEAN = struct();
    
    dum_diag = zeros(1,OPTS.nt+1);
    DIAG.OCEAN.T = dum_diag;
    DIAG.OCEAN.S = dum_diag;    
    DIAG.OCEAN.QLH = dum_diag; 
    DIAG.OCEAN.QSH = dum_diag; 
    DIAG.OCEAN.QLW_out = dum_diag; 
    DIAG.OCEAN.QSW_in = dum_diag; 
    DIAG.OCEAN.QLW_in = dum_diag;     
    DIAG.OCEAN.Qlead = dum_diag; 
    DIAG.OCEAN.Qsurf_at = dum_diag; 
    DIAG.OCEAN.Qsurf_ml = dum_diag; 
    DIAG.OCEAN.Q_mi = dum_diag;     
    DIAG.OCEAN.T_s = dum_diag; 
    DIAG.OCEAN.H_ml = dum_diag; 

end

% The diagnostics at i = 1 are the diagnostics of initialization. To allow
% this, call FSTD_Timestep_Diagnostics here

FSTD_Diagnostics;
