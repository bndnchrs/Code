%% Thermodynamics_Init
% This is a script.
% This script initializes the thermodynamic component of the FD model using
% the thermodynamics in Hibler (1979)

if nr > 1
    error('Hibler Thermodynamics does not work with multiple floe categories')
end

fprintf('Using Hibler (79) Thermodynamics \n'); 

if FSTD.DO
    
    if ~isfield(OPTS,'rho_ice')
        % Density of ice (kg/m^3)
        OPTS.rho_ice = 934;
    end
    
    if ~isfield(OPTS,'rho_water')
        % Density of water (kg/m^3)
        OPTS.rho_water = 1000;
    end
    
    if ~isfield(OPTS,'alpha_oc')
        % Ocean Albedo
        OPTS.alpha_oc = .1;
    end
    
    if ~isfield(OPTS,'sigma')
        % Stefan Bolzmann Constant
        OPTS.sigma = 5.67e-8;
    end
    
    
    if ~isfield(OPTS,'alpha_ic')
        % Ice Albedo
        OPTS.alpha_ic = .8;
    end
    
    if ~isfield(OPTS,'cp_water')
        % Specific heat of water (J/kg deg C)
        OPTS.cp_water = 3996;
    end
    
    if ~isfield(OPTS,'L_f')
        % Latent heat of freezing for water (J/kg)
        OPTS.L_f = 334000;
    end
    
    if ~isfield(OPTS,'kice')
        % Thermal conductivity of ice (?/m)
        OPTS.kice = 2.03;
    end
    
    if ~isfield(THERMO,'min_lead_frac')
        % Amount of sea surface guaranteed to be open to solar radiation in
        % melting season... allows for ice to melt
        THERMO.min_lead_frac = .05;
        
    end
    
    if ~isfield(THERMO,'T_ice')
        THERMO.T_ice = -5 + zeros(size([FSTD.H FSTD.H_max])); % Temperature of sea ice
    end
    
    if ~isfield(THERMO,'Toc')
        THERMO.Toc = 0; % The temperature of the ocean as seen by the thermodynamic package
    end
    
    % Whether we allow the thermodynamic loss of thickness to lead to a
    % loss of concentration directly
    
    if ~isfield(THERMO,'allow_adv_loss_H')
        THERMO.allow_adv_loss_H = 1; % By default, we let this happen
    end
    
    if ~isfield(THERMO,'allow_adv_loss_R')
        THERMO.allow_adv_loss_R = 1; % By default, we let this happen
    end
    
    % Floe Merging
    if THERMO.mergefloes
        
        addpath([THERMO.path_of_code 'Packages/Thermodynamics/Merging/')
        Merging_Init;
        
    end    
    
else
    
    % We didn't initialize the model!
    fprintf('NO MAIN MODEL ENABLED... WONT PROCEED FOR THERMODYNAMICS \n ')
    
end
