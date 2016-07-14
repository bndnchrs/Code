%% Thermodynamics_Init
% This is a script.
% This script initializes the thermodynamic component of the FD model using
% the thermodynamics in Horvat and Tziperman (2015).

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
    
    if ~isfield(THERMO,'lead_width')
        % The width of the lead region
        THERMO.lead_width = OPTS.r_p; 
        
    end
    
    if ~isfield(OPTS,'min_lead_frac')
        % Amount of sea surface guaranteed to be open to solar radiation in
        % melting season... allows for ice to melt
        OPTS.min_lead_frac = .05;
        
    end
    
    if ~isfield(THERMO,'panloc_r')
        % This is the size index R(panloc_r) into which pancake floes will
        % be created when they form. 
        THERMO.panloc_r = 1; 
    
    end
    
    if ~isfield(THERMO,'panloc_h')
        % This is the thickness index H(panloc_h) into which pancake floes will
        % be created when they form. 
        THERMO.panloc_h = 1; 
    
    end
    
    
    if ~isfield(THERMO,'T_ice')
        THERMO.T_ice = -5 + zeros(size(FSTD.Hmid)); % Temperature of sea ice
    end
    
    if ~isfield(THERMO,'Toc')
        THERMO.Toc = 0; % The temperature of the ocean
    end
    
    if ~isfield(THERMO,'fixQ') || THERMO.fixQ == 0
        
        THERMO.fixQ = 0; % A flag for whether we want to fix the external heating to the ocean surface
        % When on, we do not allow the ice to radiate long-wave heating and
        % only use THERMO.fixed_Q to heat the ocean
        
        % This means we need external forcing fields QLW and QSW. If
        % they don't exist, specify them as fields of zeros
        if ~isfield(EXFORC,'QLW')
            
            EXFORC.QLW = 0*(1:OPTS.nt);
        end
        
        if ~isfield(EXFORC,'QSW')
            EXFORC.QSW = 0*(1:OPTS.nt);
        end
        
        % We also need an ocean temperature
        if ~OCEAN.DO && ~isfield(THERMO,'T_oc')
            THERMO.T_oc = 0; % Degrees C
        end
        
    else
        
        THERMO.dosemtner = 0; % Don't do the semtner thermo since we don't want to have a temperature of ice
        
        if ~isfield(THERMO,'fixed_Q')
            fprintf('Specified a fixed heating but did not give the field. Quitting \n')
            error('')
        end
        
    end
    
    
    % Whether we allow the thermodynamic loss of thickness to lead to a
    % loss of concentration directly
    
    if ~isfield(THERMO,'allow_adv_loss_H')
        % This concerns letting ice concentration be removed through
        % advection from the thinnest thickness category. Since
        % thickness==volume, in order to conserve volume we need to do
        % this. 
        THERMO.allow_adv_loss_H = 1; % By default, we let this happen. 
    end
    
    if ~isfield(THERMO,'allow_adv_loss_R')
        % This concerns letting ice concentration be removed through
        % advection from the thinnest size category. This would result in
        % an additional loss of sea ice volume from the ocean that we
        % cannot account for. Therefore we do not allow this to happen. 
        THERMO.allow_adv_loss_R = 0; 
    end
    
    % Floe Merging
    if ~isfield(THERMO,'mergefloes')
        
        THERMO.mergefloes = 0; % By default, we don't merge floes
        
    end
    
    % If we do allow merging, run the code.
    if THERMO.mergefloes
        
        addpath([OPTS.path_of_code 'Packages/Thermodynamics/Merging/'])
        Merging_Init;
        
    end
    
    % This is the choice of how to handle the vertical thermodynamics. If
    % dosemtner is on, we use the semtner 2001 thermo. 
    if ~isfield(THERMO,'dosemtner')
        THERMO.dosemtner = 1; 
    end
    
    
    
else
    
    % We didn't initialize the model!
    fprintf('NO MAIN MODEL ENABLED... WONT PROCEED FOR THERMODYNAMICS \n ')
    
end
