%% get_thermo_forcing
% This routine obtains the net energy flux to both the open water and
% ice-covered regions. It does not compute the outgoing long-wave over ice,
% however. These fluxes are:

% EXFORC.Q_oc - The net heat flux exchanged between the atmosphere and
% ocean, per square meter of open water.
% EXFORC.Q_ic_noLW - The net heat flux exchanged between the atmosphere and
% ocean and ice, per square meter of sea ice, but without a calculated ice
% long-wave component as this must be calculated via the Stefan condition 
% at the ice surface

if isfield(THERMO,'fixed_Q') && THERMO.fixQ
    % We ask whether the heating at the ocean surface is fixed
    
    % The fixed_Q is a unit of energy we want to apply to the ice/ocean.
    % This is per unit of grid, unlike the other terms if fixed_Q = 0,
    % which are per unit of open water/ice. This is important to remember!
    
    EXFORC.Q_oc = THERMO.fixed_Q(FSTD.i);
    EXFORC.Q_ic_noLW = 0; 
    
else
    % Otherwise we need to specify heating at the ice and ocean surfaces
    
    % These are read-in files in the EXFORC structure
    OCEAN.LW = EXFORC.QLW(FSTD.i);
    OCEAN.SW = EXFORC.QSW(FSTD.i);
    
    % If want to worry about a sensible heat flux
    if THERMO.SHLambda ~=0
        OCEAN.SH = -OPTS.SHLambda * (OCEAN.T + 273.14);
    else
        OCEAN.SH = 0;
    end
    
    % Heat flux above water
    if ~OCEAN.DO
        % If we have no ocean model, the heat flux to a patch of open water is
        % just the sum of all the terms, with long-wave out calculated at the
        % specified temperature THERMO.T_oc.
        
        % Positive means warming/melting. 
        
        EXFORC.Q_oc = OCEAN.LW + OCEAN.SW*(1-OPTS.alpha_oc) - OPTS.sigma * (THERMO.T_oc + 273.14)^4 + OCEAN.SH;
        
    else
        
        % If we do have an ocean model, we specify OCEAN.T
        EXFORC.Q_oc = OCEAN.LW + OCEAN.SW*(1-OPTS.alpha_oc) - OPTS.sigma * (OCEAN.T + 273.14)^4 + OCEAN.SH;
        
    end
    
    % The heat flux at a patch of ice surface is different. This heat flux does not
    % contain the outgoing long-wave heat flux , which is solved using the
    % Stefan condition at the ice surface
    EXFORC.Q_ic_noLW = OCEAN.LW + OCEAN.SW*(1- OPTS.alpha_ic);
    
    % We can also specify a field called fixed_Q. This determines the heating to the ocean
    % and does not allow for re-calculation depending on the ice/ocean state.
    % The flag for this is set to be THERMO.fixQ;
    
end