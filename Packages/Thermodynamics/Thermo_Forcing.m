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
    
    if ~OCEAN.DO
        
        % These are files that are not determined in a call to Ocean_Forcing.
        
        % If we have no ocean model, the heat flux to a patch of open water is
        % just the LW calculated from a given ocean temperature, and two
        % compute incoming forcing fields.
        
        OCEAN.LW = EXFORC.QLW(FSTD.i);
        OCEAN.SW = EXFORC.QSW(FSTD.i) * (1 - OPTS.alpha_oc);
        OCEAN.LH = 0;
        OCEAN.SH = 0;
        OCEAN.LW_out = -OPTS.sigma * (THERMO.T_oc + 273.14)^4;
        
        
        
        % Positive means warming/melting.
        EXFORC.Q_oc = OCEAN.LW + OCEAN.SW + OCEAN.SH + OCEAN.LH - OCEAN.LW_out;
        
        % This heat flux does not contain the outgoing long-wave heat flux ,
        % which is solved using the Stefan condition at the ice surface
        EXFORC.Q_ic_noLW = OCEAN.LW + OCEAN.SW*(1- OPTS.alpha_ic);
        
    end
    
end