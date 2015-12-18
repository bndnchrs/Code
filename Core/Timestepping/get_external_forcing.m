%% get_external_forcing
% This routing loads in the strain rate at each timestep as well as the
% thermodynamic input

if THERMO.DO
    
    % Get the thermodynamic forcing
    Thermo_Forcing;
    
end

if OCEAN.DO
    
    % Get the ocean strain rate from ocean currents (influences ADVECT)
    Ocean_Forcing;
    
end

if MECH.DO

    % Get the strain rate
    get_ice_strain_rate;
    
end    

if WAVES.DO
    
    % Get the sea state
    get_sea_state;
    get_atten_dist; 

end