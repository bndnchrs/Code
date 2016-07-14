
function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_timestep(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT)
% Updated 12/8/2015 - Chris Horvat

% This code executes a single predefined timestep of the joint model by
% performing a minimum number of sub-cycles. The way this works is defined
% in a README.pdf file.

% Using the current thicknesses, update gridded size/thickness categories
update_grid;

% Reset counting variables and difference matrices
reset_global_variables;

% Get the grid-level forcing parameters
get_external_forcing;

%% Actual Sub Timestep
while OPTS.dt_sub > 0
    
    %% Reset all variables that are only relevant for one sub-cycle
    reset_local_variables;
    
    %% Get change from advection of ice from out of domain
    
    %% Do pancake growth, ocean heat fluxes, and salt fluxes
    

    if ADVECT.DO
        
        FD_timestep_advect;
        
        
        FSTD.diff = FSTD.diff + ADVECT.diff;
        FSTD.opening = FSTD.opening + ADVECT.opening;
        FSTD.V_max_in = FSTD.V_max_in + ADVECT.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + ADVECT.V_max_out;
        
    end
    
    %% Get Change Due To Mechanics
    
    if MECH.DO && MECH.mag ~= 0
        
        if MECH.do_thorndike
            Thorndike_timestep;
        else
            
            Mech_Timestep;
            
        end
        
        FSTD.diff = FSTD.diff + MECH.diff;
        FSTD.opening = FSTD.opening + MECH.opening;
        FSTD.V_max_in = FSTD.V_max_in + MECH.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + MECH.V_max_out;
        
    end
    
    %% Get Change Due To Thermodynamics
    
    if THERMO.DO
        
        if OCEAN.DO
            % Get the heat fluxes that are appropriate
            OCEAN = Ocean_Fluxes_Petty(OCEAN,OPTS,FSTD);
        
        end
        % This outputs diff_thermo
        
        Thermo_Timestep;
        
        FSTD.diff = FSTD.diff + THERMO.diff;
        FSTD.opening = FSTD.opening + THERMO.opening;
        FSTD.V_max_in = FSTD.V_max_in + THERMO.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + THERMO.V_max_out;
        
    end
    
    %% Do Merging of Floes Thermodynamically
    
    if THERMO.DO && THERMO.mergefloes
        
        FD_merge_floes;
        FSTD.diff = FSTD.diff + THERMO.diff_merge;
        
        % merging preserves area so there is no opening term 
        FSTD.V_max_in = FSTD.V_max_in + THERMO.V_max_in_merge;
       
        
    end
    
    %% Get Change Due to Swell
    
    if WAVES.DO == 1 && EXFORC.stormy(FSTD.i)
        
        Waves_Timestep;
        
        FSTD.diff = FSTD.diff + WAVES.diff;
        FSTD.V_max_in = FSTD.V_max_in + WAVES.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + WAVES.V_max_out;
        
        
        
    end
    
    if OCEAN.DO && THERMO.DO
     
        Ocean_Timestep_Petty; 

    end
        
    
    FSTD.dV_max = FSTD.V_max_in - FSTD.V_max_out;
    
    %% Maximal Timestep
    OPTS.dt_temp = calc_max_timestep(FSTD,OPTS,OCEAN);
    
    %% Update the main variables psi, openwater, and ocean variables
    update_psi;
    
    %% Update those local variables which will change on each timestep
    update_local_variables;
    
    update_grid;
    
    %% Check to make sure the solutions are legal
    check_FD;
    
    %% If we've thrown an error, lets leave this place
    if FSTD.eflag
        return
    end
    


    
end


%% Update variables which are changed on each timestep
update_global_variables;

% Compute diagnostic output
if DIAG.DO
    FSTD_Diagnostics;
end

%% Do Plotting

FD_timestep_plot; 


