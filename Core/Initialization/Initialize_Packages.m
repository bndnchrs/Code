function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT,PLOTS] = Initialize_Packages(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT)
%% Initialize_Packages
% This second-level program contains calls to the initialization schemes
% for each process

% Calculate the things like FSTD.conc and FSTD.psi
Init_Psi; 

%% Create a Plot structure
PLOTS = struct; 

%% Initialize the Thermodynamic Model Component

if THERMO.DO
    
    % Add the thermodynamics code to the directory
    
    if isfield(THERMO,'do_Hibler') && THERMO.do_Hibler
        
        % Do we want to do the simple Hibler thermodynamics?        
        addpath([OPTS.path_of_code 'Packages/Thermodynamics/Hibler_79/']);
        % If we do, this points to the directory with its own
        % Thermo_Timestep.m and Thermo_Init.m which relate to the Hibler
        % thermodynamics. 
    else
        
        % If not, we point our path to the typical code
        % Or the H+T (2015) Thermodynamics
        addpath([OPTS.path_of_code 'Packages/Thermodynamics/'])
        
    end
    fprintf('INITIALIZING THERMODYNAMICS \n')
    
    Thermo_Init;
    
end

%% Initialize the Mechanical Mode

if MECH.DO
    
    addpath([OPTS.path_of_code 'Packages/Mechanics/'])
    fprintf('INITIALIZING MECHANICS \n')
    Mech_Init;
    
    %     if isfield(MECH,'do_Thorndike') && MECH.do_Thorndike == 1
    %         addpath('./Thorndike_Mechanics/');
    %         fprintf('THORNDIKE MECHANICS \n')
    %         if nr > 1
    %             error('Thorndike Mechanics does not work with multiple floe categories')
    %         end
    %     end
    %
    
end

%% Initialize the Advective Mode
if ADVECT.DO
    
    addpath([OPTS.path_of_code 'Packages/Advection/']);
    fprintf('INITIALIZING ADVECTION \n')
    Advect_Init; 
    
    
end

%% Initialize the Swell Mode


if WAVES.DO
    
    fprintf('INITIALIZING WAVES FRACTURE \n')
    
    addpath([OPTS.path_of_code 'Packages/Waves/'])
    
    Waves_Init; 
    
end

%% Initialize the Ocean Model


if OCEAN.DO
    
    fprintf('INITIALIZING OCEAN MODEL \n')
    
    addpath([OPTS.path_of_code 'Packages/1D_Ocean/'])
    
    Ocean_Init;
    
end



%% Initialize the Diagnostic Mode

if DIAG.DO
    
    addpath([OPTS.path_of_code '/Packages/Diagnostics/'])
    
    fprintf('INITIALIZING DIAGNOSTICS \n')
    
    Diagnostics_Init;
    
end

end

