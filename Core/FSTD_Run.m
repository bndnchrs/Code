function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_Run(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT)
%% FSTD_Run
% Updated 12/8/2015 - Chris Horvat

% This is the main driver of a single simulation: it uses
% previously initialized conditions and variables to simulate the evolution 
% of the floe size and thickness distribution. 
% If it is not initialized externally (OPTS.driven), it contains a call to
% initialize itself

% struct OPTS
% struct FSTD
% struct THERMO
% struct MECH
% struct WAVES
% struct OCEAN
% struct DIAG

% For all the iterations we want in the model

% The diagnostics at i = 1 are the diagnostics of initialization. To allow
% this, call FSTD_Timestep_Diagnostics here

for iteration_index =  OPTS.start_it: OPTS.end_it
    
    % The index
    FSTD.i = iteration_index;
       
    %
    % Execute one timestep of the model
    if FSTD.DO % Must make sure we want to do this
        
        [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_timestep(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT);
        
    end
       
    
    % At the unit interval save_index, we save the output
    if mod(iteration_index,OPTS.save_index) == 0
        
        disp(iteration_index);
        Save_Run_Output(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT);
        
    end
    
    % If at some point we pick up an error, we quit. 
    
    if FSTD.eflag
        
        return
    
    end
        
    
end

end