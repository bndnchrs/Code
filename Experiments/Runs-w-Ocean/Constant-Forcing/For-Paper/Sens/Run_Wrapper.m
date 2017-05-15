%% Run_Wrapper
% 12/7/2015 - Chris Horvat

% Define the OPTS structure array that carries the run name and the run
% number

% Empty structure
OPTS = struct(); 

% This path points to where Drive_FSTD.m is contained
OPTS.path_of_code = '/Users/Horvat/Research/FSTD-Code/Code/';
OPTS.savepath = '/Users/Horvat/Research/FSTD-Code/Output/For-Paper/Sens/'; 

% Add to the path
addpath([OPTS.path_of_code 'Core/']); 

% The total number of runs
OPTS.numruns = 7; 
% A list of all the output names of the files we create
OPTS.names = {'zero','p5','one','onep5','two','twop5','three'};

plaws = [0 .5 1 1.5 2 2.5 3]; 

close all

for i = 1:OPTS.numruns
    
    
    OPTS.figpath = ['/Users/Horvat/Research/FSTD-Code/Output/For-Paper/Sens/' OPTS.names{i}];

    % Name of the saved file
    OPTS.NAME = OPTS.names{i}; 
    % Run number
    OPTS.run_number = i;
    OPTS.plaw = plaws(i); 
    % Do the damn thing
    [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT,PLOTS] = Drive_FSTD(OPTS); 
    
    post_plot_data;
    
end


