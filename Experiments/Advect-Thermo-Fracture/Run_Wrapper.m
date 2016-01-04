%% Run_Wrapper
% 12/7/2015 - Chris Horvat

% Define the OPTS structure array that carries the run name and the run
% number

% Empty structure
OPTS = struct(); 

% This path points to where Drive_FSTD.m is contained
OPTS.path_of_code = '/Users/Horvat/Research/FSTD-Code/Code/';
OPTS.savepath = '/Users/Horvat/Research/FSTD-Code/Output/Adv_Therm_Frac/'; 

% Add to the path
addpath([OPTS.path_of_code 'Core/']); 
mkdir(OPTS.savepath); 

% The total number of runs
OPTS.numruns = 1; 
% A list of all the output names of the files we create
OPTS.names = {'Run1'};

for i = 1:OPTS.numruns
    
    % Name of the saved file
    OPTS.NAME = OPTS.names{i}; 
    % Run number
    OPTS.run_number = i;
    % Do the damn thing
    DIAG = Drive_FSTD(OPTS); 
    
end
