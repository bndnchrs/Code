function EXFORC = load_forcing_fields(EXFORC,OPTS,time)
% FIELDS is a cell array that contains the external forcing fields we will
% use to run the model.

% These will all be taken from the 1980-2009 NCEP-2 climatology.

% The input vector time runs from 0 to some number of seconds. We will
% interpolate each field until the duration time.
dummy = 0 * time; 

%% Ocean Stats
EXFORC.UVEL = zeros(2,length(dummy)); 

% UVEL(2) needs to be equal to UVEL(1) unless mechanics is turned on. 
EXFORC.UVEL(1,:) = 1;
EXFORC.UVEL(2,:) = 1;  

%% Set Mechanics Options and External Forcing
EXFORC.nu = .1 * EXFORC.UVEL' / OPTS.Domainwidth; 
EXFORC.nu(:,1) = 0;




end
