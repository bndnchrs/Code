% This routing updates the grids which contain the details of the largest
% floe size as well as the ridging/rafting coefficients, and the volume of
% each floe category. It is updated at each timestep to reflect the changes
% to the largest thickness category.

% Most of this is not needed, we must merely update the last row of meshH.
% Since it doesn't take too much computing time, we will do it each
% timestep anyways. Consider for refining later.

% 2-D matrix of size and thickness coordinates, and volume
[FSTD.meshR,FSTD.meshH] = meshgrid(FSTD.Rint,FSTD.H);
FSTD.meshR = FSTD.meshR';
FSTD.meshH = FSTD.meshH';
FSTD.meshV = pi*FSTD.meshR.^2 .* FSTD.meshH;

%% Define the Staggered Grid

% Define the R vector that lies between subsequent size categories. This is
% the mean floe size in each bin. 
if length(FSTD.Rint) > 1
    % R(1/2) = R(1)/2+R(2)/2
    FSTD.Rmid = .5*(FSTD.Rint(2:end) + FSTD.Rint(1:end-1));
    FSTD.Rmid(end+1) = FSTD.R_max; 
    
else
    
    FSTD.Rmid = FSTD.Rint;
    
end

% Define the H vector that lies between subsequent thickness categories.
% This is the mean thickness in each bin
if length(FSTD.H) > 1
    FSTD.Hmid = .5*([FSTD.H(2:end)] + FSTD.H(1:end-1));
    FSTD.Hmid = [FSTD.Hmid FSTD.H_max]; 
else
    FSTD.Hmid= FSTD.H_max;
end

% The number of such things is equal to the number of intervals (size(FSTD.H)) minus one. This is
% because H_max represents the mean thickness of ice that is larger than
% FSTD.H(end). It is actually the same as FSTD.Hmid(end)

% Find the difference between subsequent size categories
FSTD.dR = diff(FSTD.Rint);
% Define the final value to be equal to the final-1th value.
if length(FSTD.Rint) > 1
    FSTD.dR = [FSTD.dR FSTD.dR(end)];
else
    FSTD.dR = FSTD.Rint;
end

% dH is simply the difference between successive intervals of thickness
FSTD.dH = diff(FSTD.H);
% Since H_max is the mean thickness of the thickest ice, we assume that it
% represents a median, with an interval width equal to twice the difference
% between H_max and max(H).

FSTD.dH_max = 2*(FSTD.H_max - FSTD.H(end)); 

if length(FSTD.H) > 1
    FSTD.dH = [FSTD.dH FSTD.dH_max];
else
    FSTD.dH = FSTD.H;
end

FSTD.dA = bsxfun(@times,FSTD.dH,FSTD.dR'); 

[FSTD.meshRmid,FSTD.meshHmid] = meshgrid(FSTD.Rmid,FSTD.Hmid);

FSTD.meshRmid = FSTD.meshRmid';
FSTD.meshHmid = FSTD.meshHmid';
FSTD.meshVmid = pi*FSTD.meshRmid.^2 .* FSTD.meshHmid;
