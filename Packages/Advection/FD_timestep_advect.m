%% FD_timestep_advect
% We think of the advective velocities as being u(x) + dudx * L
% Therefore we consider an advective flux at each end of the domain
% If there is convergence, this means there is an additional flux at the
% left side
% If there is divergence, the same is true at the right side

% Shear does not lead to advective changes in the ice concentration, though
% it will do so when there is mechanical forcing.

% This code is very simple. It requires v2, the velocity of the ice at the
% eastern boundary, and v1, the velocity of the ice at the western
% boundary. Additionally, the FSTD of ice at the western boundary is
% required. It is assumed that the eastern boundary.

% The outgoing advective flux is equal to the velocity of the ice,
% multiplied by the FSTD divided by the domain width (to get units of per
% unit time)
ADVECT.out = (1/OPTS.Domainwidth) * ADVECT.v2 * FSTD.psi;
% if MECH is on, ADVECT.v2 is zero. 

% If it is off, we have no way to stop advecting things in, so we have an
% advective velocity out that depends on the local conc
ADVECT.in = (1/OPTS.Domainwidth) * (ADVECT.v1 * ADVECT.FSTD_in);

% Divide by domainwidth because this is a fractional area loss
ADVECT.diff = ADVECT.in - ADVECT.out;

% Really just V_max_delta
ADVECT.V_max_in = FSTD.H_max_i*sum(ADVECT.in(:,end).*FSTD.dA(:,end));
ADVECT.V_max_out = FSTD.H_max*sum(ADVECT.out(:,end).*FSTD.dA(:,end));

% Pretty simple
ADVECT.opening = -sum(ADVECT.diff(:).*FSTD.dA(:));

