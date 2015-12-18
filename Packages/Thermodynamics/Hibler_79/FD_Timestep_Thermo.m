%% Hibler.m 
% This code simulates the thermodynamic sea ice code of Hibler '79 which is
% valid for a two-category ice, with thicknes and an area

% The equations therin are
% dh/dt = f(h/A) A + (1-A)f(0) 
% dA/dt = f(0)/h_0 * (1-A) if f(0) > 0
%         0 if f(0) < 0 
%       +
%         0 if dh/dt > 0
%         A/2h dh/dt if dh/dt < 0

% Here we will use as f(0) the growth rate of pancake ice. 

hzero = .4; % meters

%% FD_timestep_thermo
% This routing calculates the tendency due to one thermodynamic timestep

% Calculate total floe surface area per sea surface area

SAmean = integrate_FD(psi,2*bsxfun(@rdivide,[H H_max],R'),1);

if SAmean == 0;
    SAmean = Inf;
end

% Determine what portion of the heat flux goes where
partition_heat_flux;


if Q > 0 % Some small area is open to heat fom the atmosphere
    
    Al = max(Al,smallfrac);
    Q_l = Q*Al;
    Q_lat = Q_l * (SAmean/(SAmean + conc));
    
    if conc < eps
        
        Q_lat = Q_l;
        
    end
    
    if SAmean == Inf
        
        Q_lat = 0;
        
    end
    
    Q_vert = Q - Q_o - Q_l; % Heat Flux to ice top
    Q_bas = Q_l - Q_lat; 

end

%% Pancake Growth. This is the dA/dt from the Hibler Model. 

fzero = -(Q_o + Q_lat) / (L_f * rho_ice * hzero) ;

% Very simple vertical heat budget
Ttop = -2; 
deltaT = Ttop - T_ocean;

dhdt_surf = -.1*Q*conc/(rho_ice*L_f);
dhdt_bas = -Q_bas/(rho_ice*L_f);

if conc == 0
    dhdt_bas = 0;
end

dhdt_diff = conc*kice*deltaT./(rho_ice*L_f*[H_max]);
dhdt_diff = 0; 

dhdt = (dhdt_surf + dhdt_bas + dhdt_diff);

dhdt = dhdt * conc + (1-conc)*fzero; 

if fzero*(1-conc)*dt_sub + H_max < 0
    dhdt = conc*(dhdt_surf + dhdt_bas + dhdt_diff);
end

RHS = 0; 

if fzero >= 0 
    RHS = fzero * (1-conc); 
end

if dhdt < 0
    
    RHS = RHS + conc/(2*H_max) * dhdt; 
end
    
dAdt = RHS; 



%% Thickness Growth Rate

%% Ocean Temperature Changes
% Ocean Temperature changes due to direct heat fluxes into the ocean and to
% the diffusive heat flux from the ice melting. 

Q_into_ocean = Q_o * (Q_o > 0); 

dTdt_ocean = (Q_into_ocean - sum(ITD.*dhdt_diff)*(rho_ice*L_f))/ (cp_water*rho_water);

%% Now Put It Back Into the Common Framework

dhdt = dhdt; 



if sum(psi(:)) == 0 && dAdt < 0
    dAdt = 0; 
    dhdt = 0; 
end

diff_thermo = dAdt;




opening_thermo = -sum(diff_thermo(:));

T_ocean = T_ocean + dt_sub*dTdt_ocean;

V_max_in_thermo = A_max * dhdt; 

if dAdt > 0
    % If freezing, this is at 40 cm
    V_max_in_thermo = V_max_in_thermo + hzero*dAdt; 
else
    % If melting, this comes off the regular ice
    V_max_in_thermo = V_max_in_thermo + H_max*dAdt; 
end


if V_max == 0 && H_max ~= 0 && V_max_in_thermo < 0
    
    V_max_in_thermo = 0; 
    
end
    
V_max_out_thermo = 0;


%% For Diagnostics Purposes
edgegrowth = -Inf;
drdt = -Inf; 
pancakes = -Inf; 

DADTSAVE_HIBLER(i) = dAdt; 
DHDTSAVE_HIBLER(i) = dhdt; 
