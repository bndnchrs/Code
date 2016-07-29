function [latSA,Q_open,Q_o,Q_oi,Q_lead,Q_lat,Q_bas] = partition_heat_flux(FSTD,OPTS,THERMO,OCEAN,EXFORC)

% This routine partitions the incoming heat flux Q that is over open water
% into its components,
% some of which contribute to lateral melting, some to frazil formation,
% and some to thickness growth.
% For a heat flux EXFORC.Q_oc, which is the heat exchanged between the
% atmosphere and a unit square meter of ocean, Q_oc is partitioned into
% Q_open - The heat flux impinging directly onto open water
% Q_o - The heat flux that affects open water, which is the net flux at the
% ocean-atmosphere surface - the heat flux from ice to ocean
% Q_lead - The heat flux that affects the development of sea ice

% Q_l is partitioned into:
% Q_lat - the heat flux that leads to the freezing/melting of ice laterally
% Q_bas - the heat flux that leads to the freezing/melting of ice
% vertically

% Additionally, there is a flux Q_oi which is the heat flux between ocean
% and ice

% Q_o = Q_open - Q_oi
% Q_lead = A_l * EXFORC.Q_oc + Q_oi

% First we need to calculate geometric properties of the ice cover

% Ice concentration
conc = integrate_FSTD(FSTD.psi,FSTD.one,FSTD.dA,0);

% This calculates the total lateral surface area of floes in
% contact with water, per unit grid area. The lateral surface area for a given area
% of floes of size r,h is n(r,h) * 2 * pi * h * r. This is equal to f(r,h)
% * 2 * h / r.

% CORRECT CODE
latSA = integrate_FSTD(FSTD.psi,2*FSTD.meshHmid./FSTD.meshRmid,FSTD.dA,0);

% The ratio of the lateral surface are to the total surface area is
% calulated here.

lbrat = (latSA/(latSA + conc));

% If the concentration is zero, to avoid dividing by zero, we set lbrat to be
% one and send all the heat to the sides of floes.
if conc < eps
    
    lbrat = 1;
    
end

% If the concentration is really close to one, we also set it to one.
if abs(conc - 1) < 1e-3
    lbrat = 0;
end

% If the surface area is in fact zero, we set it to be Inf. This is because
% we end up dividing by it later when taking the tendencies.
if latSA == 0;
    latSA = Inf;
end

% For a given open water fraction phi and FSTD, we calculate

% Al - the lead fraction
% Ao - the true "open water" fraction that leads to pancake growth
% subject to Ao + Al = phi
[Al,Ao] = calc_lead_area(FSTD,THERMO);

% Now we need to have a minimal lead fraction at all times, since if the
% concentration is very high we may not be able to melt any ice, even with
% a net heating heat flux. If this isn't occuring, we must adjust

% if Ao+Al < OPTS.min_lead_frac && EXFORC.Q_oc > 0 % Open a little tiny bit at all times if we are heating open water
%
%     openpart = max(0,max(1-conc,OPTS.min_lead_frac)); % Either the open water or the smallest fraction, whichever smaller
%     Al = 1 - conc - openpart; % The rest is lead fraction
%
% end

if OCEAN.DO
        
    % Q_lead is now the heat to the sea ice. It has the negative sign of
    % the heat flux leaving the ocean. 
    Q_lead = -OCEAN.Q_lead;
    % Q_o is the heat flux in the ocean, so has the same sign
    Q_o = OCEAN.Q_o;
    % This is the net heat flux in the ocean region. 
    Q_open = OCEAN.Q_o + OCEAN.Q_lead;
    % This is the net heat flux to the ice from the ocean. 
     
    % This is the lateral component of the floe thermodynamics. Partition
    % the lateral transfer of heat from the surface layer and the vertical
    % transfer of heat from the mixed layer
    Q_lat = (Q_lead + OCEAN.Q_mi) * lbrat;
    
    Q_oi = OCEAN.Q_mi * (1 - lbrat);
   
    % Same, for the basal, which is just what is left over. 
    Q_bas = Q_lead - Q_lat; %  + Q_oi;
    
    
else
    
    
    
    
    Q_open = Ao*EXFORC.Q_oc; % Net heat flux to open water, per square meter of grid.
    % Units are W/m^2.
    
    Q_lead = Al*EXFORC.Q_oc; % Net heat flux to the lead region per square meter of grid.
    % Units are W/m^2
    
    Q_oi = 0;
    
    % Then we subtract to get the net heat flux that affects the region
    % of open water, per square meter of sea ice, as well as the same in
    % the lead region.
    
    Q_o = Q_open; % Take this heat from the ocean
    Q_lead = Q_lead; % Put it into the lead region
    
    
    % Now the total heat flux that goes to the leads is equal to the fraction
    % entering open water that goes into the lead region plus the heat flux
    % exchanged between the ocean and the ice.
    
    
    %%
    if THERMO.fixQ
        
        % Just divide the fixed heat flux into open water and ice components
        
        
        rat = Ao / (Ao + Al);
        
        
        
        if (Ao+Al <= 1e-6)
            % In the case in which there is no open water, so that no heat
            % comes into the ocean in this way, we still need to get the heat
            % in, since we're fixing it. To do this we set the ratio so that
            % all the heat goes to the lead, and then all of that heat goes to
            % the base of ice.
            
            rat = 0;
            lbrat = 0;
            
            
        end
        
        
        Q_o =  rat * EXFORC.Q_oc;
        Q_lead = (1 - rat) * EXFORC.Q_oc;
        
        if ~OCEAN.DO && Q_o > 0
            
            % If we don't have a way to heat the ocean, and we are heating the
            % ice, then all of the heat is used to melt ice
            Q_lead = Q_lead + Q_o;
            Q_o = 0;
            Al = Al + Ao;
            Ao = 0;
            rat = 0;
            
        end
        
        
    end
    
    %%
    
    % The heat flux that goes to the lateral edges of floes is then lbrat
    % multiplied by the total lead heat flux. This is the net heat flux that
    % goes to the sides of ice, per square meter of grid.
    Q_lat = Q_lead * lbrat;
    
    % The rest goes to the floe base. The net heat flux that goes to the base
    % of ice, per square meter of grid.
    Q_bas = Q_lead - Q_lat;
    
    
end

end
