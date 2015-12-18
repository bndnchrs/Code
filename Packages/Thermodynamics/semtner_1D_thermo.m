function [dhdt_surf,dhdt_base,F_s,F_a,T_new] = semtner_1D_thermo(Qtop,Qbase,h,T_0)
%% function [dhdt,T] = semtner_1D_thermo(Q,h)
% Q is the surface heat flux that is not the long-wave part (positive =
% melting)
% Qbase is the heat flux from the ocean (positive = melting)
% h is the vector of ice thicknesses
% T_0 is the vector of ice temperature per thickness

% dhdt is the total time rate of change of ice thickness in each thickness
% category
% F_s is the conductive flux per thickness category (positive = melting)
% F_a is the surface heat budget (0 if not melting). 
% T is the temperature of the ice per thickness category

rho_0 = 917; %kg/m^3 density of ice
L_0 = 334000; %J/kg latent heat of freezing of ice
kappa_I = 2.034; %J/(s deg C m)
sigma = 5.67e-8; % Boltzmann constant
T_B = -2; % Degrees C (temp of bottom)

F_a = zeros(1,length(h)); 
F_s = F_a; 
dhdt_surf = F_a; 

q_s = rho_0 * L_0 ; % Enthalpy of freezing/melting of surface ice
qbase = q_s * .9; 

%% For each thickness 
for i = 1:length(h)
  
    %%
    H = h(i);
    T_ice = T_0(i);
    
    % The conductive heat flux from the upper surface to the lower
    % > 0 implies cooling of surface (rare, or when near melting)
    % < 0 implies warming of surface (otherwise)
    cond_HF = @(T_i) kappa_I * (T_i - T_B)./H;
    
    % The surface heat flux from the balance of fluxes - conductive
    surf_HF = @(T_i) Qtop - cond_HF(T_i) - sigma*(T_i+273.14).^4 ;
    
    % Find the temperature which equalizes the conductive flux and the
    % surface heat balance, near T_0.     
    T_new(i) = fzero(surf_HF,T_ice); % New ice temperature
 %%   
    % If the temperature is larger than the melting point of fresh ice
    if T_new(i) > 0
       %% 
        % Set the temperature to 0
        T_new(i) = 0;
        % Calculate the budget at 0 temperature
        F_a(i) = surf_HF(T_new(i)); % >0 when not balanced (T_new = 0).
        % Includes the conductive heat flux
        F_s(i) = cond_HF(T_new(i)); % < 0 almost always.
        
        % The melting is now done so as to balance the heat budget
        % The total heat input is the residual F_a
        dhdt_surf(i) = -F_a(i)/q_s; % - F_a(i) / q_s; 
        
    else
        
        dhdt_surf(i) = 0;
        
        % F_s is the conductive heat flux 
        F_s(i) = cond_HF(T_new(i));
        
    end
    
end
    
%Qbase is positive for melting
% Lose thickness due to input from water to ice
% Lose thickess due to heat flux from ice surface to ice base
dhdt_base = -(Qbase/qbase + F_s / qbase); 
    
    
