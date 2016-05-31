
function [FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT] = FSTD_timestep(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT)
% Updated 12/8/2015 - Chris Horvat

% This code executes a single predefined timestep of the joint model by
% performing a minimum number of sub-cycles. The way this works is defined
% in a README.pdf file.

% Using the current thicknesses, update gridded size/thickness categories
update_grid;

% Reset counting variables and difference matrices
reset_global_variables;

% Get the grid-level forcing parameters
get_external_forcing;


%% Actual Sub Timestep
while OPTS.dt_sub > 0
    
    %% Reset all variables that are only relevant for one sub-cycle
    reset_local_variables;
    
    %% Get change from advection of ice from out of domain
    
    if ADVECT.DO
        
        FD_timestep_advect;
        
        
        FSTD.diff = FSTD.diff + ADVECT.diff;
        FSTD.opening = FSTD.opening + ADVECT.opening;
        FSTD.V_max_in = FSTD.V_max_in + ADVECT.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + ADVECT.V_max_out;
        
    end
    
    %% Get Change Due To Mechanics
    
    if MECH.DO && MECH.mag ~= 0
        
        if MECH.do_thorndike
            Thorndike_timestep;
        else
            
            Mech_Timestep;
            
        end
        
        FSTD.diff = FSTD.diff + MECH.diff;
        FSTD.opening = FSTD.opening + MECH.opening;
        FSTD.V_max_in = FSTD.V_max_in + MECH.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + MECH.V_max_out;
        
    end
    
    %% Get Change Due To Thermodynamics
    
    if THERMO.DO
        
        % This outputs diff_thermo
        
        Thermo_Timestep;
        
        FSTD.diff = FSTD.diff + THERMO.diff;
        FSTD.opening = FSTD.opening + THERMO.opening;
        FSTD.V_max_in = FSTD.V_max_in + THERMO.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + THERMO.V_max_out;
        
    end
    
    %% Do Merging of Floes Thermodynamically
    
    if THERMO.DO && THERMO.mergefloes
        
        FD_merge_floes;
        FSTD.diff = FSTD.diff + THERMO.diff_merge;
        
        % merging preserves area so there is no opening term
        
        FSTD.V_max_in = FSTD.V_max_in + THERMO.V_max_in_merge;
        
    end
    
    %% Get Change Due to Swell
    
    if WAVES.DO == 1 && EXFORC.stormy(FSTD.i)
        
        Waves_Timestep;
        
        FSTD.diff = FSTD.diff + WAVES.diff;
        FSTD.V_max_in = FSTD.V_max_in + WAVES.V_max_in;
        FSTD.V_max_out = FSTD.V_max_out + WAVES.V_max_out;
        
        
        
    end
    
    %% Do pancake growth, ocean heat fluxes, and salt fluxes
    
    if OCEAN.DO && THERMO.DO
        
        Ocean_Timestep;
        
        % From pancake growth
        % We keep the same timestep here: will this be a problem later?
        FSTD.diff = FSTD.diff + OCEAN.diff;
        FSTD.opening = FSTD.opening + OCEAN.opening;
        FSTD.V_max_in = FSTD.V_max_in + OCEAN.V_max_in;
        
    end
    
    FSTD.dV_max = FSTD.V_max_in - FSTD.V_max_out;
    
    %% Maximal Timestep
    % Calculate the maximal timestep possible so that the FSTD is never
    % smaller than zero anywhere.
    OPTS.dt_temp = calc_max_timestep(FSTD.psi,FSTD.diff,OPTS.dt_sub,0,FSTD.dA,1);
    % Now calculate the maximal timestep so that the volume in the highest
    % thickness category is >= 0.
    OPTS.dt_temp = calc_max_timestep(FSTD.V_max,FSTD.dV_max,OPTS.dt_temp,1,FSTD.dA(end,end),1);
    
    % If there is an error, dt_temp comes back as a string. We then error
    % ourselves out.
    if strcmp(OPTS.dt_temp,'dt')
        FSTD.eflag = 1;
        fprintf('Cut timestep is negative at timestep %d',FSTD.i);
    end
    
    
    %% Update the main variables psi, openwater, and ocean variables
    update_psi;
    
    %% Update those local variables which will change on each timestep
    
    update_local_variables;
    
    update_grid;
    
    %% Check to make sure the solutions are legal
    check_FD;
    
    %% If we've thrown an error, lets leave this place
    if FSTD.eflag
        return
    end
    
end


%% Update variables which are changed on each timestep
update_global_variables;

% Compute diagnostic output
if DIAG.DO
    FSTD_Diagnostics;
end

if mod(FSTD.i,1) == 0
    
    %% Do Plotting Here
    
    cols = [228,26,28
        55,126,184
        77,175,74
        152,78,163
        255,127,0]/256;
    
    %     subplot(121)
    %     pcolor(FSTD.Rint,FSTD.H,log10(FSTD.psi.*FSTD.dA)');
    %     set(gca,'clim',[-4 0])
    %     shading interp
    %     grid on
    %     box on
    %     set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
    %     xlabel('Size')
    %     ylabel('Thickness')
    %
    %     subplot(122)
    
    
    subplot(121)
    hold on
    plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Vtot(1:FSTD.i)/DIAG.FSTD.Vtot(1),'color',cols(1,:),'linewidth',2)
    plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.conc(1:FSTD.i)/DIAG.FSTD.conc(1),'color',cols(2,:),'linewidth',2)
    % Thickness
    plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Hmean(1:FSTD.i)/DIAG.FSTD.Hmean(1),'color',cols(3,:),'linewidth',2)
    hold off
    ylim([.75 1.25])
    xlim([0 max(DIAG.FSTD.time(FSTD.i)/86400,1)])
    
    legend('Volume','Conc','Thick')
    grid on
    box on
    xlabel('Time (days)')
    ylabel('Frac. of initial')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
    
    
    subplot(122)
    plot(DIAG.FSTD.time(1:FSTD.i)/86400,DIAG.FSTD.Rmeanarea(1:FSTD.i),'color',cols(4,:))
    hold on
    [a,b] = max(WAVES.spec);
    
    % This is R which is half the wavelength of the peak period
    R_longterm = FSTD.Rmid(b);
    % This is the initial mean floe size
    R_init = FSTD.Rmid(end);
    
    R_pred = R_longterm + exp(1).^(-DIAG.FSTD.time(1:FSTD.i)/WAVES.tau)*(R_init - R_longterm);
    plot(DIAG.FSTD.time(1:FSTD.i)/86400,R_pred,'--','color',cols(3,:));
    
    plot(DIAG.FSTD.time(1:FSTD.i)/86400,0*DIAG.FSTD.time(1:FSTD.i) + FSTD.Rmid(b),'--k');
    
    hold off
    xlim([0 max(DIAG.FSTD.time(FSTD.i)/86400,1)])
    ylim([0 max(FSTD.Rmid)])
    legend('Floe Size','Half of Peak Period')
    grid on
    box on
    xlabel('Time (days)')
    ylabel('Floe Size')
    set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
    pos = [8 4];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches');
    
end

drawnow
