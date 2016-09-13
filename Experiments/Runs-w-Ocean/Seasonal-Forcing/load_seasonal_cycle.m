function EXFORC = load_seasonal_cycle(EXFORC,time)
% FIELDS is a cell array that contains the external forcing fields we will
% use to run the model.

% These will all be taken from the 1980-2009 NCEP-2 climatology.

% The input vector time runs from 0 to some number of seconds. We will
% interpolate each field until the duration time.

dt = time(2) - time(1);
smoothdur = 30*86400 / dt; 

ncep_time = (0:11) + .5;
ncep_time = ncep_time' * (365 * 86400/12);

% Get the ncep_climatology parts
FIELDS = pull_ncep_clim;

% Mixed-layer depth
Hml  = cos(2*pi*ncep_time/ (365*86400));
Hml = 0*max(Hml,0);

FIELDS{9} = 25 * (1 + Hml);

% Deep Temperature
Tb  = 1.8*cos(2*pi*(ncep_time/(365*86400) + 1/2));
Tb = max(Tb,0);

FIELDS{10} = -1.8 + Tb;

nyears = ceil(max(time) / (365*86400));

for i = 1:length(FIELDS)
    
    
    
    FIELDS{i} = repmat(FIELDS{i},[nyears 1]);
    
end

ncep_time = (0:nyears*12-1) + .5;
ncep_time = ncep_time' * (365 * 86400/12);

for i = 1:length(FIELDS)
    FIELDS{i} = interp1(ncep_time,FIELDS{i},time,'linear','extrap');

    if i ~= 10 && i~= 9
        
        FIELDS{i} = smooth(smooth(FIELDS{i},2*smoothdur,'sgolay',1),smoothdur);
        
    end
    
    if i == 1 || i == 2
        FIELDS{i}(FIELDS{i} < 0) = 0;
    end
    
end

% For the H_ml field, we need an extra step because it is used to calculate
% w = dH/dt
FIELDS{9}(end+1) = 2*FIELDS{9}(end) - FIELDS{9}(end-1);


%
% Just to recall
% FIELDS{1} = SW;
% FIELDS{2} = LW;
% FIELDS{3} = Temp;
% FIELDS{4} = U;
% FIELDS{5} = Precip;
% FIELDS{6} = Qair;
% FIELDS{7} = Evap

EXFORC.QSW = FIELDS{1};
EXFORC.QLW = FIELDS{2};
EXFORC.TATM = FIELDS{3};
EXFORC.UATM = FIELDS{4};
EXFORC.PRECIP = FIELDS{5};
EXFORC.QATM = FIELDS{6};
EXFORC.PATM = FIELDS{7};
EXFORC.EVAP = FIELDS{8}; 
EXFORC.Hml = FIELDS{9};
EXFORC.T_b = FIELDS{10};


%% Now plot the fields

figure

%%
mos = {'J','','M','','J','','A','','O','','D',''}; 
ticks = .5:11.5
ticks = ticks / 12; 

subplot(241)

plot(time/(365*86400),EXFORC.QLW)
hold on
plot(time/(365*86400),EXFORC.QSW)

set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('Radiation')
xlim([0 1])
ylabel('W/m^2')
set(gca,'xticklabel',{})
ylim([0 max(max(EXFORC.QLW),max(EXFORC.QSW))])
set(gca,'xtick',ticks)
legend('LW','SW','location','south')

subplot(242)
plot(time/(365*86400),EXFORC.UATM,'k')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('2m Wind')
xlim([0 1])
ylabel('m/s')
set(gca,'xticklabel',{})
set(gca,'xtick',ticks)


subplot(243)
plot(time/(365*86400),EXFORC.TATM,'k')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('2 m. Temp')
xlim([0 1])
ylabel('W/m^2')
set(gca,'xticklabel',{})
set(gca,'xtick',ticks)

subplot(244)
plot(time/(365*86400),EXFORC.PRECIP)
hold on
plot(time/(365*86400),EXFORC.EVAP)
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('Precip/Evap')
xlim([0 1])
ylabel('m/s')
set(gca,'xticklabel',{})
set(gca,'xtick',ticks)
legend('P','E','location','south')

subplot(245)
plot(time/(365*86400),EXFORC.QATM,'k')
xlabel('time (years)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('Specific Humidity')
xlim([0 1])
ylabel('g/kg')
set(gca,'xtick',ticks)
set(gca,'xticklabel',mos); 

subplot(246)
plot(time/(365*86400),EXFORC.PATM/101325,'k')
xlabel('time (years)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('Surf Press')
xlim([0 1])
ylabel('atm')
set(gca,'xtick',ticks)
set(gca,'xticklabel',mos); 

subplot(247)
plot(time/(365*86400),EXFORC.T_b,'k')
xlabel('time (years)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('Deep Temp')
xlim([0 1])
ylabel('^\circ C')
set(gca,'xtick',ticks)
set(gca,'xticklabel',mos); 

subplot(248)
plot(time/(365*86400),EXFORC.Hml(1:end-1),'k')
xlabel('time (years)')
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',12)
grid on
box on
title('ML Depth')
xlim([0 1])
ylabel('m')
set(gca,'xtick',ticks)
set(gca,'xticklabel',mos); 

%%


pos = [8.5 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

drawnow

tightfig


saveas(gcf,'ex_forc.fig')
saveas(gcf,'Fig-5.pdf')

end


function FIELDS = pull_ncep_clim
% These are all taken at 80N 0E. Most comes from NCEP, but H comes from
% Peralta-Ferriz

% Shortwave forcing at surface
SW = [   0.0000000E+00
    5.9e-2
    30.68000
    150.8500
    297.1900
    345.7300
    289.4000
    174.4600
    55.94000
    3.100000
    0
    0];

% Longwave downwelling at surface
LW = [    156.5000
    155.3800
    158.5100
    181.3600
    218.9000
    254.0900
    269.3500
    263.5800
    238.0200
    197.6400
    170.0800
    161.0900  ];

% Meridional velocity at 10 m
V = [    1.319000
    1.258000
    1.622000
    1.061000
    1.239000
    1.352000
    0.8830000
    -0.1040000
    0.9050000
    1.879000
    1.014000
    1.412000];

% Zonal velocity at 10 m
U = [    1.253000
    1.193000
    1.217000
    -0.3750000
    -0.7430000
    1.0000000E-03
    0.7050000
    1.358000
    0.9260000
    0.3010000
    0.5090000
    0.5260000];

% Wind speed at 10 m
U = sqrt(smooth(U).^2 + smooth(V).^2);

% 10 m air temp
Temp = [    244.9440
    244.3920
    246.0620
    253.2550
    264.7390
    272.1360
    273.4390
    271.6490
    265.9050
    256.8530
    249.6030
    246.7280];

Temp = Temp - 273.14;

% Surface pressure
PSURF = [101834.7 % in Pa from NCEP-II
    101960.7
    102030.4
    102048.0
    101760.0
    101036.0
    100964.4
    100839.0
    101063.3
    101426.3
    101840.0
    101836.8]/1000; % Convert to kPa 

% Net precipitation
Precip = [   3.4800000E-06
    3.1699999E-06
    3.3300000E-06
    4.0099999E-06
    5.9600002E-06
    8.1400003E-06
    8.8099996E-06
    1.1090000E-05
    8.1099997E-06
    5.8900000E-06
    3.3700001E-06
    3.6199999E-06] / 999.8; % in kg / m^2 s - convert to m/s

% 2-m specific humidity

Qair = [   4.0530000E-04
    3.7699999E-04
    4.3960000E-04
    8.2620000E-04
    1.9898000E-03
    3.3108999E-03
    3.6126000E-03
    3.2358000E-03
    2.2052000E-03
    1.1133000E-03
    5.8930001E-04
    4.5759999E-04];

Evap = [    1.370000    
    1.300000    
    1.500000    
    3.130000    
    7.330000    
    13.20000    
    12.14000    
    7.250000    
    1.910000    
   0.7300000    
    1.480000    
    1.740000  ]; 

Evap = Evap / (999.8 * 2.501 * 10^6); 

FIELDS{1} = SW;
FIELDS{2} = LW;
FIELDS{3} = Temp;
FIELDS{4} = U;
FIELDS{5} = Precip;
FIELDS{6} = Qair;
FIELDS{7} = PSURF;
FIELDS{8} = Evap; 

end
