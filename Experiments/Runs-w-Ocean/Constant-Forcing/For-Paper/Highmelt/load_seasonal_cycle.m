function EXFORC = load_seasonal_cycle(EXFORC,time)
% FIELDS is a cell array that contains the external forcing fields we will
% use to run the model.

% These will all be taken from the 1980-2009 NCEP-2 climatology.

% The input vector time runs from 0 to some number of seconds. We will
% interpolate each field until the duration time.
dummy = 0 * time; 

EXFORC.QSW = dummy + 290;
EXFORC.QLW = dummy + 270;
EXFORC.TATM = dummy + 0.5;
EXFORC.UATM = dummy + 1;
EXFORC.PRECIP = dummy + .0001 / 86400;
EXFORC.QATM = dummy + 3.6E-03;
EXFORC.PATM = dummy + 101;
EXFORC.EVAP = dummy + 5E-9; 
EXFORC.Hml = dummy + 30;
EXFORC.Hml(end+1) = EXFORC.Hml(end); 
EXFORC.T_b = dummy - 1.8;


end
