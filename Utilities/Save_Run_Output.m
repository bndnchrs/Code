function Save_Run_Output(FSTD,OPTS,THERMO,MECH,WAVES,DIAG,EXFORC,OCEAN,ADVECT)
% Save_Run_Output
% This is a file that saves each run as its own file. It is called in this
% way to allow it to be nested inside a parfor loop

% Save the files
save([OPTS.savepath OPTS.NAME],'FSTD','OPTS','THERMO','MECH','WAVES','OCEAN','DIAG','EXFORC','ADVECT','-v7.3')

end