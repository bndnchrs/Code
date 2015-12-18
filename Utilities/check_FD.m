%% FD_Check
% This routine checks for typical errors in the code. 
% It also outputs required or
% requested information to the command line.

if FSTD.H_max < 0
    disp('Negative H')
    FSTD.eflag = 1; 

end

 if sum(FSTD.psi(:)) > 1 + 1e-8
     sum(FSTD.psi(:))
     disp('Too much conc')
     FSTD.i
     FSTD.eflag = 1; 
 end

if isnan(FSTD.psi)
    disp('NaNned out')
    FSTD.i
    FSTD.eflag = 1; 
end

if abs(FSTD.opening + sum(FSTD.diff(:))) > eps
    disp('Bad Opening/Closing')
    FSTD.eflag = 1; 
end

if min(FSTD.psi(:)) < 0
    disp(FSTD.i)
    disp('Less Than Zero after cutting')
    FSTD.eflag = 1; 
end
