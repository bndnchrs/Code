%% check_ridge_out
% does not allow us to continue if any of our interactions lead to either a
% loss in Volume or a gain in area. 


for r1 = 1:length(R)
    for r2 = 1:length(R)
        
        for h1 = 1:length(H)
            for h2 = 1:length(H)
                
                diagone = .5;
                diagtwo = 1;
                
                if r1 == r2 && h1 == h2
                    diagone = 1;
                    diagtwo = 2;   
                end
              
                Ain = diagone*Kfac(r1,r2,h1,h2)*R(S_R(r1,r2))^2; 
                Aout = diagtwo*(R(r1)^2 + R(r2)^2); 
                
                if Ain - Aout > 0
                    
                    disp('Areas will increase. All hope is lost.')
                    r1
                    r2
                    h1
                    h2
                    error('PROBLEM')
                end
                

                
            end
        end
    end
end