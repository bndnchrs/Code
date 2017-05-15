function create_movie_FD(file_string,startup,size_plot)

if strcmp(file_string,'close')
    %% We shut it down
    MakeQTMovie finish
    
else % we make the movie
    str = [file_string '.qt'];
    
    if nargin < 2
        startup = 1;
        size_plot = [0 0 24 12];
    else
        if nargin < 3
            size_plot = [0 0 24 12];
        end
    end
    
    %% If we are on our first iteration
    if startup == 1
        
        disp(['Opened the movie file: ' str]);
        MakeQTMovie('start',str);
        MakeQTMovie('size',768 * size_plot(3:4)/size_plot(3));
        MakeQTMovie('framerate',7); 
        set(gcf,'windowstyle','normal','position',size_plot,'paperposition',size_plot,'papersize',size_plot(3:4),'units','inches');
        set(gcf,'windowstyle','normal','position',size_plot,'paperposition',size_plot,'papersize',size_plot(3:4),'units','inches');
 %       saveas(gcf,[file_string 'tmp.jpg'])
        print('-djpeg',[file_string 'tmp.jpg'])
        MakeQTMovie('addimage',[file_string 'tmp.jpg'])
        
        % If we are not
    else if startup > 1
            set(gcf, 'Units', 'inches');
            set(gcf, 'Position', size_plot); %% [left, bottom, width, height];
            print('-djpeg',[file_string 'tmp.jpg'])
            MakeQTMovie('addimage',[file_string 'tmp.jpg'])
            
        end
        
    end
    
end