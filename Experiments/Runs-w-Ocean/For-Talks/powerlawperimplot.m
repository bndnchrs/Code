r = linspace(5,1500,300);

alpha = 0:.125:3;
perim = 0*alpha; 

f = zeros(25,300); 
p = f; 

[~, ploc] = find(alpha > 1,1);
ploc = ploc - 1; 

dr = [r(1) diff(r)]; 

for i = 1:length(alpha)
    
    f(i,:) = r.^(-alpha(i));
    
    f(i,:) = f(i,:) / sum(f(i,:)); 

    p(i,:) = 2*f(i,:) ./ (r);
    
    perim(i) = sum(p(i,:).*dr); 
    
end

plot(alpha,100*(perim-perim(ploc))./perim(ploc))
grid on
box on
set(gca,'ydir','normal','layer','top','fontname','helvetica','fontsize',14)
xlabel('\alpha')
ylabel('Floe Perimeter Error (%)')

pos = [5 3];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

saveas(gcf,'floeperim.pdf')