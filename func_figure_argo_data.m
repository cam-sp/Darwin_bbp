function func_figure_argo_data(argo,idx_argo_biomes)

addpath C:\Users\Camila\Desktop\Backup\Projects\3Doutputs\cmap_ocean

figure
tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact')
idxtile=1:3:12;

cmap=brewermap(100,'Blues');
% colormap(cmap([30,80,100],:))

gamma=1;
% wl_correct_argo=1.5e-3.*(700./443).^(-gamma);


bbp_th=1e-3;
wl_correct_argo=bbp_th.*(700./443).^(-gamma);
    
% figure
% tiledlayout(2,2)
monv=[0,3:3:12];    
for i=1:4
    nexttile(idxtile(i))
    ax = axesm ( 'Origin',  [0 -90 0], 'MapProjection','pcarree', 'Frame', 'on',...
            'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);    
    ax.XColor = 'white';
    ax.YColor = 'white';
    axis tight manual
    
    idx=argo.month>monv(i) & argo.month<=monv(i+1) & argo.idx_used ;
    scatterm(argo.lat(idx),argo.lon(idx),2,'filled','markerfacecolor',[235, 177, 42]./255,'markerfacealpha',0.5); %,[227, 172, 52]./255 cmap(40,:)
    hold on
    
    idx=argo.bbp<wl_correct_argo & argo.month>monv(i) & argo.month<=monv(i+1);
    scatterm(argo.lat(idx),argo.lon(idx),2,'filled','markerfacecolor',cmap(80,:),'markerfacealpha',0.5);
    
    geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
%     box off
    framem('FlineWidth',0.7,'FEdgeColor','k')
    set(gcf,'color','w');
    if i==1
        ttl=title('a. January-March');
    elseif i==2
        ttl=title('b. April-June');
    elseif i==3
        ttl=title('c. July-September')  ;      
    elseif i==4
        ttl=title('d. October-December') ;       
    end
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left';
end

cmap=brewermap(100,'Blues');
colormap(cmap([30,80,100],:))


edges=-4:0.05:-2;
cedges=edges(1:end-1)+0.025;   
pans=[2,5,6,8,9,11,12];
st=1;
cmap=flip(brewermap(12,'spectral'));
cmap=[cmap(1:6,:);flip(cmap(1:6,:))];

bbcrit=10.^-3;
bbcrit2=10.^-3.2;
% bbp_argo_470=bbp_argo.*(470/700)^(-1);
bbpcrit_700=bbcrit.*(700/470)^(-1);
bbpcrit2_700=bbcrit2.*(700/470)^(-1);

for i=1:4
    if i==1
        idx=idx_argo_biomes.trop;
    elseif i==2
        idx=idx_argo_biomes.oligo;
    elseif i==3
        idx=idx_argo_biomes.temp;
    elseif i==4
        idx=idx_argo_biomes.polar;        
    end
        for j=1:2
            if i==1 && j==1
                    nexttile(2)
                    idx_lat=argo.lat>=-100;       
                    st=st+1;
            elseif i>1
                if j==1
                    nexttile(pans(st))
                    idx_lat=argo.lat>=0;
                else
                    nexttile(pans(st))
                    idx_lat=argo.lat<0;
                end
                    st=st+1;
            end
            for imon=1:12
                histc=histcounts(log10(argo.bbp(idx & idx_lat & argo.month==imon)),edges);
                histc=histc./sum(histc);            
                h1=plot(cedges,cumsum(histc).*100,'color',cmap(imon,:));
                hold on
                h2=plot(cedges.*0+log10(bbpcrit_700),linspace(0,100,length(histc)),'k:','linewidth',1);
                h3=plot(cedges.*0+log10(bbpcrit2_700),linspace(0,100,length(histc)),'k--','linewidth',0.8);
            end
            ylim([0 100])
            yticks([0:25:100])
            
                if i==1
                    ttl=title('e. Tropical');
                    ylabel('% of data')
%                     legend([h1,h2,h3])
                    
                elseif i==2
                    if j==1
                        ttl=title('f. Oligotrophic N');
                    else
                        ttl=title('i. Oligotrophic S');    
                    end
                    if j==1
                       ylabel('% of data') 
                    end
                elseif i==3
                    if j==1
                        ttl=title('g. Temperate N');
                    else
                        ttl=title('j. Temperate S');    
                    end
                    if j==1
                       ylabel('% of data') 
                    end
                elseif i==4
                    if j==1
                        ttl=title('h. Polar N');
                    else
                        ttl=title('k. Polar S');    
                    end           
                    xlabel('Log_{10}b_{bp}(700)')
                    if j==1
                       ylabel('% of data') 
                    end
                end
                ttl.Units = 'Normalize'; 
                ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
                ttl.HorizontalAlignment = 'left';
            
        end
end

set(gcf,'color','w')


end