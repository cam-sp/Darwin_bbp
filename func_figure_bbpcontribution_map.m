function func_figure_bbpcontribution_map(bb_tot,bb_phyto,bb_bact,bb_detr,bb_detr_refr,lat,lon)
% contributions by other things

addpath 'C:\Users\Camila\Desktop\Backup\Projects\3Doutputs\cmap_ocean'
monstr={'January','February','March','April','May','June','July','August','September','October','November','December'};

lettervec={'a b c d e f g h i j k l m n o p q r s t u'};
lettervec=split(lettervec,' ');


x0=0;
y0=0;
width=18;
height=8;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(2,3,'TileSpacing','Compact');%,'TileSpacing','Compact','Padding','Compact'
st=0;
for imon=[1, 7]
    bbpart_2d=bb_tot(:,:,imon);%-(bb_phyto(:,:,imon)+bb_bact(:,:,imon)+bb_zoo(:,:,imon));
    for j=1:3
        if j==1
            bbphyto_2d=bb_phyto(:,:,imon);
        elseif j==2
            bbphyto_2d=bb_bact(:,:,imon);
        elseif j==3
            bbphyto_2d=bb_detr(:,:,imon);
        elseif j==4
            bbphyto_2d=bb_detr_refr(:,:,imon);
        end

        nexttile
        function_plot_maps_sat((bbphyto_2d./bbpart_2d).*100,lat, lon)
%         caxis([0 0.4])
        st=st+1;
        if j==1
            if imon==1
            ttl =title({strcat([lettervec{st},'. ','Phytoplankton'])});
            elseif imon>1
                ttl =title({strcat([lettervec{st}])});
%                 h=colorbar('horizontal');
%                 h.Ticks=0:10:50;
            end
            caxis([0 100])
        elseif j==2
            if imon==1
            ttl =title({strcat([lettervec{st},'. ','H.bacteria'])});
            elseif imon>1
                ttl =title({strcat([lettervec{st}])});
%                 h=colorbar('horizontal');
%                 ylabel(h,'% contribution to b_{bp}')
%                 h.Ticks=0:10:50;
            end
            caxis([0 100])
        elseif j==3
            if imon==1
            ttl =title({strcat([lettervec{st},'. ','Detritus'])});
            elseif imon>1
                ttl =title({strcat([lettervec{st}])});
%                 h=colorbar('horizontal');
%                 h.Ticks=0:25:100;
            end
            caxis([0 100])
        elseif j==4
            if imon==1
            ttl =title({strcat([lettervec{st},'. ','Refractory detritus'])});
            elseif imon>1
                ttl =title({strcat([lettervec{st}])});
%                 h=colorbar('horizontal');
%                 h.Ticks=0:25:100;
            end
            caxis([0 100])
        end
        ttl.Units = 'Normalize'; 
        ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
        ttl.HorizontalAlignment = 'left';  
        ttl.FontWeight='normal';
%         ttl.FontSize=8;
    end
end
    h=colorbar('horizontal','Position',[0.3 0.1 0.4 0.03]); %[left, bottom, width, height].
     ylabel(h,'% contribution to b_{bp}')
%     h=colorbar('vertical','Position',[0.95 0.3 0.03 0.4]); %[left, bottom, width, height].
%     h.YTick = [-1 -0.5 0 0.5 1];
%     h.YTickLabel = {'1/10', '1/3', '1', '3', '10'};
%     colormap(viridis(100))
%     colormap(flip(brewermap(100,'GnBu')))
%     colormap(flip(brewermap(100,'YlGnBu')))
        cmap=flip(cmocean('deep',100));
    colormap(cmap(2:end,:))
%     set(findall(gcf,'-property','FontSize'),'FontSize',9)

% saveFigPdf(fig, 'testfig')
end