function func_figure_percentdiff_maps_Chl(biom_phyto,bb_tot,idx_out,md,md_chl,bathy,lat,lon)

biom_phyto(~idx_out)=NaN;
bb_tot(~idx_out)=NaN;

addpath C:\Users\Camila\Desktop\Backup\Projects\3Doutputs\cmap_ocean
monstr={'January','February','March','April','May','June','July','August','September','October','November','December'};

% tlvec=[1:4:16,2:4:16,3:4:16,4:4:16];

tlvec=[1:3:12,2:3:12,3:3:12];

% tlvec=[1:2:8,2:2:8];
st=1;

RsqLL=zeros(12,1);
RsqAL=zeros(12,1);
maeLL=zeros(12,1);
maeAL=zeros(12,1);

minlat2=-50;
maxlat2=65;
st=1;
figure
tiledlayout(4,3,'Padding','Compact');%,'TileSpacing','Compact','Padding','Compact'

% strvec_letters={'a. ','c. ','e. ','g. ','b. ','d. ','f. ','h. '};
strvec_letters={'a. ','b. ','c. ','d. ','e. ','f. ','g. ','h. ','i. ','j. ','k. ','l. ','m. ','n. ','o. '};

fsiz=9;
for imon=1:3:12
    biom_2d=biom_phyto(:,:,imon);    
    
   p=nexttile(tlvec(st));
   function_plot_maps_sat(log10(biom_2d),lat, lon)
   hold on
    caxis([0 2.5])
    ylabel('Alk')
    ttl =title(strcat(...
        strvec_letters{st},...
        monstr{imon}...
        )...
        ,'fontweight','normal'); 
       
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
    colormap(p,viridis(100))
    cmap=flip(cmocean('deep',100));
    colormap(p,cmap(2:end,:))
       st=st+1;
    

end
    h=colorbar('horizontal'); 
    h.Ticks= [0, 1, 2];
    h.TickLabels={'1','10','100'};
    ylabel(h,'C_{phyto} [mgC m^{-3}]')

for imon=1:3:12   
    biom_2d=biom_phyto(:,:,imon);
    bbpart_2d=bb_tot(:,:,imon);
    mod_darwin=md.int+bbpart_2d.*md.slope;
    mod_darwin(bathy<500)=NaN;
    mod_darwin(mod_darwin<0)=NaN;    
    [~, ~,~, logMAE2]=func_calc_bias(mod_darwin(:),biom_2d(:));
    Rsq2=func_Rsq(log10(biom_2d(:)),log10(mod_darwin(:)));
    
    RsqAL(st)=Rsq2;
    maeAL(st)=logMAE2;
    
   p=nexttile(tlvec(st));
   function_plot_maps_sat(100.*(mod_darwin-biom_2d)./biom_2d,lat, lon)
   hold on
    caxis([-60 60])
ttl =title(strcat([...
        strvec_letters{st},...
        ' R^2=',num2str(Rsq2,'%.2f'),...
        ', MAE=',num2str(logMAE2,'%.2f')]),'fontweight','normal','fontsize',7);    

    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
    colormap(p,flip(brewermap(100,'RdBu')))
   st=st+1;
end
    h=colorbar('horizontal'); %[left, bottom, width, height].
     ylabel(h,'% Difference')
set(findall(gcf,'-property','FontSize'),'FontSize',9)




for imon=1:3:12   
    biom_2d=biom_phyto(:,:,imon);
    bbpart_2d=bb_tot(:,:,imon);
    mod_darwin=md_chl.int+bbpart_2d.*md_chl.slope;
    mod_darwin(bathy<500)=NaN;
    mod_darwin(mod_darwin<0)=NaN;    
    [~, ~,~, logMAE2]=func_calc_bias(mod_darwin(:),biom_2d(:));
    Rsq2=func_Rsq(log10(biom_2d(:)),log10(mod_darwin(:)));
    
    RsqAL(st)=Rsq2;
    maeAL(st)=logMAE2;
    
   p=nexttile(tlvec(st));
   function_plot_maps_sat(100.*(mod_darwin-biom_2d)./biom_2d,lat, lon)
   hold on
    caxis([-60 60])
ttl =title(strcat([...
        strvec_letters{st},...
        ' R^2=',num2str(Rsq2,'%.2f'),...
        ', MAE=',num2str(logMAE2,'%.2f')]),'fontweight','normal','fontsize',7);    

    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
    colormap(p,flip(brewermap(100,'RdBu')))
   st=st+1;
end
    h=colorbar('horizontal'); %[left, bottom, width, height].
     ylabel(h,'% Difference')
% set(findall(gcf,'-property','FontSize'),'FontSize',9)

end