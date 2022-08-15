function [argo, ArgoDarwin, ArgoSat]=func_get_Argo(bbp,biom,Chl,env,sat,idx_argo_biomes)

% get Argo data
addpath('C:\Users\Camila\Desktop\Backup\Projects\Argo\BGC_Argo_Mat_Toolbox-main\BGC_Argo_Mat_Toolbox-main/')
load Argo_matrix.mat

Chl_argo=Argo_matrix(:,1);
bbp_argo=Argo_matrix(:,2);
lat_argo=Argo_matrix(:,3);
lon_argo=Argo_matrix(:,4);
day_argo=Argo_matrix(:,5);
day_argo=datetime(day_argo,'convertfrom','modifiedjuliandate')+years(92);
month_argo = month(day_argo);
lon_argo_360=wrapTo360(lon_argo);


% idx_argo_used=logical(idx_argo_trop+ idx_argo_oligo+ idx_argo_temp+ idx_argo_polar_north+ idx_argo_polar_south);
idx_argo_used=idx_argo_biomes.trop+idx_argo_biomes.oligo+ idx_argo_biomes.temp+ idx_argo_biomes.polar;
%remove outliers, coastal regions and the mediterranean sea
Chl_argo(idx_argo_used==0)=NaN;
bbp_argo(idx_argo_used==0)=NaN;
Chl_argo(Chl_argo<1e-3)=NaN;


%match Darwin output with Argo data
bbp_argo_binned=NaN(size(lon_argo));
bbp_argo_binned450=NaN(size(lon_argo));
chl_argo_binned=NaN(size(lon_argo));
bbplk_argo_binned=NaN(size(lon_argo));
poc_argo_binned=NaN(size(lon_argo));
bathy_argo_binned=NaN(size(lon_argo));
biom_phyto_binned=NaN(size(lon_argo));
biom_plk_binned=NaN(size(lon_argo));
sst_binned=NaN(size(lon_argo));
for i=1:length(lon_argo)
    
    if idx_argo_used(i)==1
     [ ~, idx_lon] = min( abs( env.lon-lon_argo_360(i) ) );
     [ ~, idx_lat] = min( abs( env.lat-lat_argo(i) ) );
     var=bbp.tot700(idx_lon,idx_lat,month_argo(i));
     bbp_argo_binned(i)=var;
     
     var=bbp.tot450(idx_lon,idx_lat,month_argo(i));
     bbp_argo_binned450(i)=var;
     
     var=Chl.tot(idx_lon,idx_lat,month_argo(i));
     chl_argo_binned(i)=var;
     
%      var=bbplk_darwin700(idx_lon,idx_lat,month_argo(i));
%      bbplk_argo_binned(i)=var;

%      var=poc_mon(idx_lon,idx_lat,month_argo(i));
%      poc_argo_binned(i)=var;  
     
     var=env.bathy(idx_lon,idx_lat);
     bathy_argo_binned(i)=var;  
     
     var=squeeze(biom.phyto(idx_lon,idx_lat,month_argo(i))+biom.mixo(idx_lon,idx_lat,month_argo(i)));
     biom_phyto_binned(i)=var;
     
%      var=squeeze(biom_plk_t(idx_lon,idx_lat,month_argo(i)));
%      biom_plk_binned(i)=var.*12;
     
     sst_binned(i)=env.sst3D(idx_lon,idx_lat,month_argo(i));
    end
end


argo.Chl=Chl_argo;
argo.bbp=bbp_argo;
argo.lat=lat_argo;
argo.lon=lon_argo;
argo.lon360=lon_argo_360;
argo.day=day_argo;
argo.month=month_argo;
argo.idx_used=idx_argo_used;

ArgoDarwin.bbp700=bbp_argo_binned;
ArgoDarwin.bbp450=bbp_argo_binned450;
ArgoDarwin.Chl=chl_argo_binned;
ArgoDarwin.biom_phyto=biom_phyto_binned;
ArgoDarwin.bathy=bathy_argo_binned;
ArgoDarwin.sst=sst_binned;

% bbp_argo_binned(bathy_argo_binned<500)=NaN;
% bbp_argo_binned450(bathy_argo_binned<500)=NaN;
% chl_argo_binned(bathy_argo_binned<500)=NaN;



%match Sat output with Argo data
bbp_argo_binned450=NaN(size(lon_argo));
chl_argo_binned=NaN(size(lon_argo));
for i=1:length(lon_argo)
    
    if idx_argo_used(i)==1
     [ ~, idx_lon] = min( abs( sat.lon-lon_argo_360(i) ) );
     [ ~, idx_lat] = min( abs( sat.lat-lat_argo(i) ) );
     
     var=sat.bbp(idx_lon,idx_lat,month_argo(i));
     bbp_argo_binned450(i)=var;
     
     var=sat.Chl(idx_lon,idx_lat,month_argo(i));
     chl_argo_binned(i)=var;
 
    end
end

ArgoSat.bbp=bbp_argo_binned450;
ArgoSat.Chl=chl_argo_binned;



end