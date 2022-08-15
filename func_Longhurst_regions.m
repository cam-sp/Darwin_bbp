function [idx_darwin_biomes,idx_darwin_basins,...
    idx_argo_biomes,idx_argo_basins,...
    idx_sat_biomes,idx_sat_basins...
    ]=func_Longhurst_regions()

load Longhurst_darwin_idx.mat

idx_darwin_biomes.trop=Longhurst_darwin_idx{1};
idx_darwin_biomes.oligo=Longhurst_darwin_idx{2};
idx_darwin_biomes.temp=Longhurst_darwin_idx{3};
idx_darwin_biomes.polar=Longhurst_darwin_idx{4};


% get Ocean basins
% %1- polar north
% % 2:21 - Atlantic ocean
% % 22:30 - Indian ocean
% % 31:50 - Pacific ocean
% % 50:54 - Southern ocean
% 
% s = shaperead('Longhurst/longhurst_v4_2010/Longhurst_world_v4_2010.shp');
% 
% idx_altalntic=2:21;
% idx_indian=22:30;
% idx_pacific=31:50;
% idx_southern=50:54;
% 
% lat_mat=repmat(lat,[1,360]);
% lon_mat=repmat(lon',[160,1]);
% 
% idx_basin_Atlantic=zeros(360,160);
% for i=idx_altalntic
%     in = inpolygon(wrapTo180(lon_mat),lat_mat,s(i).X,s(i).Y);
%     in=in';
%     idx_basin_Atlantic(in==1)=1;
% end
% 
% idx_basin_Indian=zeros(360,160);
% for i=idx_indian
%     in = inpolygon(wrapTo180(lon_mat),lat_mat,s(i).X,s(i).Y);
%     in=in';
%     idx_basin_Indian(in==1)=1;
% end
% 
% idx_basin_Pacific=zeros(360,160);
% for i=idx_pacific
%     in = inpolygon(wrapTo180(lon_mat),lat_mat,s(i).X,s(i).Y);
%     in=in';
%     idx_basin_Pacific(in==1)=1;
% end
% 
% idx_basin_Southern=zeros(360,160);
% for i=idx_southern
%     in = inpolygon(wrapTo180(lon_mat),lat_mat,s(i).X,s(i).Y);
%     in=in';
%     idx_basin_Southern(in==1)=1;
% end
% 
% idx_basin_Northpolar=zeros(360,160);
% for i=1
%     in = inpolygon(wrapTo180(lon_mat),lat_mat,s(i).X,s(i).Y);
%     in=in';
%     idx_basin_Northpolar(in==1)=1;
% end
% 
% idx_darwin_basins.Atlantic=idx_basin_Atlantic;
% idx_darwin_basins.Indian=idx_basin_Indian;
% idx_darwin_basins.Pacific=idx_basin_Pacific;
% idx_darwin_basins.Southern=idx_basin_Southern;
% idx_darwin_basins.Northpolar=idx_basin_Northpolar;
% 
% save('idx_darwin_basins.mat','idx_darwin_basins');
load('idx_darwin_basins.mat')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Argo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Longhurst_argo_idx.mat

idx_argo_trop=Longhurst_argo_idx{1};
idx_argo_oligo=Longhurst_argo_idx{2};
idx_argo_temp=Longhurst_argo_idx{3};
idx_argo_polar_north=Longhurst_argo_idx{4};
idx_argo_polar_south=Longhurst_argo_idx{5};
idx_argo_polar=Longhurst_argo_idx{6};

idx_argo_biomes.trop=logical(idx_argo_trop);
idx_argo_biomes.oligo=logical(idx_argo_oligo);
idx_argo_biomes.temp=logical(idx_argo_temp);
idx_argo_biomes.polar=logical(idx_argo_polar);
idx_argo_biomes.polar_N=logical(idx_argo_polar_north);
idx_argo_biomes.polar_S=logical(idx_argo_polar_south);

% get ocean basins argo
% s = shaperead('Longhurst/longhurst_v4_2010/Longhurst_world_v4_2010.shp');
% 
% idx_altalntic=2:21;
% idx_indian=22:30;
% idx_pacific=31:50;
% idx_southern=50:54;
% 
% 
% idx_argo_basin_Atlantic=zeros(size(lat_argo));
% for i=idx_altalntic
%     in = inpolygon(lon_argo,lat_argo,s(i).X,s(i).Y);
%     in=in';
%     idx_argo_basin_Atlantic(in==1)=1;
% end
% 
% idx_argo_basin_Indian=zeros(size(lat_argo));
% for i=idx_indian
%     in = inpolygon(lon_argo,lat_argo,s(i).X,s(i).Y);
%     in=in';
%     idx_argo_basin_Indian(in==1)=1;
% end
% 
% idx_argo_basin_Pacific=zeros(size(lat_argo));
% for i=idx_pacific
%     in = inpolygon(lon_argo,lat_argo,s(i).X,s(i).Y);
%     in=in';
%     idx_argo_basin_Pacific(in==1)=1;
% end
% 
% idx_argo_basin_Southern=zeros(size(lat_argo));
% for i=idx_southern
%     in = inpolygon(lon_argo,lat_argo,s(i).X,s(i).Y);
%     in=in';
%     idx_argo_basin_Southern(in==1)=1;
% end
% 
% idx_argo_basin_Northpolar=zeros(size(lat_argo));
% for i=1
%     in = inpolygon(lon_argo,lat_argo,s(i).X,s(i).Y);
%     in=in';
%     idx_argo_basin_Northpolar(in==1)=1;
% end
% 
% idx_argo_basins.Atlantic=idx_argo_basin_Atlantic;
% idx_argo_basins.Indian=idx_argo_basin_Indian;
% idx_argo_basins.Pacific=idx_argo_basin_Pacific;
% idx_argo_basins.Southern=idx_argo_basin_Southern;
% idx_argo_basins.Northpolar=idx_argo_basin_Northpolar;
% 
% save('idx_argo_basins.mat','idx_argo_basins');

load('idx_argo_basins.mat')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Sat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do interp to reduce run time
% 
% lon_q=linspace(min(lon_sat),max(lon_sat),100);
% lat_q=flip(linspace(min(lat_sat),max(lat_sat),50));
% % lon_q_2d=repmat(lon_sat,[1,100]);
% % lat_q_2d=repmat(lat_sat',[50,1]);
% % 
% % [lon_q_2d, lat_q_2d]=meshgrid(lon_q,lat_q);
% % [lon_sat_2d, lat_sat_2d]=meshgrid(lon_sat,lat_sat);
% 
% [lat_q_2d, lon_q_2d]=meshgrid(lat_q,lon_q);
% [lat_sat_2d, lon_sat_2d]=meshgrid(lat_sat,lon_sat);
% 
% idx_sat_basin_Atlantic=zeros(size(lat_q_2d));
% for i=idx_altalntic
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_basin_Atlantic(in==1)=1;
% end
% Atl_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_basin_Atlantic,lat_sat_2d,lon_sat_2d);
% Atl_ext(Atl_ext>0)=1;
% 
% idx_sat_basin_Indian=zeros(size(lat_q_2d));
% for i=idx_indian
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_basin_Indian(in==1)=1;
% end
% Ind_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_basin_Indian, lat_sat_2d,lon_sat_2d);
% Ind_ext(Ind_ext>0)=1;
% 
% idx_sat_basin_Pacific=zeros(size(lat_q_2d));
% for i=idx_pacific
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_basin_Pacific(in==1)=1;
% end
% Pac_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_basin_Pacific, lat_sat_2d,lon_sat_2d);
% Pac_ext(Pac_ext>0)=1;
% 
% idx_sat_basin_Southern=zeros(size(lat_q_2d));
% for i=idx_southern
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_basin_Southern(in==1)=1;
% end
% S_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_basin_Southern, lat_sat_2d,lon_sat_2d);
% Pac_S(S_ext>0)=1;
% 
% idx_sat_basin_Northpolar=zeros(size(lat_q_2d));
% for i=1
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_basin_Northpolar(in==1)=1;
% end
% N_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_basin_Northpolar, lat_sat_2d,lon_sat_2d);
% N_ext(N_ext>0)=1;
% 
% idx_sat_basins.Atlantic=Atl_ext;
% idx_sat_basins.Indian=Ind_ext;
% idx_sat_basins.Pacific=Pac_ext;
% idx_sat_basins.Southern=S_ext;
% idx_sat_basins.Northpolar=N_ext;
% % 
% save('idx_sat_basins.mat','idx_sat_basins');

load('idx_sat_basins.mat')


% find biomes
% 


% s = shaperead('Longhurst/longhurst_v4_2010/Longhurst_world_v4_2010.shp');
% % 
% idx_polar1=[32,31,3,2,1]; 
% idx_temp1=[30,29,20,18,17,16,13,5,4]; 
% idx_trop1=[26,25,14,9,8,27]; 
% idx_oligo1=[24,23,21,15,10,7,6];
% 
% %remove coastal regions
% idx_open=zeros(1,54);
% for i=1:54
%    alk=split(s(i).ProvDescr);
%    if ~contains(alk{1},'Coastal')
%       idx_open(i)=i;
%    end   
% end
% idx_open=idx_open(idx_open~=0);
% 
% 
% idx_trop=idx_open(idx_trop1);
% idx_oligo=(idx_open(idx_oligo1));
% idx_temp=(idx_open(idx_temp1));
% idx_polar=(idx_open(idx_polar1));
% 
% lon_q=linspace(min(lon_sat),max(lon_sat),100);
% lat_q=flip(linspace(min(lat_sat),max(lat_sat),50));
% 
% [lat_q_2d, lon_q_2d]=meshgrid(lat_q,lon_q);
% [lat_sat_2d, lon_sat_2d]=meshgrid(lat_sat,lon_sat);
% 
% idx_sat_oligo=zeros(size(lat_q_2d));
% for i=idx_oligo
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_oligo(in==1)=1;
% end
% oligo_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_oligo,lat_sat_2d,lon_sat_2d);
% oligo_ext(oligo_ext>0)=1;
% 
% idx_sat_trop=zeros(size(lat_q_2d));
% for i=idx_trop
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_trop(in==1)=1;
% end
% trop_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_trop,lat_sat_2d,lon_sat_2d);
% trop_ext(trop_ext>0)=1;
% 
% idx_sat_temp=zeros(size(lat_q_2d));
% for i=idx_temp
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_temp(in==1)=1;
% end
% temp_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_temp,lat_sat_2d,lon_sat_2d);
% temp_ext(temp_ext>0)=1;
% 
% idx_sat_polar=zeros(size(lat_q_2d));
% for i=idx_polar
%     in = inpolygon(lon_q_2d,lat_q_2d,s(i).X,s(i).Y);
%     idx_sat_polar(in==1)=1;
% end
% polar_ext = interp2(lat_q_2d,lon_q_2d,idx_sat_polar,lat_sat_2d,lon_sat_2d);
% polar_ext(polar_ext>0)=1;
% 
% idx_sat_biomes.Oligo=oligo_ext;
% idx_sat_biomes.Trop=trop_ext;
% idx_sat_biomes.Temp=temp_ext;
% idx_sat_biomes.Polar=polar_ext;
% 
% % % 
% save('idx_sat_biomes.mat','idx_sat_biomes');

load('idx_sat_biomes.mat')


end