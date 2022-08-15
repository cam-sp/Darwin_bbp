function Env=func_get_environmenatl_data(pathname,t_vec)

%Temperature
sst_mon=ncread('woa18_all.nc','temperature');
sst=mean(sst_mon,4);
sst=sst(:,:,1);
sst(sst>40)=NaN;
Env.sst=sst;

sst3D=squeeze(sst_mon(:,:,1,:));
sst3D(sst3D>40)=NaN;
Env.sst3D=sst3D;

%depth grid boundaries
file_grid='grid.nc';
z=-ncread(file_grid,'Z');
z_face=zeros(size(z));
for i=1:length(z)
    if i==1
    z_face(i)=z(i).*2;
    else
    z_face(i)=z_face(i-1)+(z(i)-z_face(i-1)).*2;
    end
end
% Z=permute(-ncread(file_grid,'Z'),[3,2,1]);
Z=permute(z_face,[3,2,1]);
lat=ncread(file_grid,'Y');
lon=ncread(file_grid,'X');
% 3Dgrid=ncdisp('grid.nc')
deltaz=[z_face(1); z_face(2:end)-z_face(1:end-1)];
deltaz=permute(deltaz,[3,2,1]);
% sec_to_days=60*60*24;
% lat_2d=repmat(lat,[1,length(lon)])';
% lon_2d=repmat(lon',[length(lat),1])';
Env.lat=lat;
Env.lon=lon;

lon_map=lon;
lon_map(1)=0;
lon_map(end)=360;
Env.lon_map=lon_map;

Env.z=z;
Env.deltaz=deltaz;

% get par
par=NaN(360,160,12);
nitr=NaN(360,160,12);
file_extension=strcat('00000',num2str(t_vec(1)),'.nc');
for imon=1:12
    
    file_nitr=strcat(pathname,'3d.',file_extension);
    varp=ncread(file_nitr,'TRAC02');
    nitr(:,:,imon)=varp(:,:,1);
    
    file_extension=strcat('00000',num2str(t_vec(imon)),'.nc');
    file_par=strcat(pathname,'par.',file_extension);
    varp=ncread(file_par,'PAR004');
    par(:,:,imon)=varp(:,:,1);

end
Env.par=par;


%get bathy
bathy1=zeros(size(varp));
bathy1(varp==0)=1;
bathy=zeros(360,160);
for i=1:360
    for j=1:160
        id=find(bathy1(i,j,:)==1);
        if ~isempty(id)
            bathy(i,j)=z(id(1));
        end
    end
end
bathy_3D=repmat(bathy,[1,1,12]);
Env.bathy=bathy;
Env.bathy_3D=bathy_3D;
end