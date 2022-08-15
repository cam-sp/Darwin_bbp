function Chl=func_get_chl(pathname,t_vec,minchl)
% get surface chl

chl_phyto_t=zeros(360,160,12);
chl_mixos_t=zeros(360,160,12);    
for imon=1:12
         file_extension=strcat('00000',num2str(t_vec(imon)),'.nc');
         file_chl=strcat(pathname,'chl.',file_extension);

    %Chl of phyto only
    for i=71:94

        var=strcat('TRAC',num2str(i));
        c = squeeze(ncread(file_chl,var));
        chl_phyto_t(:,:,imon)=chl_phyto_t(:,:,imon)+c(:,:,1);

    end

    %Chl of mixos only
    for i=95:99%94

        var=strcat('TRAC',num2str(i));
        c = squeeze(ncread(file_chl,var));
        chl_mixos_t(:,:,imon)=chl_mixos_t(:,:,imon)+c(:,:,1);

    end

        mat=ncread(file_chl,'TRAC0a');
        chl_mixos_t(:,:,imon)=chl_mixos_t(:,:,imon)+mat(:,:,1);
        mat=ncread(file_chl,'TRAC0b');
        chl_mixos_t(:,:,imon)=chl_mixos_t(:,:,imon)+mat(:,:,1);

end
        
        chl_phyto_t(chl_phyto_t==0)=NaN;
        chl_mixos_t(chl_mixos_t==0)=NaN;

        Chl.phyto=chl_phyto_t;
        Chl.mixos=chl_mixos_t;
        Chl.tot= chl_mixos_t + chl_phyto_t;
        
        Chl.idx_minChl=Chl.tot>minchl;

end