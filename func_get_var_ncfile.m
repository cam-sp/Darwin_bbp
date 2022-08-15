function var_out=func_get_var_ncfile(pathname,filename,t_vec,var,idxstart,idx_wb)

var_out=NaN(360,160,12);
if strcmp(var,'bb') || strcmp(var,'bbplk') || strcmp(var,'bbprt')
    bbwater_file=importdata('files_model_run/optics_water_Aug2014_bandave25.txt');
    bb_wat=bbwater_file.data(:,3).*0.5;
    
    for imon=1:12
     file_extension=strcat('00000',num2str(t_vec(imon)),'.nc');
     file_name=strcat(pathname,filename,file_extension);
     st=idxstart;    
     if idx_wb<10
         numvar='00';
     else
         numvar='0';
     end
     varnew=strcat(var,numvar,num2str(idx_wb));
     c = squeeze(ncread(file_name,varnew));
     c = squeeze(c(:,:,1));
     c(c==0)=NaN;
     var_out(:,:,imon)=c;
    end  
    if strcmp(var,'bb')
        var_out=var_out-bb_wat(idx_wb);
    end
end
end