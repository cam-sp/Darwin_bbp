function bbp=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom,sensit,bb_sensit_input,imon_init) 


bb_data=importdata(strcat([pathname,'p-ini-char-bbspec.dat']));
% bb_data=importdata('p-ini-char-bbspec.dat');
bb_wb=bb_data.data./1e4;

bb_phyto=zeros(360,160,12);
bb_pico_phyto=zeros(360,160,12);
bb_nano_phyto=zeros(360,160,12);
bb_micro_phyto=zeros(360,160,12);

bb_mixo=zeros(360,160,12);
bb_nano_mixo=zeros(360,160,12);
bb_micro_mixo=zeros(360,160,12);

bb_zoo=zeros(360,160,12);
bb_nano_zoo=zeros(360,160,12);
bb_micro_zoo=zeros(360,160,12);
bb_meso_zoo=zeros(360,160,12);

bb_bact=zeros(360,160,12);

bb_pico=zeros(360,160,12);
bb_nano=zeros(360,160,12);
bb_micro=zeros(360,160,12);
bb_meso=zeros(360,160,12);

bb_diatoms=zeros(360,160,12);
bb_coccos=zeros(360,160,12);

bb_pl=zeros(360,160,12);

bbwater_file=importdata('optics_water_Aug2014_bandave25.txt');
% bb_wat=bbwater_file.data(:,3).*0.5;
% 
% num_grps_2=  [4 5 5 9 8 16 3];
% ngroups=7;
% biodm=plk_sizes.ESD;
bb1=bb_wb(idx_wb,:);
imon_vec=1:12;
if nargin>6
    if nargin==7
    bb1(1:31)=10.^(-4.5);%1e-5;
% bb1(iplk.ft_coccos)=bb1(iplk.ft_coccos).*10;
% bb1(iplk.ft_diatoms)=bb1(iplk.ft_diatoms).*10;
% bb1(iplk.ft_diazos)=bb1(iplk.ft_diazos).*3;

% load('bb_wb_new.mat');
% bb1=bb_wb_new(idx_wb,:);

% % bb_data=importdata('C:\Users\Camila\Desktop\Backup\Projects\3Doutputs\p_ini_files/p-ini-char-bbspec.dat');
% bb_data=importdata('C:\Users\Camila\Desktop\NO backup\Darwin_runs/MONTHS_iron_rpoc1_run6/p-ini-char-bbspec.dat');
% bb_wb=bb_data.data./1e4;
% bb1=bb_wb(idx_wb,:);
    elseif nargin==8
        bb1=bb_sensit_input;
        imon_vec=1:12;%imon_init;
    end
end

    for imon=imon_vec

    file_extension=strcat('00000',num2str(t_vec(imon)),'.nc');
    file_3D=strcat(pathname,'3d.',file_extension);

        st=21;
        st2=1;
        
        for i=1:50
             var=strcat('TRAC',num2str(st));
             c = squeeze(ncread(file_3D,var));
             if bb_or_biom==1
                c=squeeze(c(:,:,1)).*12.*bb1(st2);
                bb_pl(:,:,imon)=bb_pl(:,:,imon)+c;
             elseif bb_or_biom==2
                c=squeeze(c(:,:,1)).*12;
             end
%              c(c==0)=NaN; 
                          
            if ismember(i,iplk.phyto)
                bb_phyto(:,:,imon)=bb_phyto(:,:,imon)+c;
                if ismember(i,iplk.pico)
                    bb_pico_phyto(:,:,imon)=bb_pico_phyto(:,:,imon)+c;
                elseif ismember(i,iplk.nano)
                    bb_nano_phyto(:,:,imon)=bb_nano_phyto(:,:,imon)+c;
                elseif ismember(i,iplk.micro)
                    bb_micro_phyto(:,:,imon)=bb_micro_phyto(:,:,imon)+c;
                end
            elseif ismember(i,iplk.mixo)
                bb_mixo(:,:,imon)=bb_mixo(:,:,imon)+c;
                if ismember(i,iplk.nano)
                    bb_nano_mixo(:,:,imon)=bb_nano_mixo(:,:,imon)+c;
                elseif ismember(i,iplk.mixo)
                    bb_micro_mixo(:,:,imon)=bb_micro_mixo(:,:,imon)+c;
                end
            elseif ismember(i,iplk.zoo) %mixos and zoos
                    bb_zoo(:,:,imon)=bb_zoo(:,:,imon)+c;
                    if ismember(i,iplk.nano)
                        bb_nano_zoo(:,:,imon)=bb_nano_zoo(:,:,imon)+c;
                    elseif ismember(i,iplk.micro)
                        bb_micro_zoo(:,:,imon)=bb_micro_zoo(:,:,imon)+c;
                    elseif ismember(i,iplk.meso)
                        bb_meso_zoo(:,:,imon)=bb_meso_zoo(:,:,imon)+c;
                    end
             elseif ismember(i,iplk.bact) %bacteria
                    bb_bact(:,:,imon)=bb_bact(:,:,imon)+c;
             end

            if ismember(i,iplk.pico)
                bb_pico(:,:,imon)=bb_pico(:,:,imon)+c;
            elseif ismember(i,iplk.nano)
                bb_nano(:,:,imon)=bb_nano(:,:,imon)+c;
            elseif ismember(i,iplk.micro)
                bb_micro(:,:,imon)=bb_micro(:,:,imon)+c;
            elseif ismember(i,iplk.meso)
                bb_meso(:,:,imon)=bb_meso(:,:,imon)+c;
            end
            
            if ismember(i,iplk.ft_diatoms)
                bb_diatoms(:,:,imon)=bb_diatoms(:,:,imon)+c;
            elseif ismember(i,iplk.ft_coccos)
                bb_coccos(:,:,imon)=bb_coccos(:,:,imon)+c;
            end

                st=st+1;
                st2=st2+1;
                    
        end
        
        if bb_or_biom==2
            var=strcat('TRAC12');
            c = double(squeeze(ncread(file_3D,var)));
            c=squeeze(c(:,:,1)).*12;
            c(c==0)=NaN;
            bbp.detr(:,:,imon)=c;
        elseif bb_or_biom==1
%             ppart=3.98e-14;%3.98e-17;
%             bb_strams=2.06000e-18;
%             var=strcat('TRAC12');
%             c = squeeze(ncread(file_3D,var));
%             c=squeeze(c(:,:,1)).*bb_strams./(ppart.*120);
%             c(c==0)=NaN;
%             bbp.detr(:,:,imon)=c;
        end

    end
    
 bbp.pl=bb_pl;   
        
bbp.phyto=bb_phyto;
bbp.phyto_pico=bb_pico_phyto;
bbp.phyto_nano=bb_nano_phyto;
bbp.phyto_micro=bb_micro_phyto;

bbp.mixo=bb_mixo;
bbp.mixo_nano=bb_nano_mixo;
bbp.mixo_micro=bb_micro_mixo;

bbp.zoo=bb_zoo;
bbp.zoo_nano=bb_nano_zoo;
bbp.zoo_micro=bb_micro_zoo;
bbp.zoo_meso=bb_meso_zoo;

bbp.bact=bb_bact;

bbp.pico=bb_pico;
bbp.nano=bb_nano;
bbp.micro=bb_micro;
bbp.meso=bb_meso;

bbp.coccos=bb_coccos;
bbp.diatoms=bb_diatoms;

%to avoid having 0 in land areas
fn = fieldnames(bbp);
for i=1:length(fn)
   alk=bbp.(fn{i});
   alk(alk==0)=NaN;
   bbp.(fn{i})=alk;
end

end