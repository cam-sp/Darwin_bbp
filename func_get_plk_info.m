function [iplk, plk_sizes]=func_get_plk_info()

iplk.phyto=1:23;
iplk.mixo=24:31;
iplk.zoo=32:47;
iplk.bact=48:50;

iplk.ft_pico=1:4;
iplk.ft_coccos=5:9;
iplk.ft_diazos=10:14;
iplk.ft_diatoms=15:23;
iplk.ft_mixos=24:31;
iplk.ft_zoos=32:47;
iplk.ft_bact=48:50;
% get plankton sizes

sbase=-1.4130;  sincr=0.513;
logvolexp=sbase:sincr:sbase+25*sincr;
logvol=10.^logvolexp; %volume
logdm=2*((3/4*logvol/pi).^(1/3)); %diameter
%%%%
 
% 1=pro; 2=syn; 3=pico-euk; 4=cocco; 5=unicell diaz; 6=tricho; 7= diatom;
% 8= mixo; 9=zoo; 10=bact
 tgraz=8; tbact=10; tmixo=8;
 tdiat=7; tcocco=4; tzoo=9;
 grps=1:10;
 num_grps=  [1 1 2 5 4 1 9 8 16 3];
 nplank=sum(num_grps);
 ngroups=length(num_grps);
 start_size=[2 3 4 6 6 10 6 8 7 1];
 col_size=[0 0 0 0 0 3 0 0 0 0];
 ityp=0;
tsize=zeros(1,nplank);
tgrp=zeros(1,nplank);
gsize=zeros(1,nplank);
for igrp=1:10
    for is=1:num_grps(igrp)
        ityp=ityp+1;
        tsize(ityp)=start_size(igrp)-1+is;
        tgrp(ityp)=igrp;
        gsize(ityp)=tsize(ityp)+col_size(igrp);
    end 
end 
 biovol=logvol(tsize);
 biodm=logdm(tsize);
 plk_sizes.ESD=biodm;

% carbon mass mmolC
Qc =  [2.566093668483939E-012,  7.789193721295913E-012,  2.364354020783690E-011,  7.176827455596965E-011,...
   2.178474622439913E-010,  6.612603841985460E-010,  2.007208581666502E-009,  6.092738029662400E-009,  1.849407034084821E-008,...
   2.178474622439913E-010,  6.612603841985460E-010,  2.007208581666502E-009,  6.092738029662400E-009,  1.849407034084821E-008,...
   2.178474622439913E-010,  6.612603841985460E-010,  2.007208581666502E-009,  6.092738029662400E-009,  1.849407034084821E-008,...
   5.613742723010094E-008,  1.704011436062453E-007,  5.172404788573356E-007,  1.570046463929720E-006,  2.007208581666502E-009,...
   6.092738029662400E-009,  1.849407034084821E-008,  5.613742723010094E-008,  1.704011436062453E-007,  5.172404788573356E-007,...
   1.570046463929720E-006,  4.765763701139332E-006,  6.612603841985464E-010,  2.007208581666503E-009,  6.092738029662404E-009,...
   1.849407034084823E-008,  5.613742723010098E-008,  1.704011436062455E-007,  5.172404788573359E-007,  1.570046463929721E-006,...
   4.765763701139336E-006,  1.446613471441439E-005,  4.391091684330803E-005,  1.332884461596116E-004,  4.045875412494646E-004,...
   1.228096532375122E-003,  3.727799151140557E-003,  1.131546759143488E-002,  8.453810434102051E-013,  2.566093668483939E-012,...
   7.789193721295913E-012];

plk_sizes.Cmoles=Qc;
plk_sizes.Cmass=Qc.*12;
plk_sizes.Vol=biovol;

%conversion for the optics files is slightly different, so we use that one
%for conversions of cross-sections
plk_sizes.C_opts=(0.109.*biovol.^0.991).*1e-9; %mgC
 
iplk.pico=[];
iplk.nano=[];
iplk.micro=[];
iplk.meso=[];
iplk.phyto_pico=[];
iplk.phyto_nano=[];
iplk.phyto_micro=[];
iplk.mixo_pico=[];
iplk.mixo_nano=[];
iplk.mixo_micro=[];
iplk.zoo_nano=[];
iplk.zoo_micro=[];
iplk.zoo_meso=[];

for i=1:50
    if biodm(i)<=2 %pico
        iplk.pico=[iplk.pico, i];
        if ismember(i,iplk.phyto)
            iplk.phyto_pico=[iplk.phyto_pico, i];
        elseif ismember(i,iplk.mixo)
            iplk.mixo_pico=[iplk.mixo_pico, i];
        end
    elseif biodm(i)>2 && biodm(i)<20 %nano
        iplk.nano=[iplk.nano, i];
        if ismember(i,iplk.phyto)
            iplk.phyto_nano=[iplk.phyto_nano, i];
        elseif ismember(i,iplk.mixo)
            iplk.mixo_nano=[iplk.mixo_nano, i];
        elseif ismember(i,iplk.zoo)
            iplk.zoo_nano=[iplk.zoo_nano, i];     
        end
    elseif biodm(i)>20 && biodm(i)<200 %micro   
        iplk.micro=[iplk.micro, i];
        if ismember(i,iplk.phyto)
            iplk.phyto_micro=[iplk.phyto_micro, i];
        elseif ismember(i,iplk.mixo)
            iplk.mixo_micro=[iplk.mixo_micro, i];
        elseif ismember(i,iplk.zoo)
            iplk.zoo_micro=[iplk.zoo_micro, i];     
        end
    elseif biodm(i)>200 %meso    
        iplk.meso=[iplk.meso, i];
        if ismember(i,iplk.zoo)
            iplk.zoo_meso=[iplk.zoo_meso, i];     
        end
    end    
end 

end