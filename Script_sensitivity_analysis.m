%% Script for "Assessing the potential for backscattering as a proxu fro phytoplankton biomass"
%  Submitted to Global Biogeochemical Cycles

% This file performs the offline sensitivity analysis of the optical
% parameters of the model. It takes a while to run, so we have saved the
% workspace in a separate file. This file is loaded here to generate figure 10 of
% the paper. The Script for the sensitivity analysis is pasted below.

% Camila Serra-Pompei 26/01/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 10 from the paper
clearvars

%load workspace from sensitivity analysis
load('Other_extra_files/ws_uncertainties_april_2.mat')

meandet=squeeze(nanmean(nanmean(cont_detr,1),2));
meandet_def=squeeze(nanmean(nanmean(cont_detr_def,1),2));

meanpico=squeeze(nanmean(nanmean(cont_pico,1),2));
meanpico_def=squeeze(nanmean(nanmean(cont_pico_def,1),2));

cmap=brewermap(4,'Spectral');
ccol=cmap(4,:);

ccol=[42, 172, 181]./255;
ccol=[0.5 0.5 0.5];

mae_def=mae_default;

cmap1=brewermap(10,'YlGnBu');
% cmap=cmap1(4,:);
ccol=[174, 205, 214]./255;%cmap1(7,:);
cedg=[140, 174, 184]./255;
cstar=[235, 177, 42]./255;

x0=0;
y0=0;
width=16;
height=9;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(2,3)

nexttile
scatter(X(:,1),mae,5,'filled','markerfacecolor',ccol,'markeredgecolor',cedg,'markerfacealpha',0.5)
hold on
scatter(0.5,mae_def,75,'p','filled','markerfacecolor',cstar,'markeredgecolor','k')
xlabel('% change in intercept \sigma_{bb,plk,i}')
ylabel('MAE')
xticks([0 0.5 1])
xticklabels({'-50%','0','+50%'})
text(0.05,0.95,'a','Units','normalized','fontweight','bold')  

nexttile
scatter(slope_vec(ceil(X(:,2).*100)),mae,5,'filled','markerfacecolor',ccol,'markeredgecolor',cedg,'markerfacealpha',0.5)
hold on
scatter(2.387,mae_def,75,'p','filled','markerfacecolor',cstar,'markeredgecolor','k')
xlabel('slope \sigma_{bb,plk} size')
xlim([min(slope_vec), max(slope_vec)])
xticks([1.8:0.3:3])
ylabel('MAE')
text(0.05,0.95,'b','Units','normalized','fontweight','bold')  

nexttile
scatter(1./(parts_vec(ceil(X(:,3).*100)).*12.*120),mae,5,'filled','markerfacecolor',ccol,'markeredgecolor',cedg,'markerfacealpha',0.5)
hold on
scatter(1./(5e-17.*12.*120),mae_def,75,'p','filled','markerfacecolor',cstar,'markeredgecolor','k')
set(gca,'xscale','log')
xlabel('Detritus parameter "q"')
ylabel('MAE')
text(0.05,0.95,'c','Units','normalized','fontweight','bold')  

nexttile
scatter(NaN,NaN)
hold on
histogram(mae,'facecolor',ccol)
xlabel('MAE')
text(0.05,0.95,'d','Units','normalized','fontweight','bold')  
ylabel('Counts')

nexttile
scatter(meandet,mae,5,'filled','markerfacecolor',ccol,'markeredgecolor',cedg,'markerfacealpha',0.5)
hold on
scatter(meandet_def,mae_def,75,'p','filled','markerfacecolor',cstar,'markeredgecolor','k')
xlabel('b_{b,phyto}/b_{bp}')
ylabel('MAE')
text(0.05,0.95,'e','Units','normalized','fontweight','bold')  

nexttile
scatter(meanpico,mae,5,'filled','markerfacecolor',ccol,'markeredgecolor',cedg,'markerfacealpha',0.5)
hold on
scatter(meanpico_def,mae_def,75,'p','filled','markerfacecolor',cstar,'markeredgecolor','k')
xlabel('b_{b,pico}/b_{b,plk}')
ylabel('MAE')
text(0.05,0.95,'f','Units','normalized','fontweight','bold') 

set(gcf,'color','w')

% set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig10_uncertainty.eps



%% script

clearvars

t_vec=29040:240:31680;
pathname='Model_outputs/';

[iplk, plk_sizes]=func_get_plk_info(); %get plk indexes and size info
env=func_get_environmenatl_data(pathname,t_vec); %get darwin environmental vars

idx_wb=3;
bb_or_biom='biomass'; %biomass
biom=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom); %in mgC m^-3
biom.tot = biom.phyto + biom.mixo + biom.zoo + biom.bact;

% get Chl
minchl=1e-3;%0.008; %10.^(-2.1);
Chl=func_get_chl(pathname,t_vec,minchl); %in mgChl m^-3

%get bbp by detritus only (in m^-1)
filename='iops.';
var='bbprt';
bbp_detr_tot=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of detrital particles only
imon=4;
% bbp_detr_tot=bbp_detr_tot(:,:,imon);

bb_data=importdata(strcat([pathname,'p-ini-char-bbspec.dat']));
bb_wb=bb_data.data./1e4;

bt_data=importdata(strcat([pathname,'p-ini-char-btspec.dat']));
bt_wb=bt_data.data;

% idx_wb=5;
bb1=bb_wb(idx_wb,:);

%for detritus calculations
opts_init=double(readmatrix('optics_detritus.txt'));
bbratio=opts_init(1,4)./mean(opts_init(:,3));
bbin1=opts_init(:,3).*bbratio;
bbin=bbin1(idx_wb);
partsize=5e-17;
bbin=bbp_detr_tot.*partsize./(biom.detr./(12.*120)+1./120);


simnum=500;
paramnums=3;
rng('default')
X = lhsdesign(simnum,paramnums); %rows is number of simulation, columns number of input params

par_vec=linspace(0.5,1.5,100);
parts_vec=logspace(-17,-16,100);
slope_vec=linspace((1-0.25).*2.387,(1+0.25).*2.387,100);

biombb=zeros(360,160,12,50);
for iimon=1:12

        file_extension=strcat('00000',num2str(t_vec(iimon)),'.nc');
        file_3D=strcat(pathname,'3d.',file_extension);

         st=21;
    for i=1:50
         var=strcat('TRAC',num2str(st));
         c = squeeze(ncread(file_3D,var));
         c=squeeze(c(:,:,1)).*12;
         biombb(:,:,iimon,i)=c; %mgC m^-3
         st=st+1;
    end
end
 


fn = fieldnames(iplk);
tic
Rsq=zeros(1,simnum);
mae=zeros(1,simnum);
mae_tot=zeros(1,simnum);
cont_detr=zeros(360,160,simnum+1);
for i=1:simnum+1
 i
if i<=simnum
    newpval=par_vec(ceil(X(i,1).*100)); % fraction of the intercept for size bbp of plankton
    newpval(2)=slope_vec(ceil(X(i,2).*100)); % slope for size bbp of plankton
    partsize=parts_vec(ceil(X(i,3).*100)); %conversions detritus

else
    newpval(1)=1;
    newpval(2)=2.387;
    partsize=5e-17;
end

bb_sensit_input=func_recalculate_bbplk_uncertainty(iplk,plk_sizes,bb1,newpval,fn);

bi=NaN(360,160,12,50);
    for ibi=1:50
       bi(:,:,:,ibi)=biombb(:,:,:,ibi).*bb_sensit_input(ibi);
    end


[bbp_detr_tot,~,bbp_detr_refr]=func_recalculate_bb_detritus(biom,partsize,bbin,bbp_detr_tot);
% bbp.tot=bbp_detr_tot+bbp.phyto+bbp.mixo+bbp.zoo+bbp.bact;
% bbp.tot=bbp.phyto+bbp.mixo+bbp.zoo+bbp.bact;
bbp.tot=bbp_detr_tot+nansum(bi,4);
% bbp.tot3=bbp.pico+bbp.nano+bbp.micro+bbp.meso;
% bbp.tot2=bbp_detr_tot2+bbp.phyto+bbp.mixo+bbp.zoo+bbp.bact;
% get total bbp (in m^-1)
% filename='iops.';
% var='bb';
% bbp.tot=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only

if i==simnum+1
    %get total bbp (in m^-1)
    filename='iops.';
    var='bb';
    bbp.tot_original_run=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only
end
% %% some plots
%do plot fit bbp-Cphyto
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=0;
log_on=0;
md=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
if i==simnum+1
    md=func_fit_regression_plot(bbp.tot_original_run(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
end
biom_phyto=biom.phyto+biom.mixo;
biom_phyto(~idx)=NaN;

bb_tot=bbp.tot;
bb_tot(~idx)=NaN;

biom_2d=biom_phyto(:,:,imon);
mod_darwin=md.int+bb_tot(:,:,imon).*md.slope;
mod_darwin(env.bathy<500)=NaN;
mod_darwin(mod_darwin<0)=NaN;    
[~, ~,~, logMAE2]=func_calc_bias(mod_darwin(:),biom_2d(:));

if i<=simnum
    mae(i)=logMAE2;
    % cont_detr(:,:,i)=bbp_detr_tot(:,:,imon)./bbp.tot(:,:,imon);
    cont_detr(:,:,i)=nansum(bi(:,:,imon,1:31),4)./bbp.tot(:,:,imon);
    bbpico=nansum(bi(:,:,imon,iplk.phyto_pico),4);
    cont_pico(:,:,i)=bbpico./nansum(bi(:,:,imon,:),4);%bbp.tot(:,:,imon);
    cont_piconano(:,:,i)=(bbpico+ nansum(bi(:,:,imon,iplk.phyto_nano),4))./bbp.tot(:,:,imon);
    Rsq(i)=func_Rsq(log10(biom_2d(:)),log10(mod_darwin(:)));
    mae_tot(i)=md.MAElog;
else
    mae_default=logMAE2;
    cont_detr_def=nansum(bi(:,:,imon,1:31),4)./bbp.tot(:,:,imon);
    bbpico=nansum(bi(:,:,imon,iplk.phyto_pico),4);
    cont_pico_def=bbpico./nansum(bi(:,:,imon,:),4);%bbp.tot(:,:,imon);
    cont_piconano_def=(bbpico+ nansum(bi(:,:,imon,iplk.phyto_nano),4))./bbp.tot(:,:,imon);
end


end
toc