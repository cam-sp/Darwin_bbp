clearvars

t_vec=29040:240:31680;
pathname='C:\Users\Camila\Desktop\NO backup\Darwin_runs/MONTHS_iron_rpoc1_run9/';

[iplk, plk_sizes]=func_get_plk_info(); %get plk indexes and size info
env=func_get_environmenatl_data(pathname,t_vec); %get darwin environmental vars

wb_ranges=400-12.5:25:700+12.5; %waveband limits
wb=400:25:700; %mean waveength in each waveband interval
idx_wb=7; %index of waveband
bb_or_biom=1;%bb
bbp=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom); %in m^-1

%get total bbp (in m^-1)
filename='iops.';
var='bb';
bbp.tot=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only

%get bbp by detritus only (in m^-1)
filename='iops.';
var='bbprt';
bbp.detr_tot=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of detrital particles only

%get bbp by plankton only (in m^-1)
filename='iops.';
var='bbplk';
bbp.plk=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of plk only

bb_or_biom=2; %biomass
biom=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom); %in mgC m^-3
biom.tot = biom.phyto + biom.mixo + biom.zoo + biom.bact;

% get Chl
minchl=1e-3;%0.008; %10.^(-2.1);
Chl=func_get_chl(pathname,t_vec,minchl); %in mgChl m^-3

% bbp.tot2=bbp.detr_tot+bbp.phyto+bbp.mixo+bbp.zoo+bbp.bact;
% 
% bbp.plk2=bbp.phyto+bbp.mixo+bbp.zoo+bbp.bact;

% %% some plots
%do plot fit bbp-Cphyto
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=1;
log_on=0;
figure
md=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
ylabel('C_{phyto} [mgC m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
xlim([-4 -1])
xticks([-4:1:-1]);
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
yticks([-2:1:4]);
yticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});
xlim([-3.5 -1.5]) 
% 
% %do percent difference plots
idx=Chl.idx_minChl & env.bathy_3D>500;
func_figure_percentdiff_maps(biom.phyto+biom.mixo,bbp.tot,idx,md,env.bathy,env.lat,env.lon_map)

% %do contribution by each plankton to bbp plot
func_figure_bbpcontribution_map(bbp.tot,bbp.phyto+bbp.mixo,bbp.bact,bbp.detr_tot,bbp.zoo,env.lat,env.lon_map)
text(0.01,1,'January','Units','normalized','Rotation',90)  
text(0.7,1,'July','Units','normalized','Rotation',90) 

 h=text(1,1,'  My description');
 set(h,'Rotation',90);

% do plot decompose Cphyto-bbp relationship and effects of phyto sizes on bbp
idx=Chl.idx_minChl & env.bathy_3D>500;
func_figure_decompose_Cphytobbp(pathname,t_vec,bbp,biom,iplk,plk_sizes,Chl,env,idx)

%% get Argo data

idx_wb=13;
filename='iops.';
var='bb';
bbp.tot700=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only

idx_wb=3;
filename='iops.';
var='bb';
bbp.tot450=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only

[idx_darwin_biomes,idx_darwin_basins,...
    idx_argo_biomes,idx_argo_basins,...
    idx_sat_biomes,idx_sat_basins...
    ]=func_Longhurst_regions();

% get sat data
sat=func_get_sat();

%get argo data and matchings with Darwin and with sat
[argo, ArgoDarwin, ArgoSat]=func_get_Argo(bbp,biom,Chl,env,sat,idx_argo_biomes);

%plot argo data and cumsum plots
func_figure_argo_data(argo,idx_argo_biomes)



%% get mean by regions

means=func_get_mean_regions(bbp,biom,argo,ArgoDarwin,sat,ArgoSat,...
    idx_argo_biomes,idx_argo_basins,idx_sat_biomes,idx_sat_basins);

func_figure_data_comparison_2(means)


%%
figure
tiledlayout(1,2)
nexttile
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=1;
log_on=0;
md_cphyto=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
ylabel('C_{phyto} [mgC m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
ylim([-1.5 3.5])
xlim([-4 -1])
xticks([-4:1:-1]);
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
yticks([-2:1:4]);
yticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});
xlim([-3.5 -1.5]) 
text(0.9,0.05,'a','Units','normalized','fontweight','bold')   

nexttile
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=1;
log_on=0;
% md=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
md=func_fit_regression_plot(bbp.tot(idx),Chl.tot(idx).*150,rob_opts,plot_on,log_on,1);
ylabel('Chl\timesQ_{C:Chl} [mg m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
% ylim([-3.5 1.5])
% xlim([-4 -1])
% xticks([-4:1:-1]);
% xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
% yticks([-4:1:2]);
% yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'});
% xlim([-3.5 -1.5]) 
ylim([-1.5 3.5])
xlim([-4 -1])
xticks([-4:1:-1]);
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
yticks([-2:1:4]);
yticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});
xlim([-3.5 -1.5]) 
text(0.9,0.05,'b','Units','normalized','fontweight','bold')  


%%

rob_opts=1;
plot_on=0;
log_on=0;
% md=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
md_chl=func_fit_regression_plot(bbp.tot(idx),Chl.tot(idx).*150,rob_opts,plot_on,log_on,1);

% %do percent difference plots
idx=Chl.idx_minChl & env.bathy_3D>500;
func_figure_percentdiff_maps_Chl(biom.phyto+biom.mixo,bbp.tot,idx,md_cphyto,md_chl,env.bathy,env.lat,env.lon_map)

%% chl:C ratios

rob_opts=1;
plot_on=1;
log_on=0;
figure
% md=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
md_chl=func_fit_regression_plot(bbp.tot(idx),Chl.tot(idx)./biom.tot(idx),rob_opts,plot_on,log_on,1);


%%
figure
tiledlayout(1,2)
nexttile
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=1;
log_on=0;
md_cphyto=func_fit_regression_plot(Chl.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
xlabel('C_{phyto} [mgC m^{-3}]')
ylabel('b_{bp}(700) [m^{-1}]')
% ylim([-1.5 3.5])
% xlim([-4 -1])
% xticks([-4:1:-1]);
% xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
% yticks([-2:1:4]);
% yticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});
% xlim([-3.5 -1.5]) 
text(0.9,0.05,'a','Units','normalized','fontweight','bold')   

nexttile
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=0;
plot_on=1;
log_on=0;
% md=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);
md=func_fit_regression_plot(Chl.tot(idx),bbp.tot(idx),rob_opts,plot_on,log_on,1);
xlabel('Chl [mg m^{-3}]')
ylabel('b_{bp}(700) [m^{-1}]')
% ylim([-3.5 1.5])
% xlim([-4 -1])
% xticks([-4:1:-1]);
% xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
% yticks([-4:1:2]);
% yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'});
% xlim([-3.5 -1.5]) 
text(0.9,0.05,'b','Units','normalized','fontweight','bold')  

%% seasonal

gc=[16, 181, 156]./255;
fsiz=9;

figure
tiledlayout(6,2,'TileSpacing','tight','Padding','compact')

nexttile
xid=325;
yid=136;%131;
alk=25097.*squeeze(bbp.tot(xid,yid,:))-7.6;
alk2=31046.*squeeze(bbp.tot(xid,yid,:))-12.4;
bi=squeeze(biom.phyto(xid,yid,:)+biom.mixo(xid,yid,:));
plot(bi,'k','linewidth',1)
hold on
plot(alk,'color',[0.5 0.5 0.5],'linewidth',1)
plot(alk2,'color',gc,'linewidth',1)
xlim([1 12])
xticks([])
set(gca, 'box', 'off')
% text(0.05,0.9,'\bf a\rm','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('a. North Atlantic');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
% [hh, icons]=legend({'C_{phyto}','C_{phyto} b_{bp}-derived'},'Location','northoutside');
% p1 = icons(1).Position;
% p2 = icons(2).Position;
% icons(1).Position = [0.3 p1(2) 0];
% icons(2).Position = [0.3 p2(2) 0];
% icons(3).XData = [0.05 0.2];
% icons(5).XData = [0.05 0.2];
% legend boxoff

nexttile
plot(100.*(alk-bi)./bi,'color',[0.5 0.5 0.5],'linewidth',1)
hold on
plot(100.*(alk2-bi)./bi,'color',gc,'linewidth',1)
plot(bi.*0,'k--')
xlim([1 12])
xticks([])
set(gca, 'box', 'off')
% text(0.05,0.9,'b','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('b');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;

nexttile
xid=195;
yid=128;%123;
alk=25097.*squeeze(bbp.tot(xid,yid,:))-7.6;
alk2=31046.*squeeze(bbp.tot(xid,yid,:))-12.4;
bi=squeeze(biom.phyto(xid,yid,:)+biom.mixo(xid,yid,:));
plot(bi,'k','linewidth',1)
hold on
plot(alk,'color',[0.5 0.5 0.5],'linewidth',1)
plot(alk2,'color',gc,'linewidth',1)
xlim([1 12])
xticks([])
set(gca, 'box', 'off')
% text(0.05,0.9,'c','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('c. North Pacific');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;

nexttile
plot(100.*(alk-bi)./bi,'color',[0.5 0.5 0.5],'linewidth',1)
hold on
plot(100.*(alk2-bi)./bi,'color',gc,'linewidth',1)
plot(bi.*0,'k--')
xlim([1 12])
xticks([])
set(gca, 'box', 'off')
% text(0.05,0.9,'d','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('d');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;

nexttile
xid=234;
yid=82;
alk=25097.*squeeze(bbp.tot(xid,yid,:))-7.6;
alk2=31046.*squeeze(bbp.tot(xid,yid,:))-12.4;
bi=squeeze(biom.phyto(xid,yid,:)+biom.mixo(xid,yid,:));
plot(bi,'k','linewidth',1)
hold on
plot(alk,'color',[0.5 0.5 0.5],'linewidth',1)
plot(alk2,'color',gc,'linewidth',1)
xlim([1 12])
xticks([])
ylim([0 20])
ylabel('Cphyto [mgC m^{-3}]')
set(gca, 'box', 'off')
% text(0.05,0.9,'e','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('e. Tropical Pacific');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;

nexttile
plot(100.*(alk-bi)./bi,'color',[0.5 0.5 0.5],'linewidth',1)
hold on
plot(100.*(alk2-bi)./bi,'color',gc,'linewidth',1)
plot(bi.*0,'k--')
xlim([1 12])
xticks([])
ylabel('% difference [-]')
set(gca, 'box', 'off')
% text(0.05,0.9,'f','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('f');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;

nexttile
xid=178;
yid=97;
alk=25097.*squeeze(bbp.tot(xid,yid,:))-7.6;
alk2=31046.*squeeze(bbp.tot(xid,yid,:))-12.4;
bi=squeeze(biom.phyto(xid,yid,:)+biom.mixo(xid,yid,:));
plot(bi,'k','linewidth',1)
hold on
plot(alk,'color',[0.5 0.5 0.5],'linewidth',1)
plot(alk2,'color',gc,'linewidth',1)
% xlabel('Month')
% ylabel('Cphyto [mgC m^{-3}]')
xlim([1 12])
xticks([])
ylim([0 15])
% legend({'C_{phyto}','C_{phyto} b_{bp}-derived'},'Location','southoutside','NumColumns',2);
set(gca, 'box', 'off')
% text(0.05,0.9,'g','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('g. N.P. oligotrophic gyre');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
% p1 = icons(1).Position;
% p2 = icons(2).Position;
% icons(1).Position = [0.3 p1(2) 0];
% icons(2).Position = [0.3 p2(2) 0];
% icons(3).XData = [0.05 0.2];
% icons(5).XData = [0.05 0.2];
% legend boxoff

nexttile
plot(100.*(alk-bi)./bi,'color',[0.5 0.5 0.5],'linewidth',1)
hold on
plot(100.*(alk2-bi)./bi,'color',gc,'linewidth',1)
plot(bi.*0,'k--')
% xlabel('Month')
% ylabel('% difference [-]')
xlim([1 12])
ylim([-50 5])
set(gca, 'box', 'off')
% text(0.05,0.9,'h','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('h');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
xticks([])

nexttile
xid=213;%254;
yid=16;%18;
alk=25097.*squeeze(bbp.tot(xid,yid,:))-7.6;
alk2=31046.*squeeze(bbp.tot(xid,yid,:))-12.4;
bi=squeeze(biom.phyto(xid,yid,:)+biom.mixo(xid,yid,:));
plot(bi,'k','linewidth',1)
hold on
plot(alk,'color',[0.5 0.5 0.5],'linewidth',1)
plot(alk2,'color',gc,'linewidth',1)
xlabel('Month')
% ylabel('Cphyto [mgC m^{-3}]')
xlim([1 12])
leg=legend({'C_{phyto}','C_{phyto} b_{bp}-derived','C_{phyto} b_{bp}-Chl-derived'},'Location','southoutside','NumColumns',3);
leg.Position(1)=0.2;
leg.Position(2)=0.01;
set(gca, 'box', 'off')
% text(0.05,0.9,'i','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('i. Southern Ocean');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;
% p1 = icons(1).Position;
% p2 = icons(2).Position;
% icons(1).Position = [0.3 p1(2) 0];
% icons(2).Position = [0.3 p2(2) 0];
% icons(3).XData = [0.05 0.2];
% icons(5).XData = [0.05 0.2];
% legend boxoff

nexttile
plot(100.*(alk-bi)./bi,'color',[0.5 0.5 0.5],'linewidth',1)
hold on
plot(100.*(alk2-bi)./bi,'color',gc,'linewidth',1)
plot(bi.*0,'k--')
xlabel('Month')
% ylabel('% difference [-]')
xlim([1 12])
% ylim([-50 5])
set(gca, 'box', 'off')
% text(0.05,0.9,'j','Units','normalized','fontweight','bold','fontsize',fsiz)  
    ttl =title('j');
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; 
    ttl.HorizontalAlignment = 'left';  
    ttl.FontSize=fsiz;

% figure
% plot(squeeze(bbp.phyto(xid,yid,:))./squeeze(bbp.tot(xid,yid,:)))

set(gcf,'color','w')

%% seasonal locations

xid1=325;
yid1=136;

xid2=195;
yid2=128;

xid3=234;
yid3=82;

xid4=178;
yid4=97;

xid5=213;%254;
yid5=16;%18;

figure; 
surface(env.lon,env.lat,log10(bbp.tot(:,:,7)')); 
shading flat
hold on
scatter(env.lon(xid1),env.lat(yid1),80,'pw','filled','markeredgecolor','k')
scatter(env.lon(xid2),env.lat(yid2),80,'pw','filled','markeredgecolor','k')
scatter(env.lon(xid3),env.lat(yid3),80,'pw','filled','markeredgecolor','k')
scatter(env.lon(xid4),env.lat(yid4),80,'pw','filled','markeredgecolor','k')
scatter(env.lon(xid5),env.lat(yid5),80,'pw','filled','markeredgecolor','k')
axis tight

shading flat

set(gcf,'color','w')

%%

%do plot fit bbp-Cphyto
idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=0;
log_on=0;
md1=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx)+biom.mixo(idx),rob_opts,plot_on,log_on);

idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=0;
log_on=0;
figure
md2=func_fit_regression_plot(bbp.tot(idx),(Chl.phyto(idx)+Chl.mixos(idx)).*150,rob_opts,plot_on,log_on);

bbp_alg=md1.int+ md1.slope.*bbp.tot;
bbpchl_alg=md2.int+ md2.slope.*bbp.tot;
biomphy=biom.phyto+biom.mixo;
%%
lat2d=repmat(env.lat',[360,1]);
figure
tiledlayout(2,6)
bb_diff_vec=zeros(7,2);
for i=1:12
imon=i;

% bbp_diff=100.*(abs(bbp_alg(:,:,imon)-biomphy(:,:,imon)))./biomphy(:,:,imon);
bbp_diff=10.^(abs(log10(bbp_alg(:,:,imon))-log10(biomphy(:,:,imon))));
% bbp_diff=100.*(5.*bbp_diff-5)./5;
bb_diff_vec(4,1)=nanmean(bbp_diff(idx_darwin_biomes.trop==1));
bb_diff_vec(3,1)=nanmean(bbp_diff(idx_darwin_biomes.oligo==1 & lat2d>=0));
bb_diff_vec(5,1)=nanmean(bbp_diff(idx_darwin_biomes.oligo==1 & lat2d<0));
bb_diff_vec(2,1)=nanmean(bbp_diff(idx_darwin_biomes.temp==1 & lat2d>=0));
bb_diff_vec(6,1)=nanmean(bbp_diff(idx_darwin_biomes.temp==1 & lat2d<0));
bb_diff_vec(1,1)=nanmean(bbp_diff(idx_darwin_biomes.polar==1 & lat2d>=0));
bb_diff_vec(7,1)=nanmean(bbp_diff(idx_darwin_biomes.polar==1 & lat2d<0));

bbp_diff=100.*((bbpchl_alg(:,:,imon)-biomphy(:,:,imon)))./biomphy(:,:,imon);
% bbp_diff=10.^(abs(log10(bbpchl_alg(:,:,imon))-log10(biomphy(:,:,imon))));
% bbp_diff=100.*(5.*bbp_diff-5)./5;
bb_diff_vec(4,2)=nanmean(bbp_diff(idx_darwin_biomes.trop==1));
bb_diff_vec(3,2)=nanmean(bbp_diff(idx_darwin_biomes.oligo==1 & lat2d>=0));
bb_diff_vec(5,2)=nanmean(bbp_diff(idx_darwin_biomes.oligo==1 & lat2d<0));
bb_diff_vec(2,2)=nanmean(bbp_diff(idx_darwin_biomes.temp==1 & lat2d>=0));
bb_diff_vec(6,2)=nanmean(bbp_diff(idx_darwin_biomes.temp==1 & lat2d<0));
bb_diff_vec(1,2)=nanmean(bbp_diff(idx_darwin_biomes.polar==1 & lat2d>=0));
bb_diff_vec(7,2)=nanmean(bbp_diff(idx_darwin_biomes.polar==1 & lat2d<0));

nexttile
barh(flip(bb_diff_vec(:,1)))
% xlim([1 6])
% xlim([0 150])
grid on
end

set(gcf,'color','w')

%%

lat2d=repmat(env.lat',[360,1]);

bb_diff_vec=zeros(7,12);
for i=1:12
imon=i;

bbp_diff=100.*(abs(bbp_alg(:,:,imon)-biomphy(:,:,imon)))./biomphy(:,:,imon);
bb_diff_vec(4,i)=nanmean(bbp_diff(idx_darwin_biomes.trop==1));
bb_diff_vec(3,i)=nanmean(bbp_diff(idx_darwin_biomes.oligo==1 & lat2d>=0));
bb_diff_vec(5,i)=nanmean(bbp_diff(idx_darwin_biomes.oligo==1 & lat2d<0));
bb_diff_vec(2,i)=nanmean(bbp_diff(idx_darwin_biomes.temp==1 & lat2d>=0));
bb_diff_vec(6,i)=nanmean(bbp_diff(idx_darwin_biomes.temp==1 & lat2d<0));
bb_diff_vec(1,i)=nanmean(bbp_diff(idx_darwin_biomes.polar==1 & lat2d>=0));
bb_diff_vec(7,i)=nanmean(bbp_diff(idx_darwin_biomes.polar==1 & lat2d<0));

end

figure
imagesc(bb_diff_vec)
hold on
for i=1:7
    for j=1:12
   text(j,i,num2str(round(bb_diff_vec(i,j)))) 
    end
end
colorbar
set(gcf,'color','w')
colormap(brewermap(100,'Oranges'))



% 
% figure
% tiledlayout(4,2)
% iord=[4,3,5,2,6,1,7];
% for i=1:7
% nexttile
% plot(bb_diff_vec(iord(i),:))
% hold on
% plot(bb_diff_vec(i,:).*0,'k--')
% if i==1
%     nexttile
% end
% end


%%
addpath('C:\Users\Camila\Desktop\Backup\Projects\3Doutputs\cmap_ocean')
bp=(bbp.phyto(:,:,7)+bbp.mixo(:,:,7))./bbp.tot(:,:,7);

fig=figure
% axesm('pcarree','MapLatLimit',[-90 90],'MapLonLimit',[-60 300])
    ax = axesm ( 'Origin',  [0 -60 0], 'MapProjection','pcarree', 'Frame', 'on',...
            'ScaleFactor', 1, 'labelrotation', 'off', 'FLineWidth', 2);    
    ax.XColor = 'white';
    ax.YColor = 'white';

% ax = worldmap('world');
% setm(ax,'MLabelParallel',-90)
% setm(ax,'MLabelLocation',90)
lon=env.lon;
lon(1)=0;
lon(end)=360;

surfacem(env.lat,lon,bp')
% shading flat
% hold on
% contourm(env.lat,lon,bp',10,'edgecolor',[0.2 0.2 0.2])

% setm(gca,'mapprojection','pcarree')
axis off
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
cmap=flip(cmocean('deep',100));
colormap(cmap)

saveFigPdf(fig, 'testfig')

% set(gcf,'renderer','opengl'); 
% saveas(gcf, 'test', 'pdf')  
% saveas(gcf,'test','epsc') 
% saveas(gcf,'test.png') 