%% Script for "Assessing the potential for backscattering as a proxu fro phytoplankton biomass"
%  Submitted to Global Biogeochemical Cycles

% This file loads all model outputs and generates the figures that appear
% in the paper.

% Camila Serra-Pompei 26/01/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
clc

%load paths
addpath('Functions/')
addpath('Colormaps/')
addpath('Other_extra_files/')

%If you are running this for the firt time, remember to unzip the
%Model_outputs folder. Can be done by running "unzip(Model_outputs)"
pathname='Model_outputs/';

%get plk indexes and size info
[iplk, plk_sizes]=func_get_plk_info();

%get darwin environmental variables
t_vec=29040:240:31680;% vector used in the naming of the files
env=func_get_environmenatl_data(pathname,t_vec);

wb_ranges=400-12.5:25:700+12.5; %waveband limits
wb=400:25:700; %mean wavelength in each waveband interval
idx_wb=3;%7; %index of waveband
bb_or_biom='bb'; %get bbp or biomass
bbp=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom); %in m^-1

%get total bbp (in m^-1)
filename='iops.';
var='bb'; %variable to get
bbp.tot=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only
bbp_tot_700=func_get_var_ncfile(pathname,filename,t_vec,var,1,13); %for different wavelengths 
bbp_tot_475=func_get_var_ncfile(pathname,filename,t_vec,var,1,4);

%get bbp by detritus only (in m^-1)
filename='iops.';
var='bbprt';
bbp.detr_tot=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of detrital particles only

%get bbp by plankton only (in m^-1)
filename='iops.';
var='bbplk';
bbp.plk=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of plk only

bb_or_biom='biomass';
biom=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom); %in mgC m^-3
biom.tot = biom.phyto + biom.mixo + biom.zoo + biom.bact;

% get Chl
minchl=1e-3;%minimum chl threshold to account in plots
Chl=func_get_chl(pathname,t_vec,minchl); %in mgChl m^-3

% other bbps at different lambdas
idx_wb=13;
filename='iops.';
var='bb';
bbp.tot700=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only

idx_wb=3;
filename='iops.';
var='bb';
bbp.tot450=func_get_var_ncfile(pathname,filename,t_vec,var,1,idx_wb); %bb of particles only

%% some plots

% do plot contribution by each plankton to bbp plot
% figure 8 in the paper
func_figure_bbpcontribution_map(bbp.tot,bbp.phyto+bbp.mixo,bbp.bact,bbp.detr_tot,bbp.zoo,env.lat,env.lon_map)
% set(gcf,'Renderer','Painter')
% print(gcf, '-depsc2', 'Figures_paper_bbp/ffig8_contributions.eps');

% do plot decompose Cphyto-bbp relationship and effects of phyto sizes on bbp
%figure 9 in the paper
idx=Chl.idx_minChl & env.bathy_3D>500;
func_figure_decompose_Cphytobbp(pathname,t_vec,bbp,biom,iplk,plk_sizes,Chl,env,idx)
% set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig9_steps.eps


%% get Argo + satellite data

% get Longhurst regions
[idx_darwin_biomes,idx_darwin_basins,...
    idx_argo_biomes,idx_argo_basins,...
    idx_sat_biomes,idx_sat_basins...
    ]=func_Longhurst_regions();

% get sat data
sat=func_get_sat();

%get argo data and matchings with Darwin and with sat
[argo, ArgoDarwin, ArgoSat]=func_get_Argo(bbp,biom,Chl,env,sat,idx_argo_biomes);

%plot argo data and cumsum plots
%figure 3
func_figure_argo_data(argo,idx_argo_biomes)
% set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig3_argo_bbcrit.eps

%plot same as above but for sat
%SI figure
func_figure_sat_data(sat,idx_sat_biomes)
nexttile(11)
xlabel('Log_{10}b_{bp}(443)')
nexttile(12)
xlabel('Log_{10}b_{bp}(443)')

%% get mean by regions

means=func_get_mean_regions(bbp,biom,argo,ArgoDarwin,sat,ArgoSat,...
    idx_argo_biomes,idx_argo_basins,idx_sat_biomes,idx_sat_basins);

%figure 4
func_figure_data_comparison(means)
% set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig4_darwin_argo_comparison.eps


%%

%chl conversion factor (factor Q in paper)
%figure 5
chlconv=110; %conversion factor Q in paper
[md_cphyto, md]=func_figure_model_regressions(Chl, env, bbp, biom, chlconv);
% set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig5_cbbp_darwin.eps

%%

%Obtain regression parameters fro chl-based fit
rob_opts=1; %robust option for linear fitting are "on"
idx=Chl.idx_minChl & env.bathy_3D>500; %regions deeper than 500m and above chl threshold
plot_on=0; %just get the model, not the plot
log_on=0; %in linear scale
md_chl=func_fit_regression_plot(bbp.tot(idx),Chl.tot(idx).*chlconv,rob_opts,plot_on,log_on,1);

% do percent difference plots
%figure 6
idx=Chl.idx_minChl & env.bathy_3D>500;
func_figure_percentdiff_maps_Chl(biom.phyto+biom.mixo,bbp.tot,idx,md_cphyto,md_chl,env.bathy,env.lat,env.lon_map)
%  set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig6_maps.eps

%% seasonal

%Seasonality figure and locations
%figure 7
func_figure_seasonal(md_cphyto,md_chl,bbp,biom,env)

%% plot comparison of algorithms

%figure 2
func_plot_algorithm_comparison()
% set(gcf,'Renderer','Painter')
% print -depsc Figures_paper_bbp/ffig2_comparison.eps

%% plot argo data

figure
rob_opts=1;
plot_on=1;
log_on=0;
func_fit_regression_plot_argo(argo.bbp,argo.Chl.*130,rob_opts,plot_on,log_on,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra figures for the SI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do model without mixotrophs

chlconv=130;

x0=0;
y0=0;
width=10;
height=8;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

idx=Chl.idx_minChl & env.bathy_3D>500;
rob_opts=1;
plot_on=1;
log_on=0;
md_cmixo=func_fit_regression_plot(bbp.tot(idx),biom.phyto(idx),rob_opts,plot_on,log_on);
ylabel('C_{phyto} [mgC m^{-3}]')
xlabel('b_{bp}(450) [m^{-1}]')
% ylim([-1.5 3.5])
xlim([-4 -1])
xticks([-4:1:-1]);
xticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
yticks([-2:1:4]);
yticklabels({'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});
xlim([-3.5 -1.5]) 
text(0.9,0.05,'a','Units','normalized','fontweight','bold')   


%% estimate area below bbcrit

for i=1:12

bcri=bbp_tot_475(:,:,i);
bcr=bcri;
bcr(bcri>1e-3)=0;
bcr(bcri<=1e-3)=1;
bcr1(i)=nansum(nansum(bcr.*env.bins_area))./nansum(nansum(env.bins_area));

biom_bbcr=biom.phyto(:,:,i) + biom.mixo(:,:,i);
bcr=biom_bbcr;
bcr(bcri>1e-3)=0;
bcrBiom_1(i)=nansum(nansum(bcr.*env.bins_area))./nansum(nansum(biom_bbcr.*env.bins_area));

bcri=bbp_tot_475(:,:,i);
bcr=bcri;
bcr(bcri>7e-4)=0;
bcr(bcri<=7e-4)=1;
bcr2(i)=nansum(nansum(bcr.*env.bins_area))./nansum(nansum(env.bins_area));

biom_bbcr=biom.phyto(:,:,i) + biom.mixo(:,:,i);
bcr=biom_bbcr;
bcr(bcri>1e-3)=0;
bcrBiom_2(i)=nansum(nansum(bcr.*env.bins_area))./nansum(nansum(biom_bbcr.*env.bins_area));

end

%% make figures for averaging as in behrenfeld et al 2005

func_figure_average_as_beh(Chl,bbp,sat, idx_darwin_basins, idx_sat_basins, env)
