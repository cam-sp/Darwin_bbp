function func_figure_decompose_Cphytobbp(pathname,t_vec,bbp,biom,iplk,plk_sizes,Chl,env,idx)


panord=[1:2:6,2:2:6];%order of plot tiles

xmin=-5.5;
xmax=-1.5;

figure
tiledlayout(3,2)
nexttile(5)
dscatter_2(log10(bbp.tot(idx)),log10(biom.phyto(idx)+biom.mixo(idx)))
xlabel('Log_{10}b_{bp}')
ylabel('Log_{10}C_{phyto}')
xlim([xmin xmax])
text(0.05,0.95,'c','Units','normalized','fontsize',9)

nexttile(3)
dscatter_2(log10(bbp.plk(idx)),log10(biom.phyto(idx)+biom.mixo(idx)))
xlabel('Log_{10}b_{b,plk}')
ylabel('Log_{10}C_{phyto}')
xlim([xmin xmax])
text(0.05,0.95,'b','Units','normalized','fontsize',9)

nexttile(1)
dscatter_2(log10(bbp.phyto(idx)+bbp.mixo(idx)),log10(biom.phyto(idx)+biom.mixo(idx)))
xlabel('Log_{10}b_{b,phyto}')
ylabel('Log_{10}C_{phyto}')
xlim([xmin xmax])
title('Default run')
text(0.05,0.95,'a','Units','normalized','fontsize',9)

idx_wb=5; %index of waveband
bb_or_biom=1;%bb
bbps=func_get_bbp_plk(pathname,t_vec,idx_wb,iplk,plk_sizes,bb_or_biom,1); %in m^-1
bbps.tot_sensit=bbps.phyto+bbps.mixo+bbps.zoo+bbps.bact+bbp.detr_tot;%+1.*bb_strams./(ppart.*120);
bbps.plk=bbps.phyto+bbps.mixo+bbps.zoo+bbps.bact;

% idx=Chl.idx_minChl & env.bathy_3D>500;

nexttile(6)
dscatter_2(log10(bbps.tot_sensit(idx)),log10(biom.phyto(idx)+biom.mixo(idx)))
xlabel('Log_{10}b_{bp}')
ylabel('Log_{10}C_{phyto}')
xlim([xmin xmax])
text(0.05,0.95,'f','Units','normalized','fontsize',9)

nexttile(4)
dscatter_2(log10(bbps.plk(idx)),log10(biom.phyto(idx)+biom.mixo(idx)))
xlabel('Log_{10}b_{b,plk}')
ylabel('Log_{10}C_{phyto}')
xlim([xmin xmax])
text(0.05,0.95,'e','Units','normalized','fontsize',9)

nexttile(2)
dscatter_2(log10(bbps.phyto(idx)+bbps.mixo(idx)),log10(biom.phyto(idx)+biom.mixo(idx)))
xlabel('Log_{10}b_{b,phyto}')
ylabel('Log_{10}C_{phyto}')
title('Same \sigma_{bb} for all phyto')
xlim([xmin xmax])
text(0.05,0.95,'d','Units','normalized','fontsize',9)

    cmap=viridis(100);
    colormap(cmap(20:end,:))

set(gcf,'color','w')


end