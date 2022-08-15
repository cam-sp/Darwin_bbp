function func_figure_data_comparison_2(means)

msiz=10;
cmap=brewermap(4,'Spectral');
trm=0.7;
rg_color=[224, 222, 222]./255;

figure
tiledlayout(5,4,'TileSpacing','Compact','Padding','Compact')
nexttile(1,[2,2])
% scatter(chl_argo_darwin_mean,chl_argo_mean,msiz,'k','filled')
hold on
% for i=4%1:4
errorbar(means.argo.Chl,means.darwin.Chl,means.argo.Chl_std,'o','markersize',0.1,'color',rg_color)
errorbar(means.argo.Chl,means.darwin.Chl,means.darwin.Chl_std,'horizontal','o','markersize',0.1,'color',rg_color)
% end
h1=scatter(means.argo.Chl(means.idx.ibiom==1),means.darwin.Chl(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:))
h2=scatter(means.argo.Chl(means.idx.ibiom==2),means.darwin.Chl(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:))
h3=scatter(means.argo.Chl(means.idx.ibiom==3),means.darwin.Chl(means.idx.ibiom==3),msiz+5,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:))
h4=scatter(means.argo.Chl(means.idx.ibiom==4),means.darwin.Chl(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:))

plot(logspace(-3,1),logspace(-3,1),'k--')
% legend([h1,h2,h3,h4],{'Tropical','Oligotrophic','Temperate','Polar'},'Location','northwest')
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('Chl Darwin model [mg m^{-3}]')
xlabel('Chl BGC-Argo [mg m^{-3}]')
md=fitlm(log10(means.darwin.Chl),log10(means.argo.Chl));
[~,~,~, logMAE]=func_calc_bias(means.darwin.Chl,means.argo.Chl);
text(0.05,0.85,strcat(['R^2 = ',num2str(md.Rsquared.Adjusted,'%.2f')]),'Units','normalized','fontsize',8)   
text(0.05,0.71,strcat(['MAE = ',num2str(logMAE,'%.2f')]),'Units','normalized','fontsize',8)   
text(0.05,1,'a','Units','normalized','fontweight','bold')   
xlim([1e-3 1e1])
ylim([1e-3 1e1])

nexttile(3,[2,2])
% scatter(chl_argo_darwin_mean,chl_argo_mean,msiz,'k','filled')
hold on
errorbar(means.argo.bbp,means.darwin.bbp,means.argo.bbp_std,'o','markersize',0.1,'color',rg_color)
errorbar(means.argo.bbp,means.darwin.bbp,means.darwin.bbp_std,'horizontal','o','markersize',0.1,'color',rg_color)

h1=scatter(means.argo.bbp(means.idx.ibiom==1),means.darwin.bbp(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:))
h2=scatter(means.argo.bbp(means.idx.ibiom==2),means.darwin.bbp(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:))
h3=scatter(means.argo.bbp(means.idx.ibiom==3),means.darwin.bbp(means.idx.ibiom==3),msiz+5,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:))
h4=scatter(means.argo.bbp(means.idx.ibiom==4),means.darwin.bbp(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:))
% scatter(chl_argo_darwin_mean(ibas_argo_mean==2),chl_argo_mean(ibas_argo_mean==2),8,'r','filled')
plot(logspace(-4,-2),logspace(-4,-2),'k--')
% [~, objh] = legend([h1,h2,h3,h4],{'Tropical','Oligotrophic','Temperate','Sub-polar and Polar'})
% [~, objh] = legend({'one','two'}); % Instead of "h_legend" use "[~, objh]"
% objhl = findobj(objh, 'type', 'patch');
% set(objhl, 'Markersize', 5);
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('b_{bp}(700) Darwin model [m^{-1}]')
xlabel('b_{bp}(700) BGC-Argo [m^{-1}]')
% title('Comparison Darwin with Argo b_{bp}')
set(gcf,'color','w');
md=fitlm(log10(means.darwin.bbp),log10(means.argo.bbp));
[~, ~,~, logMAE2]=func_calc_bias(means.darwin.bbp,means.argo.bbp);
text(0.05,0.85,strcat(['R^2 = ',num2str(md.Rsquared.Adjusted,'%.2f')]),'Units','normalized','fontsize',8)   
text(0.05,1,'b','Units','normalized','fontweight','bold')   
text(0.05,0.71,strcat(['MAE = ',num2str(logMAE2,'%.2f')]),'Units','normalized','fontsize',8)   
ylim([2e-4 5e-3])
xlim([2e-4 5e-3])
yticks([5e-4 1e-3 5e-3])
xticks([5e-4 1e-3 5e-3])

nexttile(9,[2,2])
% scatter(chl_argo_darwin_mean,chl_argo_mean,msiz,'k','filled')
hold on
errorbar(means.argoSat.Chl,means.sat.Chl,means.argoSat.Chl_std,'o','markersize',0.1,'color',rg_color)
errorbar(means.argoSat.Chl,means.sat.Chl,means.argoSat.Chl_std,'horizontal','o','markersize',0.1,'color',rg_color)

h1=scatter(means.argoSat.Chl(means.idx.ibiom==1),means.sat.Chl(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:))
h2=scatter(means.argoSat.Chl(means.idx.ibiom==2),means.sat.Chl(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:))
h3=scatter(means.argoSat.Chl(means.idx.ibiom==3),means.sat.Chl(means.idx.ibiom==3),msiz+5,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:))
h4=scatter(means.argoSat.Chl(means.idx.ibiom==4),means.sat.Chl(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:))
plot(logspace(-4,1),logspace(-4,1),'k--')
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('Chl MODIS-GIOP [mg m^{-3}]')
xlabel('Chl BGC-Argo [mg m^{-3}]')
% title('Comparison Darwin with Argo Chl')
md=fitlm(log10(means.sat.Chl),log10(means.argoSat.Chl));
[~,~,~, logMAE]=func_calc_bias(means.sat.Chl,means.argoSat.Chl);
text(0.05,0.85,strcat(['R^2 = ',num2str(md.Rsquared.Adjusted,'%.2f')]),'Units','normalized','fontsize',8)   
text(0.05,0.71,strcat(['MAE = ',num2str(logMAE,'%.2f')]),'Units','normalized','fontsize',8)   
text(0.05,1,'c','Units','normalized','fontweight','bold')   
ylim([1e-3 1e1])
xlim([1e-3 1e1])



nexttile(11,[2,2])
% scatter(chl_argo_darwin_mean,chl_argo_mean,msiz,'k','filled')
hold on
errorbar(means.argoSat.bbp,means.sat.bbp700,means.argoSat.bbp_std,'o','markersize',0.1,'color',rg_color)
errorbar(means.argoSat.bbp,means.sat.bbp700,means.sat.bbp700_std,'horizontal','o','markersize',0.1,'color',rg_color)

h1=scatter(means.argoSat.bbp(means.idx.ibiom==1),means.sat.bbp700(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:))
h2=scatter(means.argoSat.bbp(means.idx.ibiom==2),means.sat.bbp700(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:))
h3=scatter(means.argoSat.bbp(means.idx.ibiom==3),means.sat.bbp700(means.idx.ibiom==3),msiz+5,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:))
h4=scatter(means.argoSat.bbp(means.idx.ibiom==4),means.sat.bbp700(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:))
plot(logspace(-4,-2),logspace(-4,-2),'k--')
set(gca,'yscale','log')
set(gca,'xscale','log')
ylabel('b_{bp}(700) MODIS-GIOP [m^{-1}]')
xlabel('b_{bp}(700) BGC-Argo [m^{-1}]')
% title('Comparison Darwin with Argo b_{bp}')
set(gcf,'color','w');
md=fitlm(log10(means.sat.bbp700),log10(means.argoSat.bbp));
[~, ~,~, logMAE2]=func_calc_bias(means.sat.bbp700,means.argoSat.bbp);
text(0.05,0.85,strcat(['R^2 = ',num2str(md.Rsquared.Adjusted,'%.2f')]),'Units','normalized','fontsize',8)   
text(0.05,1,'d','Units','normalized','fontweight','bold')   
text(0.05,0.71,strcat(['MAE = ',num2str(logMAE2,'%.2f')]),'Units','normalized','fontsize',8)   
ylim([2e-4 5e-3])
xlim([2e-4 5e-3])
yticks([5e-4 1e-3 5e-3])
xticks([5e-4 1e-3 5e-3])
leg=legend([h1,h2,h3,h4],{'Tropical','Oligotrophic','Temperate','Sub-polar and Polar'},'Location','southeast','NumColumns',4,'Location','southoutside')
leg.Position=[0.12 0.08 0.8057 0.0419];




end