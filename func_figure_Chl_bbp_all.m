function func_figure_Chl_bbp_all(bbp,Chl,argo,sat,means)

% chl_tot_t2=chl_tot_t;
% chl_tot_t2(chl_tot_t==0)=NaN;
% chl_tot_t2(chl_tot_t<0.008)=NaN;

% idxsub=datasample(1:length(chl_tot_t2(:)),10000);

idx= isfinite(bbp.tot700(:)) & isfinite(Chl.tot(:));
chl_tot_t2=Chl.tot(idx);
bbpart_darwin4502=bbp.tot700(idx);

xmin=-5;
xmax=-1;
ymin=-3;
ymax=2.5;

% bval=0.0003;

% alk=bbp_argo(Chl_argo>0.2);%-0.00035;
% alk(alk<0)=NaN;
% md=fitlm(log10(Chl_argo(Chl_argo>0.2)),log10(alk));
% mm=(md.Coefficients.Estimate(1)+linspace(xmin,xmax,100).*md.Coefficients.Estimate(2));

bbp_argo2=argo.bbp;%-0.00035;
% bbp_argo2(bbp_argo2<0)=NaN;

% mdtotr=fitlm(bbp_argo,Chl_argo,'RobustOpts','on');
% mold=mdtotr.Coefficients.Estimate(1)+bbp_argo.*mdtotr.Coefficients.Estimate(2);
% Rsq=func_Rsq(Chl_argo,mold)
% [~,~,~, logMAE]=func_calc_bias(mold,Chl_argo)

% b = gmregress(sqrt(bbp_argo),sqrt(Chl_argo));
% mold=b(1)+sqrt(bbp_argo).*b(2);
% Rsq2=func_Rsq(sqrt(Chl_argo),mold)
% [~,~,~, logMAE]=func_calc_bias(mold,sqrt(Chl_argo))

xx=(argo.bbp);
yy=(argo.Chl);
mdtotr=fitlm(xx,yy,'RobustOpts','on');
% mold=mdtotr.Coefficients.Estimate(1)+xx.*mdtotr.Coefficients.Estimate(2);
% b=gmregress(xx,yy);
% mold=b(1)+xx.*b(2);
% Rsq=func_Rsq(Chl_argo,mold)
% [~,~,~, logMAE]=func_calc_bias(mold,Chl_argo)



% Rsq=func_Rsq(bbp_argo2,mdtot.Coefficients.Estimate(1)+mdtot.Coefficients.Estimate(2).*Chl_argo)


figure
tiledlayout(2,3)


nexttile%(3,[2,1])
[hAxes, h3]=dscatter_2(log10(bbpart_darwin4502),log10(chl_tot_t2));
hold on
md11=fitlm(bbpart_darwin4502,chl_tot_t2,'RobustOpts','on');
mim2=md11.Coefficients.Estimate(1)+logspace(-4,-1,500).*md11.Coefficients.Estimate(2);
mim2(mim2<0)=NaN;
plot(linspace(-4,-1,500),log10(mim2),'k','linewidth',0.8);
Rsq=func_Rsq(bbpart_darwin4502,md11.Coefficients.Estimate(1)+md11.Coefficients.Estimate(2).*chl_tot_t2);
text(0.05,0.95,strcat(['Chl = ',num2str(md11.Coefficients.Estimate(2),'%.2f'),' b_{bp}(\lambda) '...
                ,num2str(md11.Coefficients.Estimate(1),'%.2f')]),'Units','normalized','fontsize',8)
text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8)        
h3.Marker='o';
h3.SizeData=5;
xticks([-4:2:4]);
xticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'});
yticks([-4:2:2]);
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'});
ylim([ymin ymax])
xlim([-4 xmax])
ylabel('Chl [mg m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
title('EcoMITgcm')
text(0.9,0.08,'a','Units','normalized','fontsize',9,'fontweight','bold') 

Chl_argo=argo.Chl;
nexttile%(1,[2,1])
idx= ~isnan(Chl_argo) & ~isnan(bbp_argo2);
[hAxes, h]=dscatter_2(log10(bbp_argo2(idx)),log10(Chl_argo(idx)))
title('BGC-Argo')
hold on
text(0.05,0.95,strcat(['Chl = ',num2str(mdtotr.Coefficients.Estimate(2),'%.2f'),' b_{bp}(\lambda) '...
                ,num2str(mdtotr.Coefficients.Estimate(1),'%.2f')]),'Units','normalized','fontsize',8)     
h.Marker='o';
h.SizeData=5;
xticks([-4:2:4]);
xticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'});
yticks([-4:2:2]);
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'});
ylim([ymin ymax])
xlim([xmin xmax])
ylabel('Chl [mg m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
mm=mdtotr.Coefficients.Estimate(1)+mdtotr.Coefficients.Estimate(2).*10.^linspace(xmin,xmax,100);
mm(mm<0)=NaN;
plot(linspace(xmin,xmax,100),...
    log10(mm),...
    'color','k','linewidth',0.8);%[0.5 0.5 0.5])
Rsq=func_Rsq(Chl_argo(idx),mdtotr.Coefficients.Estimate(1)+mdtotr.Coefficients.Estimate(2).*bbp_argo2(idx));
text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8)   
% plot(linspace(-4,2,10).*0+log10(6.3286e-04),...
%     linspace(-4,2,10),...
%     '--','color',[0.5 0.5 0.5],'linewidth',0.8);%[0.5 0.5 0.5])
text(0.9,0.08,'b','Units','normalized','fontsize',9,'fontweight','bold') 
set(gcf,'color','w');

cmap=(inferno(100));
colormap(cmap(30:end,:))



addpath C:\Users\Camila\Desktop\Backup\Projects\3Doutputs\quantreg


veclen=1:length(sat.bbp(:));
idx=datasample(veclen,100000);
bbp_sub=sat.bbp(idx);
chl_sub=sat.Chl(idx);
idx = isfinite(chl_sub(:)) & isfinite(bbp_sub(:));
md1=fitlm(bbp_sub(:),chl_sub(:),'RobustOpts','on');
% md2=fitlm(bbp_sub(:),chl_sub(:));
% b=quantreg(bbp_sub(idx)',chl_sub(idx)',0.5);

nexttile%(2,[2,1])
[hAxes, h2]=dscatter_2(log10(bbp_sub(idx)'),log10(chl_sub(idx)'))
hold on
mim=md1.Coefficients.Estimate(1)+logspace(-4,-1,500).*md1.Coefficients.Estimate(2);
mim(mim<0)=NaN;
plot(linspace(-4,-1,500),log10(mim),'k','linewidth',0.8);
h2.Marker='o';
h2.SizeData=5;
xticks([-4:2:4]);
xticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'});
    yticks([-4:2:2]);
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'});
% dscatter(log10(Chl_argo(idx)),log10(bbp_argo(idx)))
ylim([ymin ymax])
% xlim([xmin xmax])
ylabel('Chl [mg m^{-3}]')
xlabel('b_{bp}(443) [m^{-1}]')
% plot(linspace(-4,2,10).*0+log10(1.5e-3),...
%     linspace(-4,2,10),...
%     '--','color',[0.5 0.5 0.5],'linewidth',0.8);%[0.5 0.5 0.5])
text(0.05,0.95,strcat(['Chl = ',num2str(md1.Coefficients.Estimate(2),'%.2f'),' b_{bp}(\lambda) '...
                ,num2str(md1.Coefficients.Estimate(1),'%.2f')]),'Units','normalized','fontsize',8)
title('Remote sensing')
text(0.9,0.08,'c','Units','normalized','fontsize',9,'fontweight','bold') 
Rsq=func_Rsq(bbp_sub(idx),md1.Coefficients.Estimate(1)+md1.Coefficients.Estimate(2).*chl_sub(idx));
text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8) 
% yy=bbp_sub(idx);
% yyhat=md1.Coefficients.Estimate(1)+md1.Coefficients.Estimate(2).*chl_sub(idx);
% yyhat(yyhat<0)=NaN;
% Rsq = 1 - sum((yy - yyhat).^2)/sum((yy - mean(yy)).^2);
% yymean=mean(yy)
% yyhatmean=mean(yyhat)
% ((sum((yy-yymean).*(yyhat-yyhatmean)))./(sum((yy-yymean.^2)).*sum((yyhat-yyhatmean.^2))).^(1/2)).^2






cmap=brewermap(4,'spectral');
msiz=8;
trm=0.7;


nexttile%(9,[2,1])
% scatter(bbp_argo_darwin_mean,chl_argo_darwin_mean,5,'k','filled')
hold on
% scatter(bbp_sb,chl_sb,5,'filled','markerfacecolor',[0.5 0.5 0.5])
% scatter(bbp_argo_mean,chl_argo_mean,5,'filled','markerfacecolor',[0.5 0.5 0.5])
scatter(means.darwin.bbp(means.idx.ibiom==1),means.darwin.Chl(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:))
scatter(means.darwin.bbp(means.idx.ibiom==2),means.darwin.Chl(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:))
scatter(means.darwin.bbp(means.idx.ibiom==3),means.darwin.Chl(means.idx.ibiom==3),msiz,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:))
scatter(means.darwin.bbp(means.idx.ibiom==4),means.darwin.Chl(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:))
plot(logspace(-4,-1,500),mim2,'k','linewidth',0.8);
mdr3=fitlm(means.darwin.bbp,means.darwin.Chl,'RobustOpts','on');
mimr3=mdr3.Coefficients.Estimate(1)+mdr3.Coefficients.Estimate(2).*logspace(-4,-1,500);
plot(logspace(-4,-1,500),mimr3,'k--','linewidth',0.8)
text(0.05,0.95,strcat(['Chl = ',num2str(mdr3.Coefficients.Estimate(2),'%.2f'),' b_{bp}(\lambda) '...
                ,num2str(mdr3.Coefficients.Estimate(1),'%.2f')]),'Units','normalized','fontsize',8)
Rsq=func_Rsq(means.darwin.Chl,mdr3.Coefficients.Estimate(1)+mdr3.Coefficients.Estimate(2).*means.darwin.bbp);
text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8)   
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e-2 1e1])
xlim([1e-4 1e-2])
ylabel('Chl [mg m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
text(0.9,0.08,'d','Units','normalized','fontsize',9,'fontweight','bold') 




nexttile%(7,[2,1])

% figure
% scatter(bbp_argo_mean,chl_argo_mean,5,'k','filled')
hold on
p1=scatter(means.argo.bbp(means.idx.ibiom==1),means.argo.Chl(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:));
p2=scatter(means.argo.bbp(means.idx.ibiom==2),means.argo.Chl(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:));
p3=scatter(means.argo.bbp(means.idx.ibiom==3),means.argo.Chl(means.idx.ibiom==3),msiz,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:));
p4=scatter(means.argo.bbp(means.idx.ibiom==4),means.argo.Chl(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:));
mdr1=fitlm(means.argo.bbp,means.argo.Chl,'RobustOpts','on');
mimr1=mdr1.Coefficients.Estimate(1)+mdr1.Coefficients.Estimate(2).*logspace(-4,log10(max(means.argo.Chl)),500);
plot(logspace(-4,log10(max(means.argo.Chl)),500),mimr1,'k--','linewidth',0.8)
set(gca,'yscale','log')
set(gca,'xscale','log')
text(0.05,0.95,strcat(['Chl = ',num2str(mdr1.Coefficients.Estimate(2),'%.2f'),' b_{bp}(\lambda) '...
                ,num2str(mdr1.Coefficients.Estimate(1),'%.2f')]),'Units','normalized','fontsize',8)
Rsq=func_Rsq(means.argo.Chl,mdr1.Coefficients.Estimate(1)+mdr1.Coefficients.Estimate(2).*means.argo.bbp);
text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8) 
ylim([1e-2 1e1])
plot(logspace(xmin,xmax,100),...
    (mm),...
    'color','k','linewidth',0.8);%[0.5 0.5 0.5])
ylabel('Chl [mg m^{-3}]')
xlabel('b_{bp}(700) [m^{-1}]')
legend([p1,p2,p3,p4],'Tropical','Oligotrophic','Temperate','Sub-polar and polar','NumColumns',4,'Location','southoutside')
text(0.9,0.08,'e','Units','normalized','fontsize',9,'fontweight','bold') 


nexttile%(8,[2,1])
% scatter(bbp_sb,chl_sb,5,'k','filled')
hold on
p1=scatter(means.sat.bbp(means.idx.ibiom==1),means.sat.Chl(means.idx.ibiom==1),msiz,'^','filled','markerfacecolor',cmap(1,:),'markerfacealpha',trm,'markeredgecolor',cmap(1,:))
p2=scatter(means.sat.bbp(means.idx.ibiom==2),means.sat.Chl(means.idx.ibiom==2),msiz-3,'filled','markerfacecolor',cmap(2,:),'markerfacealpha',trm,'markeredgecolor',cmap(2,:))
p3=scatter(means.sat.bbp(means.idx.ibiom==3),means.sat.Chl(means.idx.ibiom==3),msiz,'s','filled','markerfacecolor',cmap(3,:),'markerfacealpha',trm,'markeredgecolor',cmap(3,:))
p4=scatter(means.sat.bbp(means.idx.ibiom==4),means.sat.Chl(means.idx.ibiom==4),msiz,'d','filled','markerfacecolor',cmap(4,:),'markerfacealpha',trm,'markeredgecolor',cmap(4,:))
% scatter(bbp_sb(ireg_sb==4 & ilat_sb==2),chl_sb(ireg_sb==4 & ilat_sb==2),5,'c','filled')
plot(logspace(-4,-1,500),(mim),'k','linewidth',0.8);
mdr2=fitlm(means.sat.bbp,means.sat.Chl,'RobustOpts','on');
mimr2=mdr2.Coefficients.Estimate(1)+mdr2.Coefficients.Estimate(2).*logspace(-4,-1,500);
plot(logspace(-4,-1,500),mimr2,'k--','linewidth',0.8)
set(gca,'yscale','log')
set(gca,'xscale','log')
text(0.05,0.95,strcat(['Chl = ',num2str(mdr2.Coefficients.Estimate(2),'%.2f'),' b_{bp}(\lambda) '...
                ,num2str(mdr2.Coefficients.Estimate(1),'%.2f')]),'Units','normalized','fontsize',8)
Rsq=func_Rsq(means.sat.Chl,mdr2.Coefficients.Estimate(1)+mdr2.Coefficients.Estimate(2).*means.sat.bbp);
text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8)  
ylabel('Chl [mg m^{-3}]')
xlabel('b_{bp}(443) [m^{-1}]')
% xlim([1e-3 0.02])
ylim([1e-2 1e1])
xlim([1e-3 1e-2])
text(0.9,0.08,'f','Units','normalized','fontsize',9,'fontweight','bold') 



end