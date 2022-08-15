function md_out=func_fit_regression_plot(xvar,yvar,rob_opts,plot_on,log_on,label_chl)
%This function fits a linear regression model to the inputs

idx=isfinite(xvar) & isfinite(yvar);
xvar=xvar(idx);
yvar=yvar(idx);

if rob_opts==0
    md=fitlm(xvar,yvar);
elseif rob_opts==1
    md=fitlm(xvar,yvar,'RobustOpts','on');
end

md_out.int=md.Coefficients.Estimate(1);
md_out.slope=md.Coefficients.Estimate(2);

mold=(md_out.int+xvar.*md_out.slope);
Rsq=func_Rsq(yvar,mold);
if log_on==0
    [~,~,~, logMAE]=func_calc_bias(mold,yvar);
elseif log_on==1
    [~,~,~, logMAE]=func_calc_bias(10.^mold,10.^yvar);    
end

md_out.Rsq=Rsq;
md_out.MAElog=logMAE;

if plot_on==1 && log_on==0 
    
%     figure
    [hAxes, h]=dscatter_2(log10(xvar),log10(yvar));
    hold on
    h.Marker='o';
    h.MarkerEdgeColor='flat';
    h.SizeData=5;

    xveclog=logspace(log10(min(xvar)),log10(max(xvar)),500);
%     xveclin=linspace(min(log10(xvar)),max(log10(xvar)),500);
    mdd=md.Coefficients.Estimate(1)+md.Coefficients.Estimate(2).*xveclog;
    mdd(mdd<0)=NaN;
    plot(log10(xveclog),log10(mdd),'k','linewidth',1)
    caxis([0 1])
    set(gcf,'color','w')
    cmap=viridis(100);
    colormap(cmap(20:end,:))

%     set(findall(gcf,'-property','FontSize'),'FontSize',9)
    if nargin==6
    text(0.05,0.95,strcat(['C_{phyto} = ',num2str(md_out.slope,'%.1f'),' b_{bp}(\lambda) '...
        ,num2str(md_out.int,'%.1f')]),'Units','normalized','fontsize',8)
    else
    text(0.05,0.95,strcat(['Chl = ',num2str(md_out.slope,'%.1f'),' b_{bp}(\lambda) '...
        ,num2str(md_out.int,'%.1f')]),'Units','normalized','fontsize',8)   
    end
    text(0.05,0.85,strcat(['R^2 = ',num2str(Rsq,'%.2f')]),'Units','normalized','fontsize',8)   
    text(0.05,0.75,strcat(['MAE = ',num2str(logMAE,'%.2f')]),'Units','normalized','fontsize',8) 
end
end
