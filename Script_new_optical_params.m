%% Script for "Assessing the potential for backscattering as a proxu fro phytoplankton biomass"
%  Submitted to Global Biogeochemical Cycles

% This file generates the new optical parameters used in the paper and
% loaded into the MITgcmBgc model. 
% It also geerates figure 1 of the paper.

% Camila Serra-Pompei 26/01/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
clc

addpath('Functions/')
addpath('Colormaps/')
addpath('Other_extra_files/')

opts_init=readmatrix('optics_plankton_new3.txt');
b=zeros(13,12);
bb_optf=zeros(13,12);
idx_start=2:15:176;
idx_end=idx_start+12;
for i=1:12
    b(:,i)=opts_init(idx_start(i):idx_end(i),4);
    bb_optf(:,i)=opts_init(idx_start(i):idx_end(i),5);
end
diam=opts_init(1:15:end,4); %micrometers
vol=10^-0.281.*diam.^3; %micrometers cube
C=(0.109.*vol.^0.991).*1e-9; %mgC Montagnes 

[~,~,bb1]=func_get_bb_bbratio(diam); % also provides bb using Vaillancourt (bb1) in m2 cell-1

%use fudge factors for functional types
bb1(6)= bb1(6).*5;  %coccos
bb1(9)= bb1(9).*2;  %mixos
bb1(11)= bb1(11)./5;  %zoos

bb2=bb1./C; % convert to units of m2 mgC-1

% bb=bb2'.*b(5,:)./b;
bb=bb2'.*b./b(5,:); %propagate to other wavelengths assuming same spectrum as bt

bb_tricho=b(:,7).*0.004; %assumes a b:bb ratio from littertature
bb(:,7)=bb_tricho;


[iplk, plk_sizes]=func_get_plk_info();

bbphy_mgC_type=bb; % bb for each references plankton type

idx_pro=1;
idx_syn=2;
idx_smeuk=3:4;
idx_diazo=10:13;
idx_tricho=14;
leg={'idx_smeuk', 'syn', 'pro','pro', 'diatoms', 'coccos', 'tricho', 'diazo', 'mixos','mean', 'zoos','bact'};
id={idx_smeuk, idx_syn, idx_pro, iplk.ft_diatoms, iplk.ft_coccos, idx_tricho, idx_diazo, iplk.ft_mixos, iplk.ft_zoos,iplk.bact};


%generate vector saying which optical group each plk is
iops=zeros(1,50);
st=1;
for i=1:12
    if i~=4 && i~=10
        iops(id{st})=i;
        st=st+1;
    end
end

darwin_bbbSlope=2.387;
bbphy_mgC=zeros(50,13);
for i=1:12
    for l=1:13
        bbphy_cell_type = bbphy_mgC_type(l,i).*C(i);
        dmratio=plk_sizes.ESD(iops==i)./diam(i);
        bbphy_cell_type = bbphy_cell_type*dmratio.^darwin_bbbSlope;
        bbphy_mgC(iops==i,l) = bbphy_cell_type./(plk_sizes.C_opts(iops==i));
    end
end

leg={'sm euk', 'syn', 'pro', 'diatoms', 'coccos', 'tricho', 'diazo', 'mixos', 'zoos','bact'};

id_wb=5;

figure
st=1;
for i=1:12
    if i~=4 && i~=10
%         iops(id{st})=i;
        plot(plk_sizes.ESD(id{st}),bbphy_mgC(id{st},id_wb).*plk_sizes.C_opts(id{st})','d-')
        st=st+1;
    end
    hold on
end
set(gca,'yscale','log')
set(gca,'xscale','log')
legend(leg)

bb_wb_new=bbphy_mgC';
save('bb_wb_new.mat','bb_wb_new');

%% total scattering

bt_data=importdata('Model_outputs/p-ini-char-btspec.dat');
bt_wb=bt_data.data;

figure
st=1;
for i=1:12
    if i~=4 && i~=10
%         iops(id{st})=i;
        plot(plk_sizes.ESD(id{st}),bt_wb(id_wb,id{st}),'d-')
        st=st+1;
    end
    hold on
end
set(gca,'yscale','log')
set(gca,'xscale','log')
legend(leg)

%% good plots paper total scattering

id_wb=7;
leg={'Pico-eukaryotes', 'Synechococcus', 'Prochlorococcus', 'Diatoms', 'Coccocolithophres',...
    'Trichodesmium', 'Diazotrophs', 'Mixotrophs', 'Zooplankton','H. bacteria'};
markst={'-*','-o','-d','-^','-v','->','-s','-p','-<','-x'};

cmap=flip(brewermap(10,'Spectral'));

x0=0;
y0=0;
width=14;
height=7;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(4,2)
nexttile(1,[3,1])
st=1;
for i=1:12
    if i~=4 && i~=10
        idx=id{st};
        plot(400:25:700,bt_wb(:,idx(1)).*plk_sizes.C_opts(idx(1))',markst{st},'color',cmap(st,:) ...
        ,'markerfacecolor',cmap(st,:),'markersize',4,'markeredgecolor',[0.5 0.5 0.5])
        st=st+1;
    end
    hold on
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Wavelength [nm]')
ylabel('\sigma_{b}(\lambda) [m^2 particle^{-1}]')
text(0.05,0.95,'a','Units','normalized','fontsize',9,'fontweight','bold')

nexttile(2,[3,1])
st=1;
for i=1:12
    if i~=4 && i~=10
        idx=id{st};
        plot(plk_sizes.ESD(id{st}),bt_wb(id_wb,id{st}).*plk_sizes.C_opts(id{st}),markst{st},'color',cmap(st,:)...
            ,'markerfacecolor',cmap(st,:),'markersize',4,'markeredgecolor',[0.5 0.5 0.5])
        st=st+1;
    end
    hold on
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Cell diameter [\mum]')
ylabel('\sigma_{b}(550) [m^2 particle^{-1}]')
xticks([1e0, 1e1, 1e2, 1e3])
text(0.05,0.95,'b','Units','normalized','fontsize',9,'fontweight','bold')
leg2=legend(leg,'Location','southeast','NumColumns',4,'Location','southoutside');
leg2.Position(1)=0.18;

set(gcf,'color','w')
% print -depsc Figures_paper_bbp/ffig1_spec.eps

%% good plots paper bb

id_wb=7;
leg={'Pico-eukaryotes', 'Synechococcus', 'Prochlorococcus', 'Diatoms', 'Coccocolithophres',...
    'Trichodesmium', 'Diazotrophs', 'Mixotrophs', 'Zooplankton','H. bacteria'};
markst={'-*','-o','-d','-^','-v','->','-s','-p','-<','-x'};

cmap=flip(brewermap(10,'Spectral'));

x0=0;
y0=0;
width=16;
height=8;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(4,2)
nexttile(1,[3,1])
st=1;
for i=1:12
    if i~=4 && i~=10
        idx=id{st};
        plot(400:25:700,bb_wb_new(:,idx(1)).*plk_sizes.C_opts(idx(1))',markst{st},'color',cmap(st,:) ...
        ,'markerfacecolor',cmap(st,:),'markersize',4,'markeredgecolor',[0.5 0.5 0.5])
        st=st+1;
    end
    hold on
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Wavelength [nm]')
ylabel('\sigma_{bb}(\lambda) [m^2 particle^{-1}]')
text(0.05,0.95,'a','Units','normalized','fontsize',9,'fontweight','bold')

nexttile(2,[3,1])
st=1;
for i=1:12
    if i~=4 && i~=10
        idx=id{st};
        plot(plk_sizes.ESD(id{st}),bbphy_mgC(id{st},id_wb).*plk_sizes.C_opts(id{st})',markst{st},'color',cmap(st,:)...
            ,'markerfacecolor',cmap(st,:),'markersize',4,'markeredgecolor',[0.5 0.5 0.5])
        st=st+1;
    end
    hold on
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Cell diameter [\mum]')
ylabel('\sigma_{bb}(550) [m^2 particle^{-1}]')
xticks([1e0, 1e1, 1e2, 1e3])
text(0.05,0.95,'b','Units','normalized','fontsize',9,'fontweight','bold')
leg2=legend(leg,'Location','southeast','NumColumns',4,'Location','southoutside');
leg2.Position(1)=0.18;

set(gcf,'color','w')


%% this file generates new data that will be pasted in the optics file
% %(this is not the real optics file)
% 
% fid=fopen('paste_for_optics_file','w');
% 
% for i=1:12
%     fprintf(fid,'%s \n','   ');
%     idx=idx_start(i):idx_end(i);
%     for j=1:13
%     fprintf(fid,'%s \n',strcat([' ',num2str(opts_init(idx(j),1),'%.0f'),...
%         ' ',...
%         num2str(opts_init(idx(j),2),'%.4f'),...
%         '     ',...
%         num2str(opts_init(idx(j),3),'%.4f'),...
%         '     ',...
%         num2str(opts_init(idx(j),4),'%.5f'),...    
%         '    ',...
%         num2str(bbphy_mgC_type(j,i),'%.9f'),... 
%         '      ',...
%         num2str(opts_init(idx(j),6),'%.6f'),...
%         ]));
%     end
% %     st=st+1;
% end
% 
% fprintf(fid,'\n');
% 
% fclose(fid);

%% generate new optics file
% t = fileread('optics_plankton_new3.txt');
% expr = '[^\n]****[^\n]*';
% matches = regexp(t,expr,'match');
% 
% tt = importdata('optics_plankton_new3.txt');
% hd=tt.textdata;
% 
% fid=fopen('optics_file_gen.txt','w');
% 
% for i=1:6
% fprintf(fid,hd{i});
% fprintf(fid,'\n');
% end
% 
% % fclose(fid);
% 
% for i=1:12
% %     fprintf(fid,'%s \n','   ');
%     fprintf(fid,matches{i});
%     fprintf(fid,'\n');
%     idx=idx_start(i):idx_end(i);
%     for j=1:13
%         if j==1 
%             if i==5
%            fprintf(fid,'%s \n',strcat(['   ',num2str(opts_init(idx(j)-1,1),'%.0f'),...
%                '   ',num2str(opts_init(idx(j)-1,2),'%.2f'),...
%                '       ', num2str(opts_init(idx(j)-1,3),'%.2f'),...
%                '       ', num2str(opts_init(idx(j)-1,4),'%.2f'),...
%                '          ', num2str(opts_init(idx(j)-1,5),'%.2f'),...
%                '         ', num2str(opts_init(idx(j)-1,6),'%.2f'),...
%                ])) 
%             elseif i == 9
%            fprintf(fid,'%s \n',strcat(['   ',num2str(opts_init(idx(j)-1,1),'%.0f'),...
%                '  ',num2str(opts_init(idx(j)-1,2),'%.2f'),...
%                '      ', num2str(opts_init(idx(j)-1,3),'%.2f'),...
%                '       ', num2str(opts_init(idx(j)-1,4),'%.2f'),...
%                '          ', num2str(opts_init(idx(j)-1,5),'%.2f'),...
%                '         ', num2str(opts_init(idx(j)-1,6),'%.2f'),...
%                ]))                
%             else
%             fprintf(fid,'%s \n',strcat(['   ',num2str(opts_init(idx(j)-1,1),'%.0f'),...
%                '   ',num2str(opts_init(idx(j)-1,2),'%.2f'),...
%                '       ', num2str(opts_init(idx(j)-1,3),'%.2f'),...
%                '        ', num2str(opts_init(idx(j)-1,4),'%.2f'),...
%                '           ', num2str(opts_init(idx(j)-1,5),'%.2f'),...
%                '          ', num2str(opts_init(idx(j)-1,6),'%.2f'),...
%                ]))
%             end
%         end
%     fprintf(fid,'%s \n',strcat([' ',num2str(opts_init(idx(j),1),'%.0f'),...
%         ' ',...
%         num2str(opts_init(idx(j),2),'%.4f'),...
%         '     ',...
%         num2str(opts_init(idx(j),3),'%.4f'),...
%         '     ',...
%         num2str(opts_init(idx(j),4),'%.5f'),...    
%         '    ',...
%         num2str(bbphy_mgC_type(j,i),'%.9f'),... 
%         '      ',...
%         num2str(opts_init(idx(j),6),'%.6f'),...
%         ]));    
%     end
% %     fprintf(fid,'\n');
% end
% 
% fprintf(fid,'\n');
% 
% fclose(fid);