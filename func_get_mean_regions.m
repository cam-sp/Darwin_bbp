function means=func_get_mean_regions(bbp,biom,argo,ArgoDarwin,sat,ArgoSat,...
    idx_argo_biomes,idx_argo_basins,idx_sat_biomes,idx_sat_basins)

lambda_argo=700;
gamma=1;%0.8;%0.8;
bbp_argo_470=argo.bbp.*(470./lambda_argo).^(-gamma);
argo.bbp470=bbp_argo_470;

lambda_sat=443;
ArgoSat.bbp700=ArgoSat.bbp.*(700./lambda_sat).^(-gamma);

idx_finite_argo=isfinite(argo.bbp) & isfinite(argo.Chl);
idx_finite_argo_Sat=isfinite(ArgoSat.bbp) & isfinite(ArgoSat.Chl);
% idx_finite_sat=isfinite(chl_sat_mon) & isfinite(bbp_sat_mon);
% bbp_argo_binned2=bbp_argo_binned+2e-4;

bbp_argo_mean=NaN(1,168);
chl_argo_mean=NaN(1,168);
bbp_argo_darwin_mean=NaN(1,168);
chl_argo_darwin_mean=NaN(1,168);
biomphyto_argo_darwin_mean=NaN(1,168);
bbp470_argo_mean=NaN(1,168);
bbp_sat_mean=NaN(1,168);
chl_sat_mean=NaN(1,168);
bbp700_sat_mean=NaN(1,168);

bbp_argo_std=NaN(1,168);
chl_argo_std=NaN(1,168);
bbp_argo_darwin_std=NaN(1,168);
chl_argo_darwin_std=NaN(1,168);
biomphyto_argo_darwin_std=NaN(1,168);
bbp470_argo_std=NaN(1,168);
bbp_sat_std=NaN(1,168);
chl_sat_std=NaN(1,168);
bbp700_sat_std=NaN(1,168);

imon_argo_mean=NaN(1,168);
ilat_argo_mean=NaN(1,168);
ibiom_argo_mean=NaN(1,168);
ibas_argo_mean=NaN(1,168);
sst_argo_darwin_mean=NaN(1,168);
st=1;
for imon=1:12 %month
    idx_mon = argo.month==imon;
%     chl_sat_mon=sat.Chl(:,:,imon);
%     bbp_sat_mon=sat.bbp(:,:,imon);
%     bbp_sat_mon_700=bbp_sat_mon.*(700/470).^(-1);
    for ilat=1:2 %north/south
        if ilat==1
            idx_lat = argo.lat>0;
%             idx_lat_sat=sat.lat_2d>0;
        elseif ilat==2
            idx_lat = argo.lat<0;
%             idx_lat_sat=sat.lat_2d<0;
        end
        for ibiom=1:4 %biome
            if ibiom==1
                idx_biom = idx_argo_biomes.trop;
%                 idx_biom_sat=idx_sat_biomes.Trop;
            elseif ibiom==2
                idx_biom = idx_argo_biomes.oligo;
%                 idx_biom_sat=idx_sat_biomes.Oligo;
            elseif ibiom==3
                idx_biom = idx_argo_biomes.temp;
%                 idx_biom_sat=idx_sat_biomes.Temp;
            elseif ibiom==4    
                idx_biom = idx_argo_biomes.polar;
%                 idx_biom_sat=idx_sat_biomes.Polar;
            end
            
            if ibiom<4
                for ibas=1:2 %basin Atlantic or pacific
                    if ibas==1
                        idx_bas=idx_argo_basins.Atlantic;
%                         idx_bas_sat=logical(idx_sat_basins.Atlantic);
                    elseif ibas==2
                        idx_bas=idx_argo_basins.Pacific;
%                         idx_bas_sat=logical(idx_sat_basins.Pacific);
                    end
%                     idx1=idx_finite_sat & idx_bas_sat & idx_lat_sat & idx_biom_sat;
                    
                    idx = idx_mon & idx_lat & idx_biom & idx_bas & idx_finite_argo;
                    idx2 = idx_mon & idx_lat & idx_biom & idx_bas & idx_finite_argo & idx_finite_argo_Sat;
                    
                    bbp_argo_mean(st)=nanmean(argo.bbp(idx));
                    bbp470_argo_mean(st)=nanmean(bbp_argo_470(idx));
                    chl_argo_mean(st)=nanmean(argo.Chl(idx));
                    
                    bbp_argoSat_mean(st)=nanmean(argo.bbp(idx2));
                    bbp470_argoSat_mean(st)=nanmean(bbp_argo_470(idx2));
                    chl_argoSat_mean(st)=nanmean(argo.Chl(idx2));
                    
                    bbp_argo_darwin_mean(st)=nanmean(ArgoDarwin.bbp700(idx));
                    chl_argo_darwin_mean(st)=nanmean(ArgoDarwin.Chl(idx));
                    biomphyto_argo_darwin_mean(st)=nanmean(ArgoDarwin.biom_phyto(idx));
                    
                    bbp_sat_mean(st)=nanmean(ArgoSat.bbp(idx));
                    bbp700_sat_mean(st)=nanmean(ArgoSat.bbp700(idx));
                    chl_sat_mean(st)=nanmean(ArgoSat.Chl(idx));
                    
                    
                    bbp_argo_std(st)=nanstd(argo.bbp(idx));
                    bbp470_argo_std(st)=nanstd(bbp_argo_470(idx));
                    chl_argo_std(st)=nanstd(argo.Chl(idx));
                    
                    bbp_argoSat_std(st)=nanstd(argo.bbp(idx2));
                    bbp470_argoSat_std(st)=nanstd(bbp_argo_470(idx2));
                    chl_argoSat_std(st)=nanstd(argo.Chl(idx2));
                    
                    bbp_argo_darwin_std(st)=nanstd(ArgoDarwin.bbp700(idx));
                    chl_argo_darwin_std(st)=nanstd(ArgoDarwin.Chl(idx));
                    biomphyto_argo_darwin_std(st)=nanstd(ArgoDarwin.biom_phyto(idx));
                    
                    bbp_sat_std(st)= nanstd(ArgoSat.bbp(idx));
                    bbp700_sat_std(st)=nanstd(ArgoSat.bbp700(idx));
                    chl_sat_std(st)=nanstd(ArgoSat.Chl(idx));
                    
                    sst_argo_darwin_mean(st)=nanmean(ArgoDarwin.sst(idx));
                    imon_argo_mean(st)=imon;
                    ilat_argo_mean(st)=ilat;
                    ibiom_argo_mean(st)=ibiom;
                    ibas_argo_mean(st)=ibas;
                    st=st+1;
                end
            elseif ibiom==4
                    idx = idx_mon & idx_lat & idx_biom & idx_finite_argo;
                    idx2 = idx_mon & idx_lat & idx_biom & idx_finite_argo & idx_finite_argo_Sat;
%                     idx1= idx_isfinite_sat & idx_lat_sat & idx_biom_sat;
                    
                    bbp_argo_mean(st)=nanmean(argo.bbp(idx));
                    bbp470_argo_mean(st)=nanmean(bbp_argo_470(idx));
                    chl_argo_mean(st)=nanmean(argo.Chl(idx));
                    
                    bbp_argoSat_mean(st)=nanmean(argo.bbp(idx2));
                    bbp470_argoSat_mean(st)=nanmean(bbp_argo_470(idx2));
                    chl_argoSat_mean(st)=nanmean(argo.Chl(idx2));
                    
                    bbp_argo_darwin_mean(st)=nanmean(ArgoDarwin.bbp700(idx));
                    chl_argo_darwin_mean(st)=nanmean(ArgoDarwin.Chl(idx));
                    biomphyto_argo_darwin_mean(st)=nanmean(ArgoDarwin.biom_phyto(idx));
                    
                    bbp_sat_mean(st)= nanmean(ArgoSat.bbp(idx));
                    bbp700_sat_mean(st)=nanmean(ArgoSat.bbp700(idx));
                    chl_sat_mean(st)=nanmean(ArgoSat.Chl(idx));   
                    
                    sst_argo_darwin_mean(st)=nanmean(ArgoDarwin.sst(idx));
                    
                    bbp_argo_std(st)=nanstd(argo.bbp(idx));
                    bbp470_argo_std(st)=nanstd(bbp_argo_470(idx));
                    chl_argo_std(st)=nanstd(argo.Chl(idx));
                    
                    bbp_argoSat_std(st)=nanstd(argo.bbp(idx2));
                    bbp470_argoSat_std(st)=nanstd(bbp_argo_470(idx2));
                    chl_argoSat_std(st)=nanstd(argo.Chl(idx2));
                    
                    bbp_argo_darwin_std(st)=nanstd(ArgoDarwin.bbp700(idx));
                    chl_argo_darwin_std(st)=nanstd(ArgoDarwin.Chl(idx));
                    biomphyto_argo_darwin_std(st)=nanstd(ArgoDarwin.biom_phyto(idx));
                    
                    bbp_sat_std(st)= nanstd(ArgoSat.bbp(idx));
                    bbp700_sat_std(st)=nanstd(ArgoSat.bbp700(idx));
                    chl_sat_std(st)=nanstd(ArgoSat.Chl(idx));   
                    
                    sst_argo_darwin_mean(st)=nanstd(ArgoDarwin.sst(idx));                   
                    imon_argo_mean(st)=imon;
                    ilat_argo_mean(st)=ilat;
                    ibiom_argo_mean(st)=ibiom;
                    ibas_argo_mean(st)=0;
                    st=st+1;                
            end
        end
    end
end

bbp_argo_mean(bbp_argo_mean==0)=NaN;
chl_argo_mean(chl_argo_mean==0)=NaN;





                                  
means.argo.bbp=bbp_argo_mean;
means.argo.bb470=bbp470_argo_mean;
means.argo.Chl=chl_argo_mean;
means.argo.bbp_std=bbp_argo_std;
means.argo.bb470_std=bbp470_argo_std;
means.argo.Chl_std=chl_argo_std;

means.argoSat.bbp=bbp_argoSat_mean;
means.argoSat.bb470=bbp470_argoSat_mean;
means.argoSat.Chl=chl_argoSat_mean;
means.argoSat.bbp_std=bbp_argoSat_std;
means.argoSat.bb470_std=bbp470_argoSat_std;
means.argoSat.Chl_std=chl_argoSat_std;

means.darwin.bbp=bbp_argo_darwin_mean;
means.darwin.Chl=chl_argo_darwin_mean;
means.darwin.biom_phyto=biomphyto_argo_darwin_mean;
means.darwin.bbp_std=bbp_argo_darwin_std;
means.darwin.Chl_std=chl_argo_darwin_std;
means.darwin.biom_phyto_std=biomphyto_argo_darwin_std;

means.sat.bbp=bbp_sat_mean;
means.sat.bbp700=bbp700_sat_mean;
means.sat.Chl=chl_sat_mean;
means.sat.bbp_std=bbp_sat_std;
means.sat.bbp700_std=bbp700_sat_std;
means.sat.Chl_std=chl_sat_std;

means.idx.imon=imon_argo_mean;
means.idx.ilat=ilat_argo_mean;
means.idx.ibiom=ibiom_argo_mean;
means.idx.ibas=ibas_argo_mean;
end