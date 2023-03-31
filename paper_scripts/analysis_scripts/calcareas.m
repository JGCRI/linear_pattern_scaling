clear all;

residuals_dir='/Users/bkravitz/Documents/github repositories/cmip6_patterns/outputs/resids_tgavs';
cd(residuals_dir);

models={'HadGEM3-GC31-MM','GFDL-ESM4','GFDL-CM4','IPSL-CM6A-LR','CNRM-CM6-1',...
 'GISS-E2-1-G','GISS-E2-1-H','BCC-CSM2-MR','BCC-ESM1','CNRM-ESM2-1',...
 'MIROC6','AWI-CM-1-1-MR','EC-Earth3-LR','MRI-ESM2-0','CESM2-WACCM',...
 'CESM2','SAM0-UNICON','GISS-E2-1-G-CC','UKESM1-0-LL','EC-Earth3',...
 'CanESM5','CanESM5-CanOE','EC-Earth3-Veg','HadGEM3-GC31-LL',...
 'MPI-ESM-1-2-HAM','NESM3','CAMS-CSM1-0','MPI-ESM1-2-LR','MPI-ESM1-2-HR',...
 'MCM-UA-1-0','NorESM2-LM','FGOALS-g3','FGOALS-f3-L','MIROC-ES2L',...
 'FIO-ESM-2-0','NorCPM1','NorESM1-F','CNRM-CM6-1-HR','ACCESS-CM2',...
 'NorESM2-MM','ACCESS-ESM1-5','IITM-ESM','CESM2-FV2','CESM2-WACCM-FV2',...
 'GISS-E2-2-G','GISS-E2-2-H','TaiESM1','AWI-ESM-1-1-LR','CIESM',...
 'CMCC-CM2-SR5','EC-Earth3-AerChem','IPSL-CM5A2-INCA','CMCC-CM2-HR4',...
 'EC-Earth3-Veg-LR','CAS-ESM2-0','EC-Earth3-CC','CMCC-ESM2','MIROC-ES2H',...
 'ICON-ESM-LR','IPSL-CM6A-LR-INCA'};

scenarios={'historical','ssp126','ssp245','ssp370','ssp434','ssp585'};

lats=-90:1:90;
lons=0.5:1:359.5;
latedges=[-90 -89.5:1:89.5 90];
r=6371000;
colatedges=pi/2-latedges*pi/180;
dphi=abs(pi/180);
aw=repmat(abs(diff(cos(colatedges)))*dphi.*(r.^2),[360 1]);

ncid=netcdf.open('sftlf_CCSM4.nc','NC_NOWRITE');
lfid=netcdf.inqVarID(ncid,'sftlf');
lfC=netcdf.getVar(ncid,lfid);
latid=netcdf.inqVarID(ncid,'lat');
lonid=netcdf.inqVarID(ncid,'lon');
lat=netcdf.getVar(ncid,latid);
lon=netcdf.getVar(ncid,lonid);
netcdf.close(ncid);
lfC=lfC/100;
xi2=repmat(lons,[length(lats) 1])';
yi2=repmat(lats,[length(lons) 1]);
x2=repmat(lon',[length(lat) 1])';
y2=repmat(lat',[length(lon) 1]);
lf=nanfill2d(interp2(y2,x2,squeeze(lfC),yi2,xi2,'*linear'));

load('procresids_tas.mat');

errormap_tas=zeros(360,181,6);
divs=zeros(360,181,6);
for k1=1:360;
for k2=1:181;
for k3=1:length(models);
for k4=1:length(scenarios);
    temp=squeeze(resids(k1,k2,k3,k4));
    whichstd=squeeze(std20yr(k1,k2,k3,k4));
    if isnan(temp)==0;
        divs(k1,k2,k4)=divs(k1,k2,k4)+1;
    end
    if abs(temp)<whichstd;
        errormap_tas(k1,k2,k4)=errormap_tas(k1,k2,k4)+1;
    end
end
end
end
end
errormap_tas=errormap_tas./divs*100;

erroraggregate_tas=zeros(360,181);
for k1=1:360;
for k2=1:181;
    K=find(errormap_tas(k1,k2,:)>90);
    erroraggregate_tas(k1,k2)=length(K);
end
end

tas_area = sum(sum(floor(erroraggregate_tas/6).*aw))./sum(sum(aw));
tas_area_land = sum(sum(floor(erroraggregate_tas/6).*aw.*lf))./sum(sum(aw.*lf));

load('procresids_pr.mat');

errormap_pr=zeros(360,181,6);
divs=zeros(360,181,6);
for k1=1:360;
for k2=1:181;
for k3=1:length(models);
for k4=1:length(scenarios);
    temp=squeeze(resids(k1,k2,k3,k4));
    whichstd=squeeze(std20yr(k1,k2,k3,k4));
    if isnan(temp)==0;
        divs(k1,k2,k4)=divs(k1,k2,k4)+1;
    end
    if abs(temp)<whichstd;
        errormap_pr(k1,k2,k4)=errormap_pr(k1,k2,k4)+1;
    end
end
end
end
end
errormap_pr=errormap_pr./divs*100;

erroraggregate_pr=zeros(360,181);
for k1=1:360;
for k2=1:181;
    K=find(errormap_pr(k1,k2,:)>90);
    erroraggregate_pr(k1,k2)=length(K);
end
end

pr_area = sum(sum(floor(erroraggregate_pr/6).*aw))./sum(sum(aw));
pr_area_land = sum(sum(floor(erroraggregate_pr/6).*aw.*lf))./sum(sum(aw.*lf));

load('procresids_hurs.mat');

errormap_hurs=zeros(360,181,6);
divs=zeros(360,181,6);
for k1=1:360;
for k2=1:181;
for k3=1:length(models);
for k4=1:length(scenarios);
    temp=squeeze(resids(k1,k2,k3,k4));
    whichstd=squeeze(std20yr(k1,k2,k3,k4));
    if isnan(temp)==0;
        divs(k1,k2,k4)=divs(k1,k2,k4)+1;
    end
    if abs(temp)<whichstd;
        errormap_hurs(k1,k2,k4)=errormap_hurs(k1,k2,k4)+1;
    end
end
end
end
end
errormap_hurs=errormap_hurs./divs*100;

erroraggregate_hurs=zeros(360,181);
for k1=1:360;
for k2=1:181;
    K=find(errormap_hurs(k1,k2,:)>90);
    erroraggregate_hurs(k1,k2)=length(K);
end
end

hurs_area = sum(sum(floor(erroraggregate_hurs/6).*aw))./sum(sum(aw));
hurs_area_land = sum(sum(floor(erroraggregate_hurs/6).*aw.*lf))./sum(sum(aw.*lf));

