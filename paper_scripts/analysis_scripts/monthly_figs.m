% monthly_figs.m
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com) and Abigail Snyder
% Last updated 20 September 2022
%
% This script will create the monthly figures (8-10 and most of the supplemental figures).
%
% A lot of the code repeats because we hard-coded each variable.  It would
% not be too difficult to make this more flexible for multiple variables,
% but we have not done that yet.

% Preamble begins here
clear all;

% Directories where the various files are.  Feel free to change this.
residuals_dir='/Users/bkravitz/Documents/github repositories/cmip6_patterns/outputs/resids_tgavs';
cd(residuals_dir);

% Models to look through - you can edit this as well.
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

% Scenarios - you can even edit this one if you want to.
scenarios={'historical','ssp126','ssp245','ssp370','ssp434','ssp585'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN PART OF THE SCRIPT BEGINS HERE %%%%%
%%%%%%%%%%%% EDIT AT YOUR OWN RISK %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mos={'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};

% This block of code is for computing an area weighting.
% Some of the figures below divide things up into land vs ocean
% or high latitudes vs low latitudes and then does an average.
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
r1=intersect(find(lats>-23.5),find(lats<23.5));
r2=union(intersect(find(lats>23.5),find(lats<66.5)),intersect(find(lats<-23.5),find(lats>-66.5)));
r3=union(find(lats>66.5),find(lats<-66.5));
weights=repmat(aw,[1 1 6]);
weights(:,:,2)=weights(:,:,2).*lf;
weights(:,:,3)=weights(:,:,3).*(1-lf);
weights(:,setdiff(1:181,r1),4)=0;
weights(:,setdiff(1:181,r2),5)=0;
weights(:,setdiff(1:181,r3),6)=0;

% Creating maps of where the residuals are less than one standard deviation of the ESM data
load('procresids_monthly_tas.mat');
errormap_tas=zeros(360,181,12,6);
divs=zeros(360,181,12,6);
for k1=1:360;
for k2=1:181;
for k3=1:length(models);
for k4=1:length(scenarios);
for k5=1:12;
    temp=squeeze(resids(k1,k2,k5,k3,k4));
    whichstd=squeeze(std20yr(k1,k2,k5,k3,k4));
    if isnan(temp)==0;
        divs(k1,k2,k5,k4)=divs(k1,k2,k5,k4)+1;
    end
    if abs(temp)<whichstd;
        errormap_tas(k1,k2,k5,k4)=errormap_tas(k1,k2,k5,k4)+1;
    end
end
end
end
end
end
errormap_tas=errormap_tas./divs*100;

for k1=1:12;
for k2=1:length(scenarios);
    temp=squeeze(errormap_tas(:,:,k1,k2));
    fname=['errormap_tas_' mos{k1} '_' scenarios{k2} '.bin'];
    writebinary(temp,fname,-1e30);
end
end

% Aggregating the error maps to see if they are scenario-dependent
erroraggregate_tas=zeros(360,181,12);
for k1=1:360;
for k2=1:181;
for k3=1:12;
    K=find(errormap_tas(k1,k2,k3,:)>90);
    erroraggregate_tas(k1,k2,k3)=length(K);
end
end
end

for k1=1:12;
    temp=squeeze(erroraggregate_tas(:,:,k1));
    fname=['erroraggregate_tas_' mos{k1} '.bin'];
    writebinary(temp,fname,-1e30);
end

% Box plot to look at climatologies of the monthly pattern scaling residuals
thelabs={'global','land','ocean','tropics','midlats','hilats'};
figure(14);
for k=1:6;
    subplot(3,2,k);
    hold on;
    w=repmat(squeeze(weights(:,:,k)),[1 1 12 length(models) length(scenarios)]);
    y=squeeze(sum(sum(resids.*w,2),1))./squeeze(sum(sum(w,2),1));
    y=reshape(y,[12 length(models)*length(scenarios)]);
    y=cat(1,y,y(1,:));
    plot(1:13,y);
    plot(1:13,nanmean(y,2),'k','LineWidth',2);
    title(thelabs{k});
    set(gca,'XTick',1:13,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J'},'FontSize',15);
    xlim([1 13]);
    ylim([-0.5 0.5]);
end


load('procresids_monthly_pr.mat');
errormap_pr=zeros(360,181,12,6);
divs=zeros(360,181,12,6);
for k1=1:360;
for k2=1:181;
for k3=1:length(models);
for k4=1:length(scenarios);
for k5=1:12;
    temp=squeeze(resids(k1,k2,k5,k3,k4));
    whichstd=squeeze(std20yr(k1,k2,k5,k3,k4));
    if isnan(temp)==0;
        divs(k1,k2,k5,k4)=divs(k1,k2,k5,k4)+1;
    end
    if abs(temp)<whichstd;
        errormap_pr(k1,k2,k5,k4)=errormap_pr(k1,k2,k5,k4)+1;
    end
end
end
end
end
end
errormap_pr=errormap_pr./divs*100;

for k1=1:12;
for k2=1:length(scenarios);
    temp=squeeze(errormap_pr(:,:,k1,k2));
    fname=['errormap_pr_' mos{k1} '_' scenarios{k2} '.bin'];
    writebinary(temp,fname,-1e30);
end
end

erroraggregate_pr=zeros(360,181,12);
for k1=1:360;
for k2=1:181;
for k3=1:12;
    K=find(errormap_pr(k1,k2,k3,:)>90);
    erroraggregate_pr(k1,k2,k3)=length(K);
end
end
end


for k1=1:12;
    temp=squeeze(erroraggregate_pr(:,:,k1));
    fname=['erroraggregate_pr_' mos{k1} '.bin'];
    writebinary(temp,fname,-1e30);
end

thelabs={'global','land','ocean','tropics','midlats','hilats'};
figure(114);
for k=1:6;
    subplot(3,2,k);
    hold on;
    w=repmat(squeeze(weights(:,:,k)),[1 1 12 length(models) length(scenarios)]);
    y=squeeze(sum(sum(86400*resids.*w,2),1))./squeeze(sum(sum(w,2),1));
    y=reshape(y,[12 length(models)*length(scenarios)]);
    y=cat(1,y,y(1,:));
    plot(1:13,y);
    plot(1:13,nanmean(y,2),'k','LineWidth',2);
    title(thelabs{k});
    set(gca,'XTick',1:13,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J'},'FontSize',15);
    xlim([1 13]);
    ylim([-0.15	0.15]);
end

load('procresids_monthly_hurs.mat');
errormap_hurs=zeros(360,181,12,6);
divs=zeros(360,181,12,6);
for k1=1:360;
for k2=1:181;
for k3=1:length(models);
for k4=1:length(scenarios);
for k5=1:12;
    temp=squeeze(resids(k1,k2,k5,k3,k4));
    whichstd=squeeze(std20yr(k1,k2,k5,k3,k4));
    if isnan(temp)==0;
        divs(k1,k2,k5,k4)=divs(k1,k2,k5,k4)+1;
    end
    if abs(temp)<whichstd;
        errormap_hurs(k1,k2,k5,k4)=errormap_hurs(k1,k2,k5,k4)+1;
    end
end
end
end
end
end
errormap_hurs=errormap_hurs./divs*100;

for k1=1:12;
for k2=1:length(scenarios);
    temp=squeeze(errormap_hurs(:,:,k1,k2));
    fname=['errormap_hurs_' mos{k1} '_' scenarios{k2} '.bin'];
    writebinary(temp,fname,-1e30);
end
end

erroraggregate_hurs=zeros(360,181,12);
for k1=1:360;
for k2=1:181;
for k3=1:12;
    K=find(errormap_hurs(k1,k2,k3,:)>90);
    erroraggregate_hurs(k1,k2,k3)=length(K);
end
end
end


for k1=1:12;
    temp=squeeze(erroraggregate_hurs(:,:,k1));
    fname=['erroraggregate_hurs_' mos{k1} '.bin'];
    writebinary(temp,fname,-1e30);
end

thelabs={'global','land','ocean','tropics','midlats','hilats'};
figure(214);
for k=1:6;
    subplot(3,2,k);
    hold on;
    w=repmat(squeeze(weights(:,:,k)),[1 1 12 length(models) length(scenarios)]);
    y=squeeze(sum(sum(resids.*w,2),1))./squeeze(sum(sum(w,2),1));
    y=reshape(y,[12 length(models)*length(scenarios)]);
    y=cat(1,y,y(1,:));
    plot(1:13,y);
    plot(1:13,nanmean(y,2),'k','LineWidth',2);
    title(thelabs{k});
    set(gca,'XTick',1:13,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D','J'},'FontSize',15);
    xlim([1 13]);
    ylim([-2 2]);
end
