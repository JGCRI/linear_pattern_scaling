% gleckler.m
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com) and Abigail Snyder
% Last updated 20 September 2022
%
% This script will produce a Gleckler plot (an exceptionally informative heat
% map) that describes the residuals for each variable and region.  It will also
% make some box plots describing model spread.  A lot of this is hard-coded.

% Preamble begins here
clear all;

% Directories where the various files are.  Feel free to change this.
residuals_dir='/Users/bkravitz/Documents/github repositories/cmip6_patterns/outputs/resids_tgavs';
patterns_dir='/Users/bkravitz/Documents/github repositories/cmip6_patterns/outputs';

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

% Loading in the already-processed residuals.
cd(residuals_dir);
load('procresids_tas.mat');
tasresids=resids./std20yr;

load('procresids_pr.mat');
prresids=resids./std20yr;

load('procresids_hurs.mat');
hursresids=resids./std20yr;

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

% Setting up the matrices that we are going to plot
gleckler_plot=zeros(18,length(models),length(scenarios));
tw=repmat(weights,[1 1 1 length(models) length(scenarios)]);
temp=permute(repmat(tasresids,[1 1 1 1 6]),[1 2 5 3 4]);
gleckler_plot(1:6,:,:)=sqrt(squeeze(sum(sum(temp.^2.*tw,2),1))./squeeze(sum(sum(tw,2),1)));
temp=permute(repmat(prresids,[1 1 1 1 6]),[1 2 5 3 4]);
gleckler_plot(7:12,:,:)=sqrt(squeeze(sum(sum(temp.^2.*tw,2),1))./squeeze(sum(sum(tw,2),1)));
temp=permute(repmat(hursresids,[1 1 1 1 6]),[1 2 5 3 4]);
gleckler_plot(13:18,:,:)=sqrt(squeeze(sum(sum(temp.^2.*tw,2),1))./squeeze(sum(sum(tw,2),1)));

thelabs={'tas global','tas land','tas ocean','tas tropics','tas midlats','tas hilats','pr global','pr land','pr ocean','pr tropics','pr midlats','pr hilats','hurs global','hurs land','hurs ocean','hurs tropics','hurs midlats','hurs hilats'};

% Making the heat map 
figure(1);
for k=1:6;
    subplot(3,2,k);
    heatmap(models,thelabs,squeeze(gleckler_plot(:,:,k)),'ColorLimits',[0 2],'ColorMap',flipud(autumn(33)),'MissingDataColor','w');
    title(scenarios{k});
end

% Making the box plots
figure(2);
for k=1:6;
    subplot(3,2,k);
    temp=squeeze(gleckler_plot(:,:,k));
    temp(find(isnan(temp)==1))=0;
    boxplot(transpose(temp),'plotstyle','compact','boxstyle','outline','medianstyle','line','datalim',[-Inf Inf],'whisker',30,'widths',0.25,'labels',thelabs);
    title(scenarios{k},'FontSize',15);
    ylabel('RMSE (# of standard devs)','FontSize',15);
    set(gca,'FontSize',15);
    ylim([0 3]);
end
