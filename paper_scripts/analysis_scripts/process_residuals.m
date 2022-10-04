% process_residuals.m
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com) and Abigail Snyder
% Last updated 20 September 2022
%
% This script will process the residual files from annual mean pattern scaling.
% It will produce matlab files as output.  There is some code to make plots,
% which we have commented out.
%
% A lot of the code repeats because we hard-coded each variable.  It would
% not be too difficult to make this more flexible for multiple variables,
% but we have not done that yet.

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

% Storage tanks for the residuals and standard deviations
resids=zeros(360,181,length(models),length(scenarios)); % = actual ESM minus linear fit
std20yr=zeros(360,181,length(models),length(scenarios));

% Processing the temperature residuals.
% This block of code reads in the residual and the patterns (slope and intercept).
% It does both because we compute the standard deviation from the original model
% output, which we can reconstruct from the residual and the patterns.
for k1=1:length(models);
for k2=1:length(scenarios);
    fname=[models{k1} '_' scenarios{k2} '_tas_annual_pattern_resids.nc'];
    try % it loops through all models, but sometimes the model/scenario combo does not exist
        cd(residuals_dir);
        ncid=netcdf.open(fname,'NC_NOWRITE');
        latid=netcdf.inqVarID(ncid,'lat');
        lonid=netcdf.inqVarID(ncid,'lon');
        varid=netcdf.inqVarID(ncid,'tas');
        lat=netcdf.getVar(ncid,latid);
        lon=netcdf.getVar(ncid,lonid);
        invar=netcdf.getVar(ncid,varid);
        netcdf.close(ncid);
        lats=-90:1:90;
        lons=0.5:1:359.5;
        xi2=repmat(lons,[length(lats) 1])';
        yi2=repmat(lats,[length(lons) 1]);
        x2=repmat(lon',[length(lat) 1])';
        y2=repmat(lat',[length(lon) 1]);
        temp=zeros(360,181,20);
        c=1;
        n=size(invar);
        for k=n(3)-19:n(3);
            temp(:,:,c)=nanfill2d(interp2(y2,x2,squeeze(invar(:,:,k)),yi2,xi2,'*linear'));
            c=c+1;
        end
        resids(:,:,k1,k2)=squeeze(mean(temp,3));
        
        cd(patterns_dir);
        fname=[models{k1} '_' scenarios{k2} '_tas_annual_pattern.nc'];
        ncid=netcdf.open(fname,'NC_NOWRITE');
        varid=netcdf.inqVarID(ncid,'slope');
        slope=netcdf.getVar(ncid,varid);
        varid2=netcdf.inqVarID(ncid,'intercept');
        intercept=netcdf.getVar(ncid,varid2);
        netcdf.close(ncid);
        
        cd(residuals_dir);
        fname=[models{k1} '_' scenarios{k2} '_ensemble_avg_tgav.nc'];
        ncid=netcdf.open(fname,'NC_NOWRITE');
        varid=netcdf.inqVarID(ncid,'tas');
        tgav=netcdf.getVar(ncid,varid);
        netcdf.close(ncid);
        
        reconstructed_last20=zeros(length(lon),length(lat),20);
        I=length(tgav);
        c=1;
        for j3=I-19:I;
            reconstructed_last20(:,:,c)=slope*tgav(j3)+intercept;
            c=c+1;
        end
        reconstructed_last20=reconstructed_last20+invar(:,:,n(3)-19:n(3));
        std20yr_noregrid=std(reconstructed_last20,0,3);
        std20yr(:,:,k1,k2)=nanfill2d(interp2(y2,x2,std20yr_noregrid,yi2,xi2,'*linear'));       
    catch
        resids(:,:,k1,k2)=NaN;
        std20yr(:,:,k1,k2)=NaN;
    end
end
end
resids=-resids; % Residuals are output as ESM minus pattern scaling, and Ben thinks the opposite way

% saving the output
save('procresids_tas.mat','resids','std20yr');

% Creating maps of where the residuals are less than one standard deviation of the ESM data
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

%figure(1);
%suptitle('Pct of models where abs(tas)<1 std');
%hold on;
for k=1:6;
%    subplot(3,2,k);
%    contourf(transpose(squeeze(errormap_tas(:,:,k))));
%    colorbar;
%    caxis([0 100]);
%    title(scenarios{k});
    writebinary(squeeze(errormap_tas(:,:,k)),['annualresidual_tas_' scenarios{k} '.bin'],-1e30);
end

% Aggregating the error maps to see if they are scenario-dependent
erroraggregate_tas=zeros(360,181);
for k1=1:360;
for k2=1:181;
    K=find(errormap_tas(k1,k2,:)>90);
    erroraggregate_tas(k1,k2)=length(K);
end
end

writebinary(erroraggregate_tas,'annualaggregate_tas.bin',-1e30);

%figure(2);
%contourf(transpose(erroraggregate_tas));
%colorbar;

% Calculating mean, minimum, and maximum (across all models) residual at each grid point
% We did not use this, so we have commented it out.
%minresid=squeeze(min(resids,[],3));
%maxresid=squeeze(max(resids,[],3));
%meanresid=squeeze(nanmean(resids,3));
%r_over_s=(maxresid-minresid)./squeeze(std(std20yr,0,3,'omitnan'));
%figure(12);
%hold on;
%for k=1:6;
%    subplot(6,3,(k-1)*3+1);
%    contourf(transpose(squeeze(meanresid(:,:,k))));
%    colorbar;
%    title('mean');
%    ylabel(scenarios{k});
%    subplot(6,3,(k-1)*3+2);
%    contourf(transpose(squeeze(maxresid(:,:,k))));
%    colorbar;
%    title('max across models');
%    ylabel(scenarios{k});
%    subplot(6,3,(k-1)*3+3);
%    contourf(transpose(squeeze(minresid(:,:,k))));
%    colorbar;
%    title('min across models');
%    ylabel(scenarios{k});
%    writebinary(squeeze(meanresid(:,:,k)),['annualmean_tas_' scenarios{k} '.bin'],-1e30);
%    writebinary(squeeze(minresid(:,:,k)),['annualmin_tas_' scenarios{k} '.bin'],-1e30);
%    writebinary(squeeze(maxresid(:,:,k)),['annualmax_tas_' scenarios{k} '.bin'],-1e30);
%end
%figure(112);
%for k=1:6;
%    subplot(3,2,k);
%    contourf(transpose(squeeze(r_over_s(:,:,k))));
%    colorbar;
%    title(scenarios{k});
%    writebinary(squeeze(r_over_s(:,:,k)),['annualrs_tas_' scenarios{k} '.bin'],-1e30);    
%end

% Repeating the above code but for precipitation
resids=zeros(360,181,length(models),length(scenarios)); % = actual ESM minus linear fit
std20yr=zeros(360,181,length(models),length(scenarios));

for k1=1:length(models);
for k2=1:length(scenarios);
    fname=[models{k1} '_' scenarios{k2} '_pr_annual_pattern_resids.nc'];
    try
        cd(residuals_dir);
        ncid=netcdf.open(fname,'NC_NOWRITE');
        latid=netcdf.inqVarID(ncid,'lat');
        lonid=netcdf.inqVarID(ncid,'lon');
        varid=netcdf.inqVarID(ncid,'pr');
        lat=netcdf.getVar(ncid,latid);
        lon=netcdf.getVar(ncid,lonid);
        invar=netcdf.getVar(ncid,varid);
        netcdf.close(ncid);
        lats=-90:1:90;
        lons=0.5:1:359.5;
        xi2=repmat(lons,[length(lats) 1])';
        yi2=repmat(lats,[length(lons) 1]);
        x2=repmat(lon',[length(lat) 1])';
        y2=repmat(lat',[length(lon) 1]);
        temp=zeros(360,181,20);
        c=1;
        n=size(invar);
        for k=n(3)-19:n(3);
            temp(:,:,c)=nanfill2d(interp2(y2,x2,squeeze(invar(:,:,k)),yi2,xi2,'*linear'));
            c=c+1;
        end
        resids(:,:,k1,k2)=squeeze(mean(temp,3));
        
        cd(patterns_dir);
        fname=[models{k1} '_' scenarios{k2} '_pr_annual_pattern.nc'];
        ncid=netcdf.open(fname,'NC_NOWRITE');
        varid=netcdf.inqVarID(ncid,'slope');
        slope=netcdf.getVar(ncid,varid);
        varid2=netcdf.inqVarID(ncid,'intercept');
        intercept=netcdf.getVar(ncid,varid2);
        netcdf.close(ncid);
        
        cd(residuals_dir);
        fname=[models{k1} '_' scenarios{k2} '_ensemble_avg_tgav.nc'];
        ncid=netcdf.open(fname,'NC_NOWRITE');
        varid=netcdf.inqVarID(ncid,'tas');
        tgav=netcdf.getVar(ncid,varid);
        netcdf.close(ncid);
        
        reconstructed_last20=zeros(length(lon),length(lat),20);
        I=length(tgav);
        c=1;
        for j3=I-19:I;
            reconstructed_last20(:,:,c)=slope*tgav(j3)+intercept;
            c=c+1;
        end
        reconstructed_last20=reconstructed_last20+invar(:,:,n(3)-19:n(3));
        std20yr_noregrid=std(reconstructed_last20,0,3);
        std20yr(:,:,k1,k2)=nanfill2d(interp2(y2,x2,std20yr_noregrid,yi2,xi2,'*linear'));       
    catch
        resids(:,:,k1,k2)=NaN;
        std20yr(:,:,k1,k2)=NaN;
    end
end
end
resids=-resids;

save('procresids_pr.mat','resids','std20yr');

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

%figure(3);
%suptitle('Pct of models where abs(pr)<1 std');
%hold on;
for k=1:6;
%    subplot(3,2,k);
%    contourf(transpose(squeeze(errormap_pr(:,:,k))));
%    colorbar;
%    caxis([0 100]);
%    title(scenarios{k});
    writebinary(squeeze(errormap_pr(:,:,k)),['annualresidual_pr_' scenarios{k} '.bin'],-1e30);
end

erroraggregate_pr=zeros(360,181);
for k1=1:360;
for k2=1:181;
    K=find(errormap_pr(k1,k2,:)>90);
    erroraggregate_pr(k1,k2)=length(K);
end
end

writebinary(erroraggregate_pr,'annualaggregate_pr.bin',-1e30);

%figure(4);
%contourf(transpose(erroraggregate_pr));
%colorbar;

%minresid=squeeze(min(resids,[],3));
%maxresid=squeeze(max(resids,[],3));
%meanresid=squeeze(nanmean(resids,3));
%r_over_s=(maxresid-minresid)./squeeze(std(std20yr,0,3,'omitnan'));
%figure(14);
%hold on;
%for k=1:6;
%    subplot(6,3,(k-1)*3+1);
%    contourf(transpose(squeeze(meanresid(:,:,k))));
%    colorbar;
%    title('mean');
%    ylabel(scenarios{k});
%    subplot(6,3,(k-1)*3+2);
%    contourf(transpose(squeeze(maxresid(:,:,k))));
%    colorbar;
%    title('max across models');
%    ylabel(scenarios{k});
%    subplot(6,3,(k-1)*3+3);
%    contourf(transpose(squeeze(minresid(:,:,k))));
%    colorbar;
%    title('min across models');
%    ylabel(scenarios{k});
%    writebinary(squeeze(meanresid(:,:,k)),['annualmean_pr_' scenarios{k} '.bin'],-1e30);
%    writebinary(squeeze(minresid(:,:,k)),['annualmin_pr_' scenarios{k} '.bin'],-1e30);
%    writebinary(squeeze(maxresid(:,:,k)),['annualmax_pr_' scenarios{k} '.bin'],-1e30);
%end
%figure(114);
%for k=1:6;
%    subplot(3,2,k);
%    contourf(transpose(squeeze(r_over_s(:,:,k))));
%    colorbar;
%    title(scenarios{k});
%    writebinary(squeeze(r_over_s(:,:,k)),['annualrs_pr_' scenarios{k} '.bin'],-1e30);    
%end

% Repeating the above code but for relative humidity
resids=zeros(360,181,length(models),length(scenarios)); % = actual ESM minus linear fit
std20yr=zeros(360,181,length(models),length(scenarios));

for k1=1:length(models);
for k2=1:length(scenarios);
    fname=[models{k1} '_' scenarios{k2} '_hurs_annual_pattern_resids.nc'];
    try
        cd(residuals_dir);
        ncid=netcdf.open(fname,'NC_NOWRITE');
        latid=netcdf.inqVarID(ncid,'lat');
        lonid=netcdf.inqVarID(ncid,'lon');
        varid=netcdf.inqVarID(ncid,'hurs');
        lat=netcdf.getVar(ncid,latid);
        lon=netcdf.getVar(ncid,lonid);
        invar=netcdf.getVar(ncid,varid);
        netcdf.close(ncid);
        lats=-90:1:90;
        lons=0.5:1:359.5;
        xi2=repmat(lons,[length(lats) 1])';
        yi2=repmat(lats,[length(lons) 1]);
        x2=repmat(lon',[length(lat) 1])';
        y2=repmat(lat',[length(lon) 1]);
        temp=zeros(360,181,20);
        c=1;
        n=size(invar);
        for k=n(3)-19:n(3);
            temp(:,:,c)=nanfill2d(interp2(y2,x2,squeeze(invar(:,:,k)),yi2,xi2,'*linear'));
            c=c+1;
        end
        resids(:,:,k1,k2)=squeeze(mean(temp,3));
        
        cd(patterns_dir);
        fname=[models{k1} '_' scenarios{k2} '_hurs_annual_pattern.nc'];
        ncid=netcdf.open(fname,'NC_NOWRITE');
        varid=netcdf.inqVarID(ncid,'slope');
        slope=netcdf.getVar(ncid,varid);
        varid2=netcdf.inqVarID(ncid,'intercept');
        intercept=netcdf.getVar(ncid,varid2);
        netcdf.close(ncid);
        
        cd(residuals_dir);
        fname=[models{k1} '_' scenarios{k2} '_ensemble_avg_tgav.nc'];
        ncid=netcdf.open(fname,'NC_NOWRITE');
        varid=netcdf.inqVarID(ncid,'tas');
        tgav=netcdf.getVar(ncid,varid);
        netcdf.close(ncid);
        
        reconstructed_last20=zeros(length(lon),length(lat),20);
        I=length(tgav);
        c=1;
        for j3=I-19:I;
            reconstructed_last20(:,:,c)=slope*tgav(j3)+intercept;
            c=c+1;
        end
        reconstructed_last20=reconstructed_last20+invar(:,:,n(3)-19:n(3));
        std20yr_noregrid=std(reconstructed_last20,0,3);
        std20yr(:,:,k1,k2)=nanfill2d(interp2(y2,x2,std20yr_noregrid,yi2,xi2,'*linear'));       
    catch
        resids(:,:,k1,k2)=NaN;
        std20yr(:,:,k1,k2)=NaN;
    end
end
end
resids=-resids;

save('procresids_hurs.mat','resids','std20yr');

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

%figure(5);
%suptitle('Pct of models where abs(hurs)<1 std');
%hold on;
for k=1:6;
%    subplot(3,2,k);
%    contourf(transpose(squeeze(errormap_hurs(:,:,k))));
%    colorbar;
%    caxis([0 100]);
%    title(scenarios{k});
    writebinary(squeeze(errormap_hurs(:,:,k)),['annualresidual_hurs_' scenarios{k} '.bin'],-1e30);
end

erroraggregate_hurs=zeros(360,181);
for k1=1:360;
for k2=1:181;
    K=find(errormap_hurs(k1,k2,:)>90);
    erroraggregate_hurs(k1,k2)=length(K);
end
end

writebinary(erroraggregate_hurs,'annualaggregate_hurs.bin',-1e30);

%figure(6);
%contourf(transpose(erroraggregate_hurs));
%colorbar;

%minresid=squeeze(min(resids,[],3));
%maxresid=squeeze(max(resids,[],3));
%meanresid=squeeze(nanmean(resids,3));
%r_over_s=(maxresid-minresid)./squeeze(std(std20yr,0,3,'omitnan'));
%figure(16);
%hold on;
%for k=1:6;
%    subplot(6,3,(k-1)*3+1);
%    contourf(transpose(squeeze(meanresid(:,:,k))));
%    colorbar;
%    title('mean');
%    ylabel(scenarios{k});
%    subplot(6,3,(k-1)*3+2);
%    contourf(transpose(squeeze(maxresid(:,:,k))));
%    colorbar;
%    title('max across models');
%    ylabel(scenarios{k});
%    subplot(6,3,(k-1)*3+3);
%    contourf(transpose(squeeze(minresid(:,:,k))));
%    colorbar;
%    title('min across models');
%    ylabel(scenarios{k});
%    writebinary(squeeze(meanresid(:,:,k)),['annualmean_hurs_' scenarios{k} '.bin'],-1e30);
%    writebinary(squeeze(minresid(:,:,k)),['annualmin_hurs_' scenarios{k} '.bin'],-1e30);
%    writebinary(squeeze(maxresid(:,:,k)),['annualmax_hurs_' scenarios{k} '.bin'],-1e30);
%end
%figure(116);
%for k=1:6;
%    subplot(3,2,k);
%    contourf(transpose(squeeze(r_over_s(:,:,k))));
%    colorbar;
%    title(scenarios{k});
%    writebinary(squeeze(r_over_s(:,:,k)),['annualrs_hurs_' scenarios{k} '.bin'],-1e30);    
%end
