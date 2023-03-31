% sample.m
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com) and Abigail Snyder
% Last updated 20 September 2022
%
% This script will produce some sample figures of pattern scaling for Figure 1.

% Preamble begins here
clear all;

% Directories where the various files are.  Feel free to change this.
residuals_dir='/Users/bkravitz/Documents/github repositories/cmip6_patterns/outputs/resids_tgavs';
patterns_dir='/Users/bkravitz/Documents/github repositories/cmip6_patterns/outputs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN PART OF THE SCRIPT BEGINS HERE %%%%%
%%%%%%%%%%%% EDIT AT YOUR OWN RISK %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname1='CESM2-WACCM_ssp585_pr_annual_pattern_resids.nc';
fname2='CESM2-WACCM_ssp585_pr_annual_pattern.nc';
fname3='CESM2-WACCM_ssp585_ensemble_avg_tgav.nc';

cd(residuals_dir);
ncid=netcdf.open(fname1,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'pr');
resid=netcdf.getVar(ncid,varid);
netcdf.close(ncid);

cd(patterns_dir);
ncid=netcdf.open(fname2,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'slope');
slope=netcdf.getVar(ncid,varid);
varid2=netcdf.inqVarID(ncid,'intercept');
intercept=netcdf.getVar(ncid,varid2);
netcdf.close(ncid);

cd(residuals_dir);
ncid=netcdf.open(fname3,'NC_NOWRITE');
varid=netcdf.inqVarID(ncid,'tas');
tgav=netcdf.getVar(ncid,varid);
netcdf.close(ncid);

% Reconstructing the pattern scaled timeseries
J=size(resid);
reconstructed_last20=zeros(J(1),J(2),20);
I=length(tgav);
c=1;
for j3=I-19:I;
    reconstructed_last20(:,:,c)=slope*tgav(j3)+intercept;
    c=c+1;
end

resid=-resid; % Residuals are output as ESM minus pattern scaling, and Ben thinks the opposite way

%figure(1);
%hold on;
%subplot(2,2,1);
%contourf(transpose(slope*86400));
%colorbar;
%title('Slope');
%subplot(2,2,2);
%contourf(transpose(intercept*86400));
%colorbar;
%title('Intercept');
%subplot(2,2,3);
%contourf(transpose(squeeze(mean(reconstructed_last20*86400,3))));
%colorbar;
%title('Generated timeseries from PS');
%subplot(2,2,4);
%contourf(transpose(squeeze(mean(resid(:,:,I-19:I)*86400,3))));
%colorbar;
%title('Residual');

writebinary(slope*86400,'sampleslope.bin',-1e30);
writebinary(intercept*86400,'sampleintercept.bin',-1e30);
writebinary(squeeze(mean(reconstructed_last20*86400,3)),'samplegenerated.bin',-1e30);
writebinary(squeeze(mean(resid(:,:,I-19:I)*86400,3)),'sampleresidual.bin',-1e30);
