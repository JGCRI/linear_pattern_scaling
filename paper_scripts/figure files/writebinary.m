function writebinary(variable,fname,missingvalue)

I=find(isnan(variable)==1);
variable(I)=missingvalue;
fid=fopen(fname,'wb');
fwrite(fid,variable,'single','b');
fclose(fid);