import os
import sys

# Subroutine that converts read in text file lines into a list
def linestolist(linesraw,linelist):
    delimiter = ''
    for k in range(len(linelist)):
        temp = linelist[k]
        temp4 = list(temp)
        for j in range(len(temp4)):
            if temp4[j]=='\n':
                del temp4[j]
        temp5 = delimiter.join(temp4)
        temp6 = temp5.split(',')
        linelist[k] = temp6

os.system('ls *.bin > thefiles.txt')
filesfile=open('thefiles.txt')
fileslines=filesfile.readlines()
fileslineslist=list(fileslines)
linestolist(fileslines,fileslineslist)
filesfile.close()
os.system('rm thefiles.txt')

for j in range(len(fileslineslist)):
    whichfile=fileslineslist[j][0]
    name=whichfile[:-4]
    
    writetext=["DSET "+name+".bin\n","OPTIONS big_endian\n","UNDEF -1.e30\n","TITLE "+name+"\n","XDEF 360 LINEAR -0.5 1\n","YDEF 181 LINEAR -90 1\n","ZDEF 1 LEVELS 1000\n","TDEF 1 LINEAR 01jan0000 1mo\n","VARS 1\n","var 0 t,y,x var\n","ENDVARS\n"]
    fname=name+".ctl"
    
    fid=open(fname,'w')
    fid.writelines(writetext)
    fid.close()