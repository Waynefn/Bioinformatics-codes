# Find binding sites by pdb info.
# Please use version 2.
# By: Li Zhixin
import os,sys

def readBinding(path,ligndType,fpdb,outputf):
    file=open(path+'\\'+fpdb)
    line=file.readline()
    while line:
        if line.find("BINDING SITE FOR RESIDUE %s"%ligndType)>0:
            #print("%s|%s|%s"%(fpdb[-8:-4],line[57:59].strip(),line[59:65].strip()))
            outputf.write("%s %s %s\n"%(fpdb[-8:-4],line[57:59].strip(),line[59:65].strip()))
        line=file.readline()

def dealWithType(path,ligndType):
    print("\nProcessing "+path,end='')
    allFiles=os.listdir(path)
    total=len(allFiles)
    count=0
    outputf=open(path+'\\'+ligndType+"_binding.txt",'w')
    outputf.write(ligndType+'\n')
    for f in allFiles:
        if f[-3:]=='pdb':
            count+=1
            if (count*10)%total==0:
                print(".",end='')
            readBinding(path,ligndType,f,outputf)
    outputf.close()


allFiles=os.listdir()
for dirr in allFiles:
    if len(dirr)==3:
        dealWithType(dirr,dirr)
