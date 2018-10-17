# Find PDB file binding sites description
# By: Li Zhixin
import os,sys
Amino=("GLY","ALA","VAL","LEU","ILE","PRO","PHE","TYR","TRP","SER","THR","CYS","MET","ASN","GLN","ASP","GLU","LYS","ARG","HIS")  


def writeSiteLine(fpdb,line,outputf):
    for i in range(0,49,11):
        cut=line[i:i+11]
        a=cut[0:3]
        if a in Amino:
            b=cut[3:5].strip()
            c=cut[5:].strip()
            #print("%s|%s|%s|%s"%(fpdb[-8:-4],b,c,a))
            outputf.write("%s|%s|%s|%s\n"%(fpdb[-8:-4],b,c,a))


def readBinding(path,ligndType,fpdb,outputf):
    file=open(path+'\\'+fpdb)
    #file=open(fpdb)
    line=file.readline()
    IDENTIFIER_LIST=[]
    IDENTIFIER=''
    while line:
        if line.find("SITE_IDENTIFIER:")>0:
            IDENTIFIER=line.split()[-1]
        elif line.upper().find("BINDING SITE FOR RESIDUE %s"%ligndType)>0:
            IDENTIFIER_LIST.append(IDENTIFIER)
        elif line[:4]=="SITE":
            break
            #print("%s|%s|%s"%(fpdb[-8:-4],line[57:59].strip(),line[59:65].strip()))
            #outputf.write("%s %s %s\n"%(fpdb[-8:-4],line[57:59].strip(),line[59:65].strip()))
        line=file.readline()
    while line[:4]=="SITE":
        if line[11:14] in IDENTIFIER_LIST:
            writeSiteLine(fpdb,line[18:61],outputf)
            #print("%s:%s\n"%(fpdb[-8:-4],line[18:61]))
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
            if count%100==0:
                print(".",end='')
            readBinding(path,ligndType,f,outputf)
    outputf.close()


allFiles=os.listdir()
for dirr in allFiles:
    if len(dirr)==3:
        dealWithType(dirr,dirr)

