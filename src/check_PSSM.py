# this scipt check the completeness of pssm result
# first, it check the list of the fasta
# second, it find pssm in two path
# if found, copy pssm to .\pssm and delete fasta file
# else, rename *.fasta to *.pssm.error
# By: Li Zhixin, Mar 21, 2018.
from __future__ import print_function
import os,sys



def checkPSSMFinish(pssm):
    f=open(pssm)
    line=f.readline()
    lastline=''
    while line:
        lastline=line
        line=f.readline()
    if lastline[:10]=='PSI Gapped':
        return True
    return False

def findPSSM(fasta):
    pssmname=fasta.split('.')[0]+'.pssm'
    if os.path.exists('pssm\\'+pssmname) and checkPSSMFinish('pssm\\'+pssmname):
        print('OK in \\pssm')
        os.remove(fasta)
    elif os.path.exists('..\\..\\pssm\\'+pssmname)and checkPSSMFinish('..\\..\\pssm\\'+pssmname):
        print('OK in ..\\pssm')
        os.system('copy {} {}'.format('..\\..\\pssm\\'+pssmname,'pssm\\'))
        os.remove(fasta)
    else:
        print('Error file {}'.format(fasta))
        os.system('rename {} {}'.format(fasta,pssmname+'.error'))
        


pathnow=sys.path[0]
os.chdir(pathnow)
allFiles=os.listdir(pathnow)
try:
    os.mkdir('pssm')
except WindowsError:
    pass
for fasta in allFiles:
    if fasta[-5:]!='fasta': continue
    findPSSM(fasta)
input('finish')
