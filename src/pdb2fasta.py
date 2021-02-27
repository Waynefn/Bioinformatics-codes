# Convert pdb to fasta with start index and X (Non/Any amino acid)
# Change: Consider HETATM tag, fix bugs
# By: Li Zhixin

import os,sys
lib3=("GLY","ALA","VAL","LEU","ILE","PRO","PHE","TYR","TRP","SER","THR","CYS","MET","ASN","GLN","ASP","GLU","LYS","ARG","HIS")    
lib1=('G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','D','E','K','R','H') 
    
def three2one(three):
    global lib3,lib1
    if three in lib3:
        return lib1[lib3.index(three)]
    return 'X'
def one2three(one):
    global lib3,lib1
    return lib3[lib1.index(one)]
fasta=open('fasta.fasta','w')

allPDB=os.listdir()
for pp in allPDB:
    if pp[-3:]!='pdb': continue
    name=pp[:4]
    pdb=open(pp)
    line=pdb.readline()
    chid=''
    laststart=134124125
    count=0
    while line:
        if line[:4]=='ATOM' or line[:6]=='HETATM':
            ch=line[21]
            start=line[22:26]
            amno=line[17:20]
            if ch!=chid:
                chid=ch
                count=0
                print('.',end='')
                fasta.write('\n>{pdb}:{ch}|{offset}\n'.format(pdb=name,ch=ch,offset=int(start)))

            if start!=laststart:
                count+=1
                fasta.write(three2one(amno))
                if count==80:
                    fasta.write('\n')
                    count=0
            laststart=start
        
        line=pdb.readline()
fasta.close()
        
