#!/usr/bin/python
#*********************************************************************************
# nwvSymmetry.py
# Created: Jan 8, 2013
# Last edited: Jan 25, 2013
# Bradley Kearney
# Calculates RMSDs and other values for virus capsid subunits.
#*********************************************************************************

import sys
import os
from pdb import *
import numpy as np
import math
import re
import SVD

def question(prompt="",valid=[]):
    flag=0
    while not flag:
        print prompt
        input1=raw_input()
        try:
            if not len(valid):
                return input1
            if input1[0] not in valid:
                print "Invalid input."
            else:
                return input1
        except:
            if len(input1)==0 and '' in valid:
                return ''
            else:
                print "Make a selection."

class Analyzer:
    
    def __init__(self, pdb_filenames):
        self.pdbs = []
        for filename in pdb_filenames:
            self.pdbs.append(PDB(filename))
        self.pdbs.sort(key=str)
        for i, pdb in enumerate(self.pdbs):
            pdb.number = i
            pdb.shortname = pdb.filename.split('/')[-1]
        self.angle=math.radians(180)
        self.d=[0,0,0]
        self.s=[0,0,0]
        self.set=0
        self.rms=0

    def center_of_mass(self, chains):
        masses={'C':12.01, 'N':14.01, 'O':16.01, 'S':32.01}
        backbone=["CA","N","C"]
        totalmass=0
        sumx=0
        sumy=0
        sumz=0
        for i in chains:
            for j in i:
                    totalmass+=masses[j['atom name'][0]]
                    sumx+=masses[j['atom name'][0]]*j['x']
                    sumy+=masses[j['atom name'][0]]*j['y']
                    sumz+=masses[j['atom name'][0]]*j['z']
                
        a=(sumx/totalmass)
        b=(sumy/totalmass)
        c=(sumz/totalmass)
        return [a,b,c]

    def get_center_of_mass(self,chains):
        return self.center_of_mass(chains)

    def getangle(self,R):
        return "%.2f"%np.rad2deg(math.acos(( R[0][0] + R[1][1] + R[2][2] - 1)/2))


    def getaxis(self):
        if not(self.set):
            return [0,0,0]
        angle=self.angle
        print angle
        d=np.array(self.d)
        s=np.array(self.s)
        b=math.tan(angle/2)*s
        print b
        c=(np.cross(b,d)-np.cross(b,np.cross(b,d)))/(np.dot(2*b,b))
        print "Distance from origin: %s"%c
        return c

    def getvector(self,R,flag=1,COM=[0,0,0]):
        R21=float(R[2][1])
        R12=float(R[1][2])
        R02=float(R[0][2])
        R20=float(R[2][0])
        R10=float(R[1][0])
        R01=float(R[0][1])
        xx = (R21 - R12)/math.sqrt(np.power((R21 - R12),2)+np.power((R02 - R20),2)+np.power((R10 - R01),2))
        yy = (R02 - R20)/math.sqrt(np.power((R21 - R12),2)+np.power((R02 - R20),2)+np.power((R10 - R01),2))
        zz = (R10 - R01)/math.sqrt(np.power((R21 - R12),2)+np.power((R02 - R20),2)+np.power((R10 - R01),2))
        xx=abs(xx)
        yy=abs(yy)
        zz=abs(zz)
        if flag:
            print [xx,yy,zz]
        return[xx,yy,zz]

    def printchain(self,chain,title="output/default.pdb",xyz=[0,0,0],COM=[0,0,0],opt=0,nowrite=0,RMSD=0,add=1,dorun=0):
        title='output/'+title
        rCOM=self.getaxis()
        RMSD=self.rms
        prev=0
        if not opt:
            prev=0
        
        prev=math.sqrt((rCOM[0]-COM[0])**2+(rCOM[1]-COM[1])**2+(rCOM[2]-COM[2])**2)
        current=prev
        counter=0
        prev=current
        rCOM[0]=rCOM[0]-xyz[0]
        rCOM[1]=rCOM[1]-xyz[1]
        rCOM[2]=rCOM[2]-xyz[2]
        current=math.sqrt((rCOM[0]-COM[0])**2+(rCOM[1]-COM[1])**2+(rCOM[2]-COM[2])**2)
        if(prev>current):
            #The axis is going in the negative direction, so propogate by subtracting.
            while(prev>=current): 
                prev=current
                rCOM[0]=rCOM[0]-xyz[0]
                rCOM[1]=rCOM[1]-xyz[1]
                rCOM[2]=rCOM[2]-xyz[2]
                current=math.sqrt((rCOM[0]-COM[0])**2+(rCOM[1]-COM[1])**2+(rCOM[2]-COM[2])**2)
        else:
            #The axis is going in the negative direction, so propogate by adding.
            prev=current
            rCOM[0]=rCOM[0]+xyz[0]
            rCOM[1]=rCOM[1]+xyz[1]
            rCOM[2]=rCOM[2]+xyz[2]
            current=math.sqrt((rCOM[0]-COM[0])**2+(rCOM[1]-COM[1])**2+(rCOM[2]-COM[2])**2)
            while(prev>=current):
                prev=current
                rCOM[0]=rCOM[0]+xyz[0]
                rCOM[1]=rCOM[1]+xyz[1]
                rCOM[2]=rCOM[2]+xyz[2]
                current=math.sqrt((rCOM[0]-COM[0])**2+(rCOM[1]-COM[1])**2+(rCOM[2]-COM[2])**2)
        file = open(str(title), 'w')
        clusterline = 'ATOM   1052  O   ORG O   1      35.155  48.405 170.828  1.00             \n' #Generic PDB line.
        temp = Patom(clusterline, 9998)
        temp['x']=float(COM[0])
        temp['y']=float(COM[1])
        temp['z']=float(COM[2])
        temp['B']=float(RMSD)
        file.write(str(temp))
        COM=rCOM
        for i in chain:
            i['B']=float(RMSD)
            file.write(str(i))
            file.write('\n')
        if math.isnan(COM[0]) or nowrite:
            file.close()
            return
        temp = Patom(clusterline, 9998)
        temp['x']=float(COM[0])
        temp['y']=float(COM[1])
        temp['z']=float(COM[2])
        temp['B']=float(RMSD)
        file.write(str(temp))
        i=0
        while i<250:
            temp=Patom(clusterline,9999+i)
            temp['x']=COM[0]+i*xyz[0]
            temp['y']=COM[1]+i*xyz[1]
            temp['z']=COM[2]+i*xyz[2]
            temp['B']=float(RMSD)
            file.write(str(temp))
            i=i+1
        i=0
        file.close()

    def backtocoord(self,chain,newcoords):
        chain2=deepcopy(chain)
        for i,j in zip(chain2,newcoords):
            i['x']=j[0]
            i['y']=j[1]
            i['z']=j[2]
        return chain2

    def symop(self,matrixA,matrixB,title=""):
        a=SVD.SVDSuperimposer()
        a.set(matrixA[1],matrixB[1])
        x=[]
        x.append(matrixA[0])
        x.append(matrixB[0])
        COM=self.get_center_of_mass(x)
        print "COMPARISON "+title
        print "Centroid: "
        print COM
        a.run()

        R=np.transpose(a.get_rot())
        print "Rotation Matrix:"
        print R
        print title+" Rotation Angle: "+self.getangle(R)
        angle=self.getangle(R)
        print title+" RMSD: "+a.get_rms()
        print '-----'

    def coordmatrix(self,chain):
        temp=[]
        for i in chain:
            temp.append([i['x'],i['y'],i['z']])
        return np.array((temp),'f')
        
    def main(self,filetype):
        self.set=1
        for pdb in self.pdbs:
            print "Autoloading detected PDB file: %s"%pdb
            print ", ".join([str(x) for x in pdb.list_chains()] )
            commonAA=pdb.commonAA(123,531)#Only use amino acids 123-531 for superimposing and computing RMSD.
        f = open ( 'input.txt' , 'r')
        l=[]
        matrices=[]
        matrices.append([])
        for line in f:
            if line.strip()=="":
                matrices.append(l)
                l=[]
            else:
                l.append(map(float,line.split(' ')))
        flag=1
        while flag:
            towrite1=[]
            matrix_static=[]
            matrix_moving=[]
            towrite2=[]
            counter1=1
            counter2=1
            
            input2=question("What subunit is the static reference? (eg 5A)")
            input2=input2.strip()
            matrixID2= int(re.findall(r'\d+',input2)[0])
            chainID2= re.findall(r'\D+',input2)[0].upper()
            input2=input2.upper()
            input1=question("What subunit is moving? (eg 1A)")
            input1=input1.strip()
            matrixID1=int(re.findall(r'\d+',input1)[0])
            chainID1= re.findall(r'\D+',input1)[0].upper()
            input1=input1.upper()
            for pdb in self.pdbs:
                for atom in pdb.backbone:
                    if atom['chain']==chainID1 and atom['res num'] in commonAA:
                        a=matrices[matrixID1]*np.transpose(atom.coord())
                        temp=Patom(str(atom),1)
                        temp['x']=float(a[0])
                        temp['y']=float(a[1])
                        temp['z']=float(a[2])
                        temp['line num']=counter1
                        counter1+=1
                        towrite1.append(temp)
                    if atom['chain']==chainID2 and atom['res num'] in commonAA:
                        a=matrices[matrixID2]*np.transpose(atom.coord())
                        temp=Patom(str(atom),1)
                        temp['x']=float(a[0])
                        temp['y']=float(a[1])
                        temp['z']=float(a[2])
                        temp['line num']=counter2
                        counter2+=1
                        towrite2.append(temp)
            file = open("output/%s%s.pdb"%(filetype,input1),'w')
            for i in towrite1:
                file.write(str(i)+'\n')
            file.close()
            file = open("output/%s%s.pdb"%(filetype,input2),'w')
            for i in towrite2:
                file.write(str(i)+'\n')
            file.close()
            matrix_static=self.coordmatrix(towrite1)
            matrix_moving=self.coordmatrix(towrite2)
            a=SVD.SVDSuperimposer()
            static=[towrite1,matrix_static]
            moving=[towrite2,matrix_moving]
            self.symop(static,moving,"%s-to-%s"%(input1,input2))
            continue_check=question("Continue? (Y/n)",['','y','Y','n','N'])
            if continue_check not in ['','y','Y']:
                flag=0

def run():
    argv = sys.argv[1:]
    pdb_filenames=[]
    filenames = os.listdir(os.getcwd())
    for f in filenames:
        if(f[-4:] == '.pdb'):
            pdb_filenames.append(f)
    a = Analyzer(pdb_filenames)
    reader=pdb_filenames[0]
    if 'capsid' in reader:
        if 'procapsid' in reader:
            input1="Procapsid"
        else:
            input1="Capsid"
    else:
        input1=question("Is this a (c)apsid or (p)rocapsid structure?",['c','C','p','P'])
        if input1[0] in ['C','c']:
            input1="Capsid"
        else:
            input1="Procapsid"
    a.main(input1)
    return

if(__name__ == '__main__'):
    run()
        
        
