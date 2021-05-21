######## INPUT ############################################################
#xyzfile="xmolout-orig.xyz"
xyzfile=“xyz.xyz"
firstframe= int(raw_input("Please enter first frame: "))
lastframe= int(raw_input("Please enter last frame: "))
selectframes=range(firstframe,lastframe)
#selectframes=range(550,600) 
#places an atom or the mid-point of a group of atoms in the center
centerindex=[210] 
#following values obtained from limiting distances where jmol gives a bond:
#i.e. dX is half the threshold distance for the X-X bond
dH=0.91/2.
dO=1.81/2.
dSi=2.85/2.
dNa=2.39/2.
###### END OF INPUT #####################################################

########## IMPORTS #######################################################
from numpy import *
from os import system
from copy import *


#########################################################################
def CHANGESTRING(file,oldstring,newstring):
 """
 Change a certain string inside an existing file
 """
 st="sed 's/"+oldstring+"/"+newstring+"/g' "+\
     file+" > XXX.tmp;mv XXX.tmp "+file
 system(st)

#########################################################################
def DISTANCE(c1,c2,L):
  """
  Miniumum distance between coordinates c1 and c2 in periodic system
  """
  vector=c1-c2
  vector-=L*around(vector/L)
  d=sqrt(sum(vector*vector))
  return d

#########################################################################
def EXTRACTNEIGHBORSFROMLIST(atom,left,cutoffs,L):
  """
  Check neighbors of atom in list left. Remove them from left and return them
  as list extract.
  """
  indexleft=0
  extract=[]
  c1,el1=atom[0],atom[1]
  cutoff1=cutoffs[el1]
  while indexleft<len(left):
    secatom=left[indexleft]
    c2,el2=secatom[0],secatom[1]
    cutoff2=cutoffs[el2]
    d=DISTANCE(c1,c2,L)
    if d<cutoff1+cutoff2:
      extract+=[secatom]
      del left[indexleft]
    else:
      indexleft+=1
  return extract,left

#########################################################################
def MOLECLIST(atomlist,L,cutoffs):
 """
 Deliver back a list of molecules by checking with atoms are connected 
 """
 moleclist=[]
 left=deepcopy(atomlist)
 while len(left)>0:
   mol=[]
   mol+=[left[0]]
   del left[0]
   iat=0
   while iat<len(mol):
     atom=mol[iat]
     neighbors,left=EXTRACTNEIGHBORSFROMLIST(atom,left,cutoffs,L)
     mol+=neighbors
     iat+=1
   moleclist+=[mol] 
 return moleclist 

#########################################################################
def MIRRORCOORDINATES(mol,L,N):
  mirrors=[]
  for at1 in mol:
    c1,el1,index1=at1[0],at1[1],at1[2]
    translations=[]
    for at2 in mol: 
      c2=at2[0]
      vector=c2-c1
      trans=list(around(vector/L)) 
      if trans not in translations:
        translations+=[trans]
        c3=c1+trans*L
        index3=index1
        if trans!=[0.0,0.0,0.0]: index3+=N #change index indicating a mirror-image
        at3=[c3,el1,index3]
        mirrors+=[at3]
  return mirrors 

#########################################################################
def WRITEMOLECLIST(g,moleclist,counter):
  """
  Write info about which atoms form molecules
  """
  g.write("********** counter="+str(counter)+" ****************\n")
  for mol in moleclist:
    for at in mol:
      g.write(at[1])
    for at in mol:
      g.write(" "+str(at[2]))
    g.write("\n")

#########################################################################
def WRITEMIRROR2MOV(mirrorlist,framecount,f):
  """
  reorder atoms and write to file
  """
  for m in mirrorlist:
    m[0]=list(m[0])
  listatindex=map(lambda x: x[2],mirrorlist)
  reorderedmirror=zip(*sorted(zip(listatindex,mirrorlist)))[1]
  f.write(str(len(mirrorlist))+"\n")
  f.write(str(framecount)+"\n")
  for at in reorderedmirror:
    f.write(at[1]+" "+str(at[0]).replace("[","").replace("]","").replace(",","")+" "+str(at[2])+"\n") 

#########################################################################
def WRITEFRAME(N,f,g,counter,L,strcoordinates,cutoffs,centerindex):
 """
 write the frame 
 """
 #determine center
 shift=array([0.,0.,0.])
 for i in centerindex:
   x,y,z=strcoordinates[i-1].split()[1:4]
   xyz=array([float(x),float(y),float(z)])
   shift+=xyz
 shift/=max(1,len(centerindex))
 shift-=L/2.
 atomlist=[] 
 for atindex in range(N):
    element,x,y,z=strcoordinates[atindex].split()[0:4]
    xyz=array([float(x),float(y),float(z)])
    xyz-=shift
    xyz-=L*floor(xyz/L)
    atom=[xyz,element,atindex+1]
    atomlist+=[atom]
 moleclist=MOLECLIST(atomlist,L,cutoffs)
 WRITEMOLECLIST(g,moleclist,counter)
 mirrorlist=[]
 for mol in moleclist:
   mol2=MIRRORCOORDINATES(mol,L,N)
   mirrorlist+=mol2
 WRITEMIRROR2MOV(mirrorlist,framecount,f)

######## MAIN PROGRAM ###################################
f=open(xyzfile,"r");lines=f.readlines();f.close()
N=int(lines[0])
words=lines[1].split()
L=array([13,13,13])
nframes=len(lines)/(N+2)

print "number of atoms:", N
print "periodic box", L
print "number of frames", nframes
cutoffs={"H":dH,"O":dO,"Si":dSi,"Na":dNa}
f=open("movie.xyz","w");g=open("moleclist.txt","w")
lcount=0;framecount=0
for l in range(nframes):
  if framecount in selectframes:
    WRITEFRAME(N,f,g,framecount,L,lines[lcount+2:lcount+2+N],cutoffs,centerindex)
  framecount+=1 
  lcount=lcount+2+N
f.close();g.close()

st="cp jmol.base jmol.run";print st;system(st)
#Adjust jmol.run
newstring=str(L).replace("[","").replace("]","").replace(",","")
CHANGESTRING("jmol.run","boxdimensions",newstring)
newstring=""
for aa in cutoffs.keys():
  for bb in cutoffs.keys():
    newstring+="connect "+str(cutoffs[aa]+cutoffs[bb])+" (_"+aa+") (_"+bb+")\\\n "
CHANGESTRING("jmol.run","connections",newstring)
