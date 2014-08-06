__author__ = 'yutongpang'
import IOEF
import Photocurrent
import numpy as np
import math
import pylab


Structurearray=np.array([[1.45,500*10**-9],[1.8+0.2j,100*10**-9],[1.2,500*10**-9],[1,500*10**-9]])
wavelength=500*10**-9
p1=IOEF.IOEF(Structurearray,0,wavelength,"TE")

p1M=IOEF.IOEF(Structurearray,0,wavelength,"TM")
layer=1
a=math.pi

#p1.plotEFfunc(layer)
numberofpoints=int(Structurearray[layer][1]*10**9)
Ej=np.empty([numberofpoints,])
EjM=np.empty([numberofpoints,])
for x in range(0,numberofpoints):
    xi=x*10**-9
    Ej[x]=p1.EFfunc(layer,xi)
for x in range(0,numberofpoints):
    xi=x*10**-9
    EjM[x]=p1M.EFfunc(layer,xi)
x=np.arange(0,numberofpoints)*10**-9
pylab.plot(x,Ej**2, label="TE")
pylab.plot(x,EjM**2, label="TM")
pylab.ylabel('|E|^2 (normalized)')
pylab.xlabel('thickmess (m)')
pylab.legend(loc='upper right')
pylab.show()

"""
#initate IQE class
activelayer=1;
p2=IQE.IQEfunc(Structurearray,0,wavelength)
Sp=p1.Sjpfun(activelayer)
Spp=p1.Sjppfun(activelayer)
p2.rjp=Sp[1,0]/Sp[0,0]
p2.tjp=1/Sp[0,0]
p2.rjpp=Spp[1,0]/Spp[0,0]
p2.tjpp=1/Spp[0,0]
p2.layer=activelayer
p2.dthickness=100*10**-9
########################
#JPhton=p2.JPhotonfuncd()
#print(JPhton)
Qj=np.empty([100,])
for x in range(0,100):
    Qj[x]=p2.Qjfunc(activelayer,Ej[x]) #energy aborbed at x
#print(Qj)
"""
#import experiment data
list_of_files = [('data/p3htpcbmnew.txt', 'experiment data')]
datalist = [ ( pylab.loadtxt(filename), label ) for filename, label in list_of_files ]
for data,label in datalist:
    W=data[0:,0]
    n=data[0:,1]
    k=data[0:,2]
print(k[1])
#calculate Jphoton Spectrum
numberofwavelength=W.shape[0]
JPhoton=np.empty([numberofwavelength,])
for i in range(0,numberofwavelength):
    Structurearray=np.array([[1,500*10**-9],[n[i]+k[i]*1j,100*10**-9],[1.2,500*10**-9],[1,500*10**-9]])
    #initate IQE class
    activelayer=1
    p2=Photocurrent.Photocurrent(Structurearray,0,W[i]*10**-9,"TE")
    p2.layer=activelayer
    p2.dthickness=100*10**-9
    JPhoton[i]=p2.JPhotonfunc0()
########################
print(JPhoton)
pylab.plot(W,JPhoton)
pylab.ylabel('Jsc (A/m^2)')
pylab.xlabel('wavelength (nm)')
pylab.show()

