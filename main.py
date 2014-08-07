__author__ = 'yutongpang'
import IOEF
import Photocurrent
import numpy as np
import math
import pylab
import matplotlib.pyplot as plt
import constant


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
"""
x=np.arange(0,numberofpoints)*10**-9
pylab.plot(x,Ej**2, label="TE")
pylab.plot(x,EjM**2, label="TM")
pylab.ylabel('|E|^2 (normalized)')
pylab.xlabel('thickmess (m)')
pylab.legend(loc='upper right')
pylab.show()
"""
"""
#initate IQE class
activelayer=1;
p2=IQE.IQEfunc(Structurearray,0,wavelength)
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
    W=data[70:,0]
    n=data[70:,1]
    k=data[70:,2]
print(k[1])
#calculate Jphoton Spectrum
numberofwavelength=W.shape[0]
JPhoton=np.empty([numberofwavelength,])
JPhotond=np.empty([numberofwavelength,])
IPCE=np.empty([numberofwavelength,])
IPCEd=np.empty([numberofwavelength,])

for i in range(0,numberofwavelength):
    Structurearray=np.array([[1,500*10**-9],[n[i]+k[i]*1j,100*10**-9],[1.2,500*10**-9],[1,500*10**-9]])
    #initate IQE class
    activelayer=1
    p2=Photocurrent.Photocurrent(Structurearray,0,W[i]*10**-9,"TE")
    p2.layer=activelayer
    Nphoton=p2.Nfunc(W[i])
    p2.N=Nphoton
    p2.dthickness=100*10**-9
    JPhoton[i]=p2.JPhotonfunc0()
    JPhotond[i]=p2.JPhotonfuncd()
    IPCE[i]=p2.ipce0()
    IPCEd[i]=p2.ipced()
########################
"""
pylab.plot(W,IPCE,label='x=0')
pylab.plot(W,IPCEd,label='x=d')
pylab.ylabel('IPCE(%)')
pylab.xlabel('wavelength (nm)')
pylab.legend(loc='upper right')
pylab.show()
"""
"""
am15 = [('data/AM1.5.txt', 'experiment data')]
am15data = [ ( pylab.loadtxt(filename), label ) for filename, label in am15 ]
for data,label in am15data:
    am15W=data[0:,0]
    am15I=data[0:,1]
pylab.plot(am15W,am15I)
pylab.ylabel('Spectral Irradiance W m-2 nm -1')
pylab.xlabel('wavelength (nm)')
pylab.show()
"""
"""
pylab.plot(W,JPhoton,label='x=0')
pylab.plot(W,JPhotond,label='x=d')
pylab.ylabel('Jsc (A/m^2)')
pylab.xlabel('wavelength (nm)')
pylab.legend(loc='upper right')
pylab.show()
"""
"""
#Demonstrate how to do two plots on the same axes with different left right scales.
fig,ax1=plt.subplots()
lns1=ax1.plot(W,n,'b-',label=r'$\eta$')
ax1.set_xlabel('wavelength(nm)')
ax1.set_ylabel(r'$\eta$',color='b')
for t1 in ax1.get_yticklabels():
    t1.set_color('b')
ax2 = ax1.twinx()
lns2=ax2.plot(W,k,'r.',label=r'$\kappa$')
ax2.set_ylabel(r'$\kappa$',color='r')
for t1 in ax2.get_yticklabels():
    t1.set_color('r')
lns=lns1+lns2
labs=[l.get_label() for l in lns]
ax1.legend(lns,labs,loc='upper right')
plt.show()
"""

