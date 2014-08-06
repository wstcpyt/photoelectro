__author__ = 'yutongpang'
#internal optical electric field
import numpy as np
import math
import pylab
#for TE
class IOEF():

    def __init__(self, id, phi0, wavelength, TME):
        self.id=id
        self.phi0=phi0
        self.wavelength=wavelength
        self.TME=TME

    def interfacematrix(self):
        Im=np.array([])
        n=self.id.size


#Transfer Matrix Sp
    def Sjpfun(self,j):
        Po=1;
        for v in range(1,j):
            Po0=self.Ijkfunc(v-1,v)
            Po1=self.Ljkfunc(v)
            Po=Po0*Po1*Po
        Sjp=Po*self.Ijkfunc(j-1,j)
        return Sjp
#Transfer matrix Spp
    def Sjppfun(self,j):
        m=self.id.shape[0]-2
        Po=1
        for v in range(j+1,m+1):
            Po0=self.Ijkfunc(v-1,v)
            Po1=self.Ljkfunc(v)
            Po=Po0*Po1*Po
        Sjpp=Po*self.Ijkfunc(m,m+1)
        return Sjpp
## calculatge Ijk, Ljk matrix
    def Ijkfunc(self,j,k):
        rjk=self.rjkfunc(j,k)
        tjk=self.tjkfunc(j,k)
        Ijk=np.matrix([[1,rjk],[rjk,1]])/tjk
        return Ijk

    def Ljkfunc(self,j):
        dj=self.id[j][1]
        W=self.wavelength
        yitaj=2*math.pi*self.qjfunc(j)/W
        Lj11=np.exp(yitaj*dj*-1j)
        Lj22=np.exp(yitaj*dj*1j)
        Lj=np.matrix([[Lj11,0],[0,Lj22]])
        return Lj

    def rjkfunc(self,j,k):
        qj=self.qjfunc(j)
        qk=self.qjfunc(k)
        nj=self.id[j][0]
        nk=self.id[k][0]
        if (self.TME=="TE"):
            rjk=(qj-qk)/(qj+qk)
            return rjk
        if (self.TME=="TM"):
            rjk=(qj*nk**2-qk*nj**2)/(qj*nk**2+qk*nj**2)
            return rjk



    def tjkfunc(self,j,k):
        qj=self.qjfunc(j)
        qk=self.qjfunc(k)
        nj=self.id[j][0]
        nk=self.id[k][0]
        if(self.TME=="TE"):
            tjk=2*qj/(qj+qk)
            return tjk
        if(self.TME=="TM"):
            tjk=2*nj*nk*qj/(qj*nk**2+qk*nj**2)
            return tjk

    def qjfunc(self,j):
        n=self.id
        phi0=self.phi0
        qj=(n[j][0]**2-np.sin(phi0)*n[0][0].real**2)**0.5
        return qj
    #electric field function

    def EFfunc(self,j,x):
        dj=self.id[j][1]
        W=self.wavelength
        yitaj=2*math.pi/W*self.qjfunc(j)
        top=self.Sjppfun(j)[0,0]*np.exp(yitaj*(dj-x)*-1j)+self.Sjppfun(j)[1,0]*np.exp(yitaj*(dj-x)*1j)
        bottom=self.Sjpfun(j)[0,0]*self.Sjppfun(j)[0,0]*np.exp(dj*yitaj*-1j)+self.Sjpfun(j)[0,1]*self.Sjppfun(j)[1,0]*np.exp(yitaj*dj*1j)
        Ej= top/bottom
        return Ej
    #plot function
    def plotEFfunc(self,layer):
        numberofpoints=int(self.id[layer][1]*10**9)
        Ej=np.empty([numberofpoints,])
        for x in range(0,numberofpoints):
            xi=x*10**-9
            Ej[x]=self.EFfunc(layer,xi)
        x=np.arange(0,numberofpoints)*10**-9
        pylab.plot(x,Ej**2)
        pylab.ylabel('|E|^2 (normalized)')
        pylab.xlabel('thickmess (m)')
        pylab.show()















