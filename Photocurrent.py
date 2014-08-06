__author__ = 'yutongpang'
#internal quantum efficenc
import numpy as np
import IOEF
import math
import constant as cons
class Photocurrent(IOEF.IOEF):
    theta1=1.0 #is the quantum efficiency of the exciton generation
    theta2=1.0 #is the efficiency of the exciton dissociation at the interface
    D=8*10**-8 #is the diffusion constant
    t=2.9*10**-10 #is the mean life time of the exciton
    N=2.365*10**18 #incident photon flux for 1W/m**2*s
    Beta=1/(np.sqrt(D*t))
    #wait for setting
    rjp=0
    tjp=0
    rjpp=0
    tjpp=0
    layer=0
    dthickness=0
    ###############
#class method

    def initialsettingTransferMatrix(self):
        Sp=self.Sjpfun(self.layer)
        Spp=self.Sjppfun(self.layer)
        self.rjp=Sp[1,0]/Sp[0,0]
        self.tjp=1/Sp[0,0]
        self.rjpp=Spp[1,0]/Spp[0,0]
        self.tjpp=1/Spp[0,0]


    def Tjfunc(self):
        phi0=self.phi0
        n=self.id
        qj=(n[self.layer][0]**2-np.sin(phi0)*n[0][0]**2)**0.5
        yipsinon=2*math.pi*qj/self.wavelength
        tjplus=(np.absolute(self.tjp/(1-self.rjp*self.rjpp*np.exp(1j*2*yipsinon*self.dthickness))))**2
        Tj=(self.id[self.layer][0].real/self.id[0][0].real)* tjplus
        return Tj

    def alphajfunc(self):
        kj=self.id[self.layer][0].imag
        alphaj=4*math.pi*kj/self.wavelength
        return  alphaj

    def deltajppfunc(self):
        deltajpp=np.arctan(self.rjpp.imag/self.rjpp.real)
        return deltajpp

    def rojppfunc(self):
        rojpp=np.absolute(self.rjpp)
        return rojpp

    def C1func(self):
        ropp=self.rojppfunc()
        alpha=self.alphajfunc()
        d=self.dthickness
        C1=ropp**2*np.exp(-2*alpha*d)
        return C1

    def C2func(self):
        yita=self.id[self.layer][0].real
        ropp=self.rojppfunc()
        alpha=self.alphajfunc()
        beta=self.Beta
        d=self.dthickness
        top=(beta**2-alpha**2)*2*ropp*np.exp(-alpha*d)
        bottom=(beta**2+(4*math.pi*yita/self.wavelength)**2)
        C2=top/bottom
        return C2

    def Afunc(self):
        yita=self.id[self.layer][0].real
        alpha=self.alphajfunc()
        deltapp=self.deltajppfunc()
        d=self.dthickness
        beta=self.Beta
        C1=self.C1func()
        C2=self.C2func()
        top=(np.exp(beta*d-np.exp(-alpha*d)))+C1*(np.exp(beta*d)-np.exp(alpha*d))+C2*(np.exp(beta*d)*np.cos(4*math.pi*yita*d/self.wavelength+deltapp)-np.cos(deltapp))
        bottom=(np.exp(-beta*d)-np.exp(beta*d))
        A=top/bottom
        return A

    def Bfunc(self):
        n=self.id
        qj=(n[self.layer][0]**2-np.sin(self.phi0)*n[0][0]**2)**0.5
        yita=self.id[self.layer][0].real
        alpha=self.alphajfunc()
        deltapp=self.deltajppfunc()
        d=self.dthickness
        beta=self.Beta
        C1=self.C1func()
        C2=self.C2func()
        top=(np.exp(-beta*d-np.exp(-alpha*d)))+C1*(np.exp(-beta*d)-np.exp(alpha*d))+C2*(np.exp(-beta*d)*np.cos(4*math.pi*yita*d/self.wavelength+deltapp)-np.cos(deltapp))
        bottom=(np.exp(-beta*d)-np.exp(beta*d))
        B=-top/bottom
        return B

#public method
    def Qjfunc(self,j,Ejx):
        n=self.id[j][0].real
        k=self.id[j][0].imag
        alphaj=4*math.pi*k/self.wavelength
        Qj=(cons.c*cons.yita0*alphaj*n*Ejx**2)/2
        return Qj
#Jphoton at x=0
    def JPhotonfunc0(self):
        self.initialsettingTransferMatrix()
        yita=self.id[self.layer][0].real
        alpha=self.alphajfunc()
        T=self.Tjfunc()
        beta=self.Beta
        A=self.Afunc()
        B=self.Bfunc()
        C1=self.C1func()
        C2=self.C2func()
        d=self.dthickness
        delta=self.deltajppfunc()
        theta=self.theta1*self.theta2
        JPhoton=cons.e*theta*alpha*T*self.N*(-beta*A+beta*B-alpha+alpha*C1+4*math.pi*yita*C2*np.sin(4*math.pi*yita*d/self.wavelength+delta)/self.wavelength)/(beta**2-alpha**2)
        return JPhoton
    def JPhotonfuncd(self):
        yita=self.id[self.layer][0].real
        alpha=self.alphajfunc()
        T=self.Tjfunc()
        beta=self.Beta
        A=self.Afunc()
        B=self.Bfunc()
        C1=self.C1func()
        C2=self.C2func()
        d=self.dthickness
        delta=self.deltajppfunc()
        test1=np.exp(-beta*d)
        test2=np.exp(beta*d)
        theta=self.theta1*self.the2
        JPhoton=cons.e*theta*alpha*T*self.N*(beta*A*np.exp(-beta*d)-beta*B*np.exp(beta*d)+alpha*np.exp(-alpha*d)-alpha*C1*np.exp(alpha*d)-4*math.pi*yita*C2*np.sin(delta)/self.wavelength)/(beta**2-alpha**2)
        return JPhoton








