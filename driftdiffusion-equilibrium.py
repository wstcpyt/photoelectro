##       1D Drift Diffusion Model 

import numpy as np

#% Defining the Fundamental and Material Constants %

q     = 1.602E-19# C or [J/eV]
kb    = 1.38E-23# [J/K]
eps   = 1.05E-12*13.6/11.7# This includes the eps  = 11.7 for Si [F/cm]
T     = 300# [K]
NC = 2.2e18
NV = 1.8e19
Eg  = 1.15 # [eV]
##ni    = 1.5E10# Intrinsic carrier concentration [1/cm^3]
ni = np.sqrt(NC*NV*np.exp(-Eg*q/(kb*T)))
Vt    = kb*T/q# [eV]
RNc   = 2.8E19# 
TAUN0 = 0.1E-6# Electron SRH life time
TAUP0 = 0.1E-6# Hole SRH life time
mun0   = 100# Electron Mobility in cm2/V-s
mup0   = 25# Hole Mobility in cm2/V-s
                        
                        
dEc = Vt*np.log(RNc/ni);

#% Define Doping Values %

Na = 1E16#; 1E18            % [1/cm^3]
Nd = 1E17#;   1E14          % [1/cm^3]

#% Calculate relevant parameters for the simulation %

Vbi = Vt*np.log(Na*Nd/(ni*ni))#;
W   = np.sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd)) #    % [cm]
Wn  = W*np.sqrt(Na/(Na+Nd))                 #   % [cm]
Wp  = W*np.sqrt(Nd/(Na+Nd))                  #  % [cm]
Wone = np.sqrt(2*eps*Vbi/(q*Na))              # % [cm]
E_p = q*Nd*Wn/eps                          # % [V/cm]
Ldn = np.sqrt(eps*Vt/(q*Nd))#;
Ldp = np.sqrt(eps*Vt/(q*Na))#;
Ldi = np.sqrt(eps*Vt/(q*ni))


#% Calculate relevant parameters in an input file %

#% Write to a file
##np.savetxt('input_params.txt', np.array([ Na, Nd, Vbi, W, Wn, Wp, E_p, Ldn, Ldp]))

#%Material_Constants    %Define some material constants
      
#% Setting the size of the simulation domain based 
#% on the analytical results for the width of the depletion regions
#% for a simple pn-diode %


x_max = 0#;
if(x_max < Wn):
    x_max = Wn#;

if(x_max < Wp):
    x_max = Wp

x_max = 20*x_max

#% Setting the grid size based on the extrinsic Debye lengths %

dx = Ldn#;
if(dx > Ldp):
    dx=Ldp

dx = dx/20

x_max = 2.5e-4



#% Calculate the required number of grid points and renormalize dx %

n_max = x_max/dx
n_max = np.round(n_max)

dx = dx/Ldi#;    % Renormalize lengths with Ldi

#% Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni %

dop=np.zeros(n_max)
for i in range(0,np.int(n_max)):
    if(i <= n_max/2):
        dop[i] = - Na/ni
    elif(i > n_max/2):
        dop[i] = Nd/ni

##% Initialize the potential based on the requirement of charge
##% neutrality throughout the whole structure

fi = np.zeros(n_max)
n = np.zeros(n_max)
p =np.zeros(n_max)

for i in range(0, np.int(n_max)):
    zz = 0.5*dop[i]
    if(zz > 0):
        xx = zz*(1 + np.sqrt(1+1/(zz*zz)))
    elif(zz <  0):
        xx = zz*(1 - np.sqrt(1+1/(zz*zz)))

    fi[i] = np.log(xx)
    n[i] = xx
    p[i] = 1/xx



delta_acc = 1E-5#;               % Preset the Tolerance

                                                               
##EQUILIBRIUM  SOLUTION PART                    %%

##%(A) Define the elements of the coefficient matrix for the internal nodes and
##%    initialize the forcing function

dx2 = dx*dx
a = np.zeros(n_max)
c = np.zeros(n_max)
b = np.zeros(n_max)
f = np.zeros(n_max)

for i in range( 0, np.int(n_max)):
    a[i] = 1/dx2;
    c[i] = 1/dx2;
    b[i] = -(2/dx2+np.exp(fi[i])+np.exp(-fi[i]));
    f[i] = np.exp(fi[i]) - np.exp(-fi[i]) - dop[i] - fi[i]*(np.exp(fi[i])+np.exp(-fi[i]))



##%(B) Define the elements of the coefficient matrix and initialize the forcing
##%    function at the ohmic contacts 

a[0] = 0
c[0] = 0
b[0] = 1
f[0] = fi[1]
a[n_max-1] = 0
c[n_max-1] = 0
b[n_max-1] = 1
f[n_max-1] = fi[n_max-1]

##%(C)  Start the iterative procedure for the solution of the linearized Poisson
##%     equation using LU decomposition method:

flag_conv = 0#;		           % convergence of the Poisson loop
k_iter= 0#;
alpha=np.zeros(n_max)
beta=np.zeros(n_max)
delta=np.zeros(n_max)
v=np.zeros(n_max)
while(flag_conv!=1):
    k_iter = k_iter + 1
    
    alpha[0] = b[0]
    for i in range(1,np.int(n_max)):
        beta[i]=a[i]/alpha[i-1]
        alpha[i]=b[i]-beta[i]*c[i-1]
    
##% Solution of Lv = f %    

    v[0] = f[0]
    for i in range( 1,np.int(n_max)):
        v[i] = f[i] - beta[i]*v[i-1]

     
##% Solution of U*fi = v %    

    temp = v[n_max-1]/alpha[n_max-1]
    delta[n_max-1] = temp - fi[n_max-1]
    fi[n_max-1]=temp
    for i in range((np.int(n_max)-2),0, -1):#       %delta%
        temp = (v[i]-c[i]*fi[i+1])/alpha[i]
        delta[i] = temp - fi[i]
        fi[i] = temp
    
    delta_max = 0
    for i in range( 0, np.int(n_max)):
        xx = abs(delta[i])
        if(xx > delta_max):
            delta_max=xx

#        %sprintf('delta_max = %d',delta_max)      %'k_iter = %d',k_iter,'
        
 
 #    %delta_max=max(abs(delta));
   
##% Test convergence and recalculate forcing function and 
##% central coefficient b if necessary
##    

    if(delta_max < delta_acc):
        flag_conv = 1
    else:
        for i in range( 1, np.int(n_max)-1):
            b[i] = -(2/dx2 + np.exp(fi[i]) + np.exp(-fi[i]))
            f[i] = np.exp(fi[i]) - np.exp(-fi[i]) - dop[i] - fi[i]*(np.exp(fi[i]) + np.exp(-fi[i]))

xx1 = np.zeros(n_max)
Ec = np.zeros(n_max)
ro = np.zeros(n_max)
el_field1 = np.zeros(n_max)
el_field2 = np.zeros(n_max)

xx1[0] = dx*1e4#;
for i in range( 1,np.int(n_max)-1):
    Ec[i] = dEc - Vt*fi[i]#;     %Values from the second Node%
    ro[i] = -ni*(np.exp(fi[i]) - np.exp(-fi[i]) - dop[i])
    el_field1[i] = -(fi[i+1] - fi[i])*Vt/(dx*Ldi)
    el_field2[i] = -(fi[i+1] - fi[i-1])*Vt/(2*dx*Ldi)
    n[i] = np.exp(fi[i])
    p[i] = np.exp(-fi[i])
    xx1[i] = xx1[i-1] + dx*Ldi*1e4


Ec[0] = Ec[1]
Ec[n_max-1] = Ec[n_max-2]
xx1[n_max-1] = xx1[n_max-2] + dx*Ldi*1e4
el_field1[0] = el_field1[1]
el_field2[0] = el_field2[1]
el_field1[n_max-1] = el_field1[n_max-2]
el_field2[n_max-1] = el_field2[n_max-2]
nf = n*ni
pf = p*ni
ro[0] = ro[1]
ro[n_max-1] = ro[n_max-2]

import matplotlib.pyplot as plt

plt.plot(xx1, Vt*fi,'r')
plt.xlabel('x [um]')
plt.ylabel('Potential [eV]')
plt.title('Potential vs Position - at Equilibrium')
plt.show()

plt.plot(xx1, el_field1,'r')

plt.plot(xx1, el_field2,'r')
plt.xlabel('x [um]')
plt.ylabel('Electric Field [V/cm]')
plt.title('Field Profile vs Position - at Equilibrium')
plt.show()

plt.plot(xx1, nf,'g', label='n')
##semilogy(xx1, nf,'g','LineWidth',2)

plt.plot(xx1, pf,'r', label='p')
##semilogy(xx1, pf,'r','LineWidth',2)
plt.xlabel('x [um]');
plt.ylabel('Electron & Hole Densities [1/cm^3]');
plt.title('Electron & Hole Densities vs Position - at Equilibrium');
plt.legend()
plt.yscale('log')
##%axis([0 6.75 0 10.2e17])
plt.show()

##%plot(xx1, ro,'r','LineWidth',2)
plt.plot(xx1, q*ro,'r')
plt.xlabel('x [um]')
##%ylabel('Total Charge Density [1/cm^3]');
plt.ylabel('Total Charge Density [C/cm^3]')
plt.title('Total Charge Density vs Position - at Equilibrium')
##%axis([0.5 5 -3e17 8e17])
plt.show()

plt.plot(xx1, Ec,'r')
plt.xlabel('x [um]')
##%ylabel('Total Charge Density [1/cm^3]');
plt.ylabel('Conduction Band Energy (eV)')
plt.title('Conduction Band vs Position - at Equilibrium')

plt.show()
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%                 END OF EQUILIBRIUM  SOLUTION PART                    %%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

