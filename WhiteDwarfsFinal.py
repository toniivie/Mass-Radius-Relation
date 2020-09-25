import matplotlib.pyplot as plt
import numpy as np
import math as m
"""Code including both carbon-12 and iron-56 (both pure) for mass-radius relation graph and the three different white dwarfs of different compositions"""
    
N = 1000 #Where N is the number of steps the Runge-Kutta method will take.
rho_mass=1.6726231e-27 #constants (physical)
e_mass=9.1093897e-31
hbar=1.05457266e-34
G=6.67259e-11
c=2.99792458e8
C_12Ye=0.5
FeYe=0.464
C_12p0=9.79e8/C_12Ye
Fep0=9.79e8/FeYe
solarrad=6.95e8
solarmass=1.98e30
CR0=(9*m.pi*(C_12Ye**2)*(hbar**3)/(4*G*(rho_mass**2)*(e_mass**2)*c))**0.5        #R0 and M0 are functions for carbon(C_12) and iron(Fe)
CM0=(4*m.pi*(C_12p0)*(CR0**3))/3
FeR0=(9*m.pi*(FeYe**2)*(hbar**3)/(4*G*(rho_mass**2)*(e_mass**2)*c))**0.5
FeM0=(4*m.pi*(Fep0)*(FeR0**3))/3
totmassC_12=[]   #total mass for carbon
totradC_12=[] # total radius for carbon
totmassFe=[] #total mass for iron
totradFe=[] #total radius for iron

def dv(n, x, y):
    """defining the derivatives including the equation for the mass"""
    dyy=[0 for i in range(0,n+1)]
    if x==0.0:        #this checks to see if position is at the centre of star so that the density's derivative is 0, if not then the density derivative is equal to equation
        dyy[1]=-0.0
    else:
        dyy[1]=-y[1]*y[2]/((x**2)*gamma(y[1]))
    dyy[2]= 3*(x**2)*y[1]        #derivative equation for mass  
    return dyy

def RungeKutta(n, x, y, h):
    "Uses each step of solution for differential equation into final step 4th Runge-Kutta (takes from x to x+h)"
    y0=y[:]
    K1=dv(n, x, y)
    for i in range(1,n+1):                                         #only first ordinary differntial equations can be solved by Runge Kutta method  (RKM)
        y[i]=y0[i]+0.5*h*K1[i]
    K2=dv(n, x +0.5*h, y)
    for i in range(1,n+1):
        y[i]=y0[i]+h*(0.2071067811*K1[i]+0.2928932188*K2[i])
    K3=dv(n, x +0.5*h, y)
    for i in range(1,n+1):
        y[i]=y0[i]-h*(0.7071067811*K2[i]-1.7071067811*K3[i])
    K4=dv(n, x+h, y)
    for i in range(1,n+1):
        a=K1[i]+0.5857864376*K2[i]+3.4142135623*K3[i]+K4[i]
        y[i]=y0[i]+0.16666666667*h*a

    x+=h
    return (x,y)

def gamma(y):
    """Definition for gamma function with y variable"""
    gamma = ((y**(2/3))/(3*(1+(y**(2/3)))**0.5))
    return gamma  

ES=[0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000] # the different densities needing to be integrated over.
for t in ES:
    for i in np.arange(10*t, 100*t, t):     # this makes sure the intervals are go up in equal number of steps by a factor of 10
        rho=float(i) #this sets the initial density rho (the central density)
        y=[0, rho, 0.0] # Setting the intial boundary conditions
        x=0.0
        radiival=[0]    # creates the arrays which will hold the values for each Runge-Kutta method run
        massval=[0]
        densityvalues=[rho]
        for k in range(0, 100*N):
            (x,y)=RungeKutta(2,x,y,1.0/N) # this runs the Runge-Kutta program and then saves the values at each position in arrays
            radiival.append(x)
            massval.append(y[2])
            densityvalues.append(y[1])
            if y[1].imag !=0:              # this checks if the density becomes 0 and  represents the surface of the star
                totmassC_12.append(massval[-1]*CM0/solarmass)       # Picks out the  last massval value value storing in arrays with expression including mass and radius
                totradC_12.append(radiival[-1]*CR0/solarrad)        # This converts data into units of solar masses and solar radii
                totmassFe.append(massval[-1]*FeM0/solarmass) 
                totradFe.append(radiival[-1]*FeR0/solarrad)
                break                                #this break makes loop exits loop once the condition on the surface has been reached, the radius and mass are addedd into arrays

plt.figure(figsize=(10,10), frameon = True) # to print graph
plt.plot(totradC_12, totmassC_12, color = 'r', label = 'Pure Carbon-12')
plt.plot(totradFe, totmassFe, color = 'g', label = 'Pure Iron-56')
plt.errorbar(0.0124, 0.48, yerr = 0.02, xerr = 0.0005, label = '40 Eridani B', color = 'c')
plt.errorbar(0.0074, 1.053, yerr = 0.028, xerr = 0.0006, label = 'Sirius B', color = 'b')#plots astronomical data

plt.errorbar(0.0115, 0.50, yerr = 0.05, xerr = 0.0012, label = 'Stein 2051',color = 'k')
plt.xlabel('Radius of white dwarf (Solar Radii)')
plt.ylabel('Mass of white dwarf (Solar mass)')
plt.grid(which='both',color='black', linestyle ='--')
plt.title('Mass-Radius relationship for White Dwarfs')
plt.legend()
