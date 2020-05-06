import matplotlib.pyplot as plt
import time
from tqdm import tqdm
# U=k\x\
import numpy as np
from numba import njit
f0=0.01 #0.0026
Omega1=1
Omega2=1.2
k=0.6
A=((3.*np.pi)/(4.*np.sqrt(2.)))
#   global Omega1, Omega2, f0, A
@njit
def ThetaDif(t,Theta,I):
    return  (2./3.)*A*(k**(2./3.))* (I ** (-1./3.))-A* f0 *(k**(-1./3.))*(I ** (-1./3.))*np.abs(np.cos(Theta))*(np.cos(Omega1*t)+np.cos(Omega2*t))
@njit
def IDif(t,Theta,I):
        return -A* f0 *(k**(-1./3.))*(I ** (2./3.))*np.sin(Theta)*(np.cos(Omega1*t)+np.cos(Omega2*t))
Nstepperiod=500
@njit
def DifSolving(Theta,I,N):
    List1=[]
    List2=[]
    t=0.
    dt=5.*2.*np.pi/(Omega1*Nstepperiod)
    for i in range(N):
        k11=dt*ThetaDif(t,Theta,I)
        k12=dt*IDif(t,Theta,I)
        k21=dt*ThetaDif(t+dt/2.,Theta+k11/2.,I+k12/2.)
        k22=dt*IDif(t+dt/2.,Theta+k11/2.,I+k12/2.)
        k31=dt*ThetaDif(t+dt/2.,Theta+k21/2.,I+k22/2.)
        k32=dt*IDif(t+dt/2.,Theta+k21/2.,I+k22/2.)
        k41=dt*ThetaDif(t,Theta+k31,I+k32)
        k42=dt*IDif(t,Theta+k31,I+k32)
        Theta+=(1./6.)*(k11+2.*k21+2.*k31+k41)
        I+=(1./6.)*(k12+2.*k22+2.*k32+k42)
        if Theta>np.pi:
            Theta=Theta-2.*np.pi
        if Theta<-np.pi:
            Theta=Theta+2.*np.pi
        if (i+1)%Nstepperiod==0:
            List1.append(Theta)
            List2.append(I)
        t+=dt
    return [List1,List2]
def PhasePortret(Th0,I0,step):
    # int((2.-I0)/step))
    for i in tqdm(range(int((1.-I0)/step))):
        Graph=DifSolving(Th0,I0+step*i,2000*Nstepperiod)
        ax.scatter(Graph[0], Graph[1], marker=".",s=0.6)
        ax.yaxis.set_label_coords(-1., 1.)
        ax.xaxis.set_label_coords(-np.pi, np.pi)
        ax.set_ylabel('f0 = '+str(f0)+'  k = '+str(k), rotation=0, fontsize=20, labelpad=20)
        ax.yaxis.set_label_coords(0.5, 1.01)
        ax.set_xticks(np.linspace(-3, 3, 21))
        ax.set_yticks(np.linspace(0.1, 1., 19))
fig, ax = plt.subplots(figsize=(10,40))
PhasePortret(0.,0.01, 0.005)
PhasePortret(np.pi,0.01, 0.008)
fig.savefig('DiffSolv.png', dpi=170)
