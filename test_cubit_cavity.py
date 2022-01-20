import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import numpy as np
from qutip import *
from scipy import linalg

N = 2
wq =0.25 #qubit frequency (delta)
t1 = 1
t2 = 8
# A = 0.25
# Gf = 0.0 #phase relaxation
# Ge = 0.0 #energy relaxation
# Am=0.02 #qubit field amplitude
# Wrab = np.sqrt((wq-w)*(wq-w)+Am*Am)
# Trubbi=2
#Trubbi=2*2*np.pi/np.sqrt((wq-w)*(wq-w)+Am*Am)

# T = Trubbi/2
# dt = 0.1
#Nstep = int(2 * np.pi /(w*dt))
# tmax = Trubbi*2*2*np.pi/np.sqrt((wq-w)*(wq-w)+Am*Am)

tmax = 10
#Nstep=int(tmax/(2*np.pi/w))*6
Ntime= 1000
tlist = np.linspace(0, tmax, Ntime)
S1=Qobj([[1,0],[0,0]])
sm = destroy(2)
sz = sigmaz()
psi0 = basis(2,0)


def hamiltonian_t(t, args):
        H0 = args['H0']
        H1 = args['H1']
        #H2  = args['H2']
        return H0 + H1*Theta(t) #+H2
def Theta(t):
    if ((t>t1) and (t<t2)):
        return 1
    else:
        return 0
Num_A = 250
ListPsiPsi=[]
listTheta=[]
Alist = np.linspace(0, 5, Num_A)
for j in range(Num_A):
        ListPsiPsi.append(Ntime*[0])
for i in range(Num_A):
    A = Alist[i]
    print("progress ", i ,"of", Num_A )
    H0 = - 0.5*wq*sigmaz()
    H1 = - 0.5*A*sigmax()
    H_args = {'H0': H0, 'H1': H1} #'H2': H2

    output = sesolve(hamiltonian_t, psi0, tlist, [], H_args, progress_bar=True) # sm.dag() * sm, H_args) 
    for j in range(Ntime):
        ListPsiPsi[i][j]=linalg.blas.cdotc(output.states[j][0],output.states[j][0]).real
        #print(ListPsiPsi)
        listTheta.append(Theta(tlist[j]))
    
    
fig, ax = plt.subplots(figsize=(10,20))
X,Y = np.meshgrid(tlist, Alist)
Z = np.asarray(ListPsiPsi, dtype=float)
cp = ax.contourf(X, Y, Z, cmap='jet')

ax.set_title('$|Ψ|^2$', fontsize=20, pad=20) #$|Sx|^2$
ax.set_ylabel('A', rotation=0, fontsize=20, labelpad=20)
ax.set_xlabel('t', rotation=0, fontsize=20, labelpad=20)
ax.tick_params(labelsize=15)
ax.yaxis.set_label_coords(0.01, 1.02)
ax.xaxis.set_label_coords(1.05, -0.05)
plt.show()


plt.savefig('saved_figure_test22_0.png', dpi=200)
