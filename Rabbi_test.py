import matplotlib.pyplot as plt
import numpy as np
from qutip import *

N = 2
wq = 1.0 #qubit frequency
w = 1.0 #frequency of the field exciting the qubit
Gf = 0.0 #phase relaxation
Ge = 0.0 #energy relaxation
Am=0.02 #qubit field amplitude
kappa = 0.0          # cavity dissipation rate
gamma = 0.0           # atom dissipation rate
Wrab = np.sqrt((wq-w)*(wq-w)+Am*Am)
Trubbi=2
#Trubbi=2*2*np.pi/np.sqrt((wq-w)*(wq-w)+Am*Am)

# T = Trubbi/2
# dt = 0.1
#Nstep = int(2 * np.pi /(w*dt))
tmax = Trubbi*2*2*np.pi/np.sqrt((wq-w)*(wq-w)+Am*Am)
tlist = np.linspace(0, tmax, 1000)
S1=Qobj([[1,0],[0,0]])
sm = destroy(2)
psi0 = basis(2,0)

def hamiltonian_t(t, args):
    H0 = args['H0']
    H1 = args['H1']
    #H2  = args['H2']
    w  = args['w']
    return H0 + H1*np.cos(w*t) #+H2

H0=0.5*wq*sigmaz()
H1=0.5*Am*sigmax()
#H2= - 0.5j*Gf*qeye(2)- 0.5j*Ge*S1

H_args = {'H0': H0, 'H1': H1,'w': w} #'H2': H2

# collapse operators
c_op_list = []
n_th = 0.0 # zero temperature

c_op_list = []

rate = kappa * (1 + n_th)
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * a)

rate = kappa * n_th
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * a.dag())

rate = gamma
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * sm)
        
output = mesolve(hamiltonian_t, psi0, tlist, c_op_list, sm.dag() * sm, H_args) 

func_analitic_rubi =Am*Am*np.sin(tlist*Wrab/4.)*np.sin(tlist*Wrab/4.)/(Wrab*Wrab)

fig, ax = plt.subplots(figsize=(10,5))
ax.plot(tlist, output.expect[0], label="Qubit Rubbi")
ax.plot(tlist, func_analitic_rubi, label="Qubit Rubbi analitic")
#ax.plot(tlist, output.expect[1], label="Atom excited state")
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Occupation probability')
ax.set_title('Vacuum Rabi oscillations')

