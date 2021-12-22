import matplotlib.pyplot as plt
import numpy as np
from qutip import *

N = 2
wq = 1.0 #qubit frequency
w = 1.01 #frequency of the field exciting the qubit
Gf = 0.0 #phase relaxation
Ge = 0.03 #energy relaxation
Am=0.6 #qubit field amplitude

T = 2 
dt = 0.1
Nstep = int(2 * np.pi /dt) 
tmax = Nstep*T
tlist = np.linspace(0, tmax, Nstep)
S1=Qobj([[1,0],[0,0]])
sm = destroy(2)
psi0 = basis(2,0)

def hamiltonian_t(t, args):
    H0 = args['H0']
    H1 = args['H1']
    return H0+H1*np.cos(w*t)+H2

H0=0.5*wq*sigmaz()
H1=0.5*Am*sigmax()
H2= - 0.5j*Gf*qeye(2)- 0.5j*Ge*S1

H_args = {'H0': H0, 'H1': H1,'H1': H2}

# collapse operators
c_op_list = []
n_th = 0.5 # zero temperature

# relaxation
rate = Ge * (1 + n_th)
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * sm)

# excitation
rate = Ge * n_th
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * sm.dag())

# dephasing 
rate = Gf
if rate > 0.0:
    c_op_list.append(np.sqrt(rate) * sz)
        
output = mesolve(hamiltonian_t, psi0, tlist, c_op_list, sm.dag() * sm, H_args) 

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(tlist, output.expect[0], label="Cavity")
#ax.plot(tlist, output.expect[1], label="Atom excited state")
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Occupation probability')
ax.set_title('Vacuum Rabi oscillations');


