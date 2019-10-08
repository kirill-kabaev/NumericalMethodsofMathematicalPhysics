from mpl_toolkits.mplot3d import Axes3D  

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
from ipywidgets import widgets
from IPython.display import display

Grafic = "filename.png"



def update_plot(p, q, r, s, Thetha, Phi):
        N = 64
        M = 32
        lx=1
        ly=1
        c=0.5
        g=[]
        f=[]
        h=[]
        l=[]
        pStr=str(p)
        qStr=str(q)
        rStr=str(r)
        sStr=str(s)
        fig = plt.figure(figsize=(12, 10))
        ax = fig.gca(projection='3d')
        #ax.set_title('p = '+pStr+'    '+'q = '+qStr+'    '+'r = '+rStr+'    '+'s = '+sStr)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        x = np.linspace(0, lx, N)
        y = np.linspace(0, ly, M)

        for i in range(N):
            f.append([])
            for j in range(M):
                f[i].append(c*(np.sin(math.pi*(x[i]/lx) ** q)** p) * (np.sin(math.pi*(y[j]/ly) ** r) ** s)) #x - i -N, y -j- M
        f = np.asarray(f, dtype=float)

        def IDFT(f):
            f = np.asarray(f, dtype=float)
            M = f.shape[0]
            m = np.arange(M)
            wy = m.reshape((M, 1))
            H = (1/M)*np.exp(2j * np.pi * (wy*m) / M)
            return np.dot(H, f)
        def IFFT(f):
            f = np.asarray(f, dtype=float)
            M = f.shape[0]                 #  индексы k
            if (M % 2) > 0:
                raise ValueError("size of x must be a power of 2")
            elif M <= 32:  
                return IDFT(f)
            else:
                X_even = IFFT(f[::2])
                X_odd = IFFT(f[1::2])        
                factor =  (1/M)*np.exp(2j * np.pi * np.arange(M)/M)
                return np.concatenate([X_even + factor[:int(M / 2)] * X_odd,
                                    X_even + factor[int(M / 2):] * X_odd])


        furie = []
        for r in range(N): 
            h=IFFT(f[r])
            furie.append([]) 
            for c in range(M): 
                furie[r].append(h[c])
        furie = np.asarray(furie, dtype=float)
        #метод прогонки
        Lambda=[]
        d = ly / M
        for k in range(M):
            Lambda = np.append( Lambda, ((math.pi*k)/ly)**2)
        Lambda = np.asarray(Lambda, dtype=float)

        e = [0]*M
        e = np.asarray(e, dtype=float)
        for i in range(M-1):
            e[i+1] = 1/(2 + Lambda[i+1]*(d**2) - e[i])
            Ck=[[0]*M]*N
            fi=[0]*M
        for j in range(N-1): 
            fi = np.asarray(fi, dtype=float)
            Ck = np.asarray(Ck, dtype=float)
            n = M
            d = ly / M
            for i in range(M-1):
                fi[i+1] = e[i+1]*(fi[i]- furie[j][i+1]*(d**2))
            for i in reversed(range(M-1)):
                Ck[j][i] = e[i]*Ck[j][i+1] + fi[i]


        def DFT(f):
            f = np.asarray(f, dtype=float)
            M = f.shape[0]
            m = np.arange(M)
            wy = m.reshape((M, 1))
            H = np.exp(-2j * np.pi * (wy*m) / M)

            #return(np.size(M))
            return np.dot(H, f)
        def FFT(f):
            f = np.asarray(f, dtype=float)
            M = f.shape[0]                 #  индексы k
            if (M % 2) > 0:
                raise ValueError("size of x must be a power of 2")
            elif M <= 32:  
                return DFT(f)
            else:
                X_even = FFT(f[::2])
                X_odd = FFT(f[1::2])        
                factor =  np.exp(-2j * np.pi * np.arange(M)/M)
                return np.concatenate([X_even + factor[:int(M / 2)] * X_odd,
                                               X_even + factor[int(M / 2):] * X_odd])
        #Востонавливаем решение
        U=[]
        s=[]
        for i in range(N):
            s= FFT(Ck[i])
            U.append([])
            for j in range(M):
                U[i].append(s[j])    
        U = np.asarray(U, dtype=float)

        #Построение графика
        x, y = np.meshgrid(x, y)
        f = f.reshape((M, N))
        furie = furie.reshape((M, N))
        Ck = Ck.reshape((M, N))
        U = U.reshape((M, N))

        ax.view_init(Thetha, Phi)  #поворот графика
        
        surf = ax.plot_surface(x, y, U.real, cmap=cm.coolwarm,   #если .real поменять на .imag, то построит мнимую часть, карты: cm.coolwarm,'viridis'
                                      linewidth=0)                  
        plt.show()
        fig.savefig(Grafic, dpi=170)

   
p_height = widgets.FloatSlider(min=0.5, max=2, step=0.1, description= 'p')
q_height = widgets.FloatSlider(min=0.5, max=2, step=0.1, description= 'q')
r_height = widgets.FloatSlider(min=0.5, max=2, step=0.1, description= 'r')
s_height = widgets.FloatSlider(min=0.5, max=2, step=0.1, description= 's')
Thetha_height = widgets.IntSlider(min=0, max=360, value =30, step=5, description= 'Thetha')
Phi_height = widgets.IntSlider(min=0, max=360, value = 60, step=5, description= 'Phi')

widgets.interactive(update_plot, p=p_height, q=q_height, r=r_height, s=s_height, Thetha=Thetha_height, Phi=Phi_height)       
   
