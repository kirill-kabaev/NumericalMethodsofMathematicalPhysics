import matplotlib.pyplot as plt
import numpy as np
import time
n=400
p=300
r=150

w=[]
mas=[]
time11=time.time() 
def l(x):
    return x*x-0.4+0.6j
re=np.arange(-2., 2., float(4)/float(n))
im=np.arange(-2., 2., float(4)/float(p))

graf=[]
z=[]
j=0
q=0
for b in range(n*p):
    if q<(n):        
        z.append(re[q]+im[j]*1j)         
    else: 
        j=j+1
        q=0
        z.append(re[q]+im[j]*1j)
    q=q+1
    k=(l(z[b]))
    i=1
    while abs(k)<=2 and i<=r:
        i=i+1
        f=l(k)
        k=f
    mas.append(i)
    for v in range(p+1):
        v=v+1
        if (b==(v*n-1)):
            graf.append(mas)
            mas=[]
time12=time.time() 
plt.xlabel('ReZ')
plt.ylabel('ImZ')
print('nashe vremya=', (time12-time11))
im=plt.imshow(graf,origin='lower', cmap='hot', interpolation='none', extent=[-2,2,-2,2])
#im=plt.imshow(graf,origin='lower', cmap='gnuplot_r', interpolation='none', extent=[-2,2,-2,2])
#im=plt.imshow(graf,origin='lower', cmap='gist_rainbow', interpolation='none', extent=[-2,2,-2,2])
#im=plt.imshow(graf,origin='lower', cmap='RdGy', interpolation='none', extent=[-2,2,-2,2])
#im=plt.imshow(graf,origin='lower', cmap='flag', interpolation='none', extent=[-2,2,-2,2])
plt.show()
