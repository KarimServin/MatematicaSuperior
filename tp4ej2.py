import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import scipy as sc


dx=0.1 

dominio=np.array([0.00,0.40,0.80,1.00,1.15,1.30,1.50,1.70,1.90,2.00,2.10,2.30,2.40,2.50,2.60,2.70,3.00,3.30,3.60,4.00,4.50,5.00,5.50,6.00])

volt=np.array([-70.00,-70.00,-69.72,-65.78,-56.94,-48.28,-34.49,-15.21, 10.96, 29.44, 39.64, 14.19,-16.24,-45.10,-65.76,-78.98,-87.38,-84.70,-80.08,-75.12,-71.00,-70.00,-70.00,-70.00])

for j in range(1,5):
    corrimiento= j*10*np.ones_like(dominio)
    volt=np.concatenate((volt,volt))
    newconcate=np.concatenate((dominio,dominio+corrimiento))
    dominio=newconcate
    lineabase=70*np.ones_like(volt)
    volt=+lineabase


def F(t,x,v):
    result= np.array([1,1])
    result[0]= float(x[1])
    result[1]=float(4.9*v[int(t/dx)]-3*x[1]-2*x[0])
    return result


def Runge_Kutta_2(x0,t0,h,v):
    t1=t0+h
    k1=F(t0,x0,v)
    print(t0)
    k2=F(t1,x0+k1*h,v)
    v1=x0+h/2*(k2+k1)
    e_t = (np.sinh(t1)-v1[0])/np.sinh(t1)

    return v1,t1,e_t


def Runge_Kutta_4(x0,t0,h,v):
    t1=t0+h/2
    k1=F(t0,x0,v)
    print(t0)
    k2=F(t1,x0+k1*h/2,v)
    k3=F(t1,x0+k2*h/2,v)
    k4=F(t0,x0+k3*h,v)
    v1=x0+h/6*(k1+2*k2+2*k3+k4)

    e_t = (np.sinh(t1)-v1[0])/np.sinh(t1)
    return v1,t1,e_t


def Predictor_Corrector(Y,y,y_prima,t,h):
    Y_p=y[int(t/dx)]+3*h/20*y_prima[int((t+5*h/6)/dx)]
    Y_c=Y[int(t/dx)]+[]




x_v = np.array([])
y_v= np.array([])
y_deriv=np.array([])


h=dx;

vp=[0.0]
tp=[0]
vn=[0,0]
e_2=[0]
t0=0

vn2=[0,0]
t02=0
v2=[0.0]
t2=[0]
e_4=[0]

for j in range(0,int(10/dx)):
    
    
    vn,t1,e12=Runge_Kutta_2(vn,t0,h,y_v)


    t0 = t1
    vp.append(vn[0])
    tp.append(t1)
    e_2.append(e12*100)
    vn2,tn2,e_14= Runge_Kutta_4(vn2,t02,h,y_v)
    t02=tn2
    v2.append(vn[0])
    t2.append(tn2)
    e_4.append(e_14*100)
    vn2,tn2,e_14= Runge_Kutta_4 (vn2,t02,h,y_v)
    t02=tn2
    v2.append(vn2[0])
    t2.append(tn2)
    e_4.append(e_14*100)





plt.plot(x_v,y_v,label='calamar')
plt.plot(tp,vp,'-.',label='R-K 2')
plt.plot(t2,v2,'--',label='R-K 4')
plt.legend()
plt.xlim(0,10)
plt.grid()
plt.show()