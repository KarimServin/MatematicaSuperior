
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from sympy.abc import x, a, b, i, j, n;


def Runge_Kutta_2(x0, t0, h, v):
    t1 = t0+h
    k1 = F(t0, x0, v)
    print(t0)
    k2 = F(t1, x0 + k1*h, v)
    v1 = x0 + h/2*(k2+k1)

    e_t = (np.sinh(t1)-v1[0])/np.sinh(t1)

    return v1, t1,e_t


def Runge_Kutta_4(x0, t0, h, v):
    t1 = t0 + h/2
    k1 = F(t0, x0, v)
    print(t0)
    k2 = F(t1, x0 + k1*h/2, v)
    k3 = F(t1, x0 + k2*h/2, v)
    k4 = F(t0, x0 + k3*h, v)
    v1 = x0 + h/6*(k1 + 2*k2 + 2*k3 + k4)

    e_t = (sinh(t1)-v1[0])/sinh(t1)

    return v1, t1, e_t



dx = 0.01


tiempo=[0.00,0.40,0.80,1.00,1.15,1.30,1.50,1.70,1.90,2.00,2.10,2.30,2.40,2.50,2.60,2.70,3.00,3.30,3.60,4.00,4.50,5.00,5.50,6.00]
voltaje=[-70.00,-70.00,-69.72,-65.78,-56.94,-48.28,-34.49,-15.21, 10.96, 29.44, 39.64, 14.19,-16.24,-45.10,-65.76,-78.98,-87.38,-84.70,-80.08,-75.12,-71.00,-70.00,-70.00,-70.00]


285.0
285.0
300.0
600.0
285.0
600.0
600.0
600.0
275.0
275.0
550.0
550.0
550.0
1100.0
1100.0
385.0
275.0
275.0
275.0
275.0
275.0
275.0
550.0
550.0
550.0
275.0
550.0
275.0
500.0
500.0
500.0
200.0
200.0
200.0
450.0
450.0
450.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
220.0
450.0
450.0
220.0
220.0
220.0
220.0
220.0
200.0
200.0
285.0
285.0
300.0
600.0
285.0
600.0
600.0
600.0
275.0
275.0
550.0
550.0
550.0
1100.0
1100.0
385.0
275.0
275.0
275.0
275.0
275.0
275.0
550.0
550.0
550.0
275.0
550.0
275.0
500.0
500.0
500.0
200.0
200.0
200.0
450.0
450.0
450.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
45.0
220.0
450.0
450.0
220.0
220.0
220.0
220.0
220.0
200.0
200.0










def pol3natural(xi, yi):
    # ASplines 3
    n = len(xi)

    # h funciona como array Calculo las diferencias en X entre cada punto
    h = np.zeros(n - 1, dtype=float) #crea array de ceros de tamano n-1
    for j in range(0, n - 1, 1):
        h[j] = xi[j + 1] - xi[j]

    '''
        A matriz con los u_i y los h_i
        B el vector con los  v_i
        Z son los valores z_i
    '''
    A = np.zeros(shape=(n - 2, n - 2), dtype=float)
    B = np.zeros(n - 2, dtype=float)
    Z = np.zeros(n, dtype=float)

    # ARMO EL SISTEMA
    # Estos valores van fuera del for ya que son las unicas que no puedo armar con un ciclo
    A[0, 0] = 2 * (h[0] + h[1])
    A[0, 1] = h[1]
    B[0] = 6 * ((yi[2] - yi[1]) / h[1] - (yi[1] - yi[0]) / h[0])  # <---- V_0
    A[n - 3, n - 4] = h[n - 3]
    A[n - 3, n - 3] = 2 * (h[n - 3] + h[n - 2])
    B[n - 3] = 6 * ((yi[n - 1] - yi[n - 2]) / h[n - 2] -
                    (yi[n - 2] - yi[n - 3]) / h[n - 3])

    for i in range(1, n - 3, 1):
        # Vease la matriz creada en el informe
        A[i, i - 1] = h[i]
        A[i, i] = 2 * (h[i] + h[i + 1])  # <-- U_i
        A[i, i + 1] = h[i + 1]
        B[i] = 6 * ((yi[i + 2] - yi[i + 1]) / h[i + 1] -
                    (yi[i + 1] - yi[i]) / h[i])  # <---- V_i

    # INGRESO LOS VALORES DE Z como deben ir de acuerdso a que es natural
    r = np.linalg.solve(A, B)
    for j in range(1, n - 1, 1):
        Z[j] = r[j - 1]
    Z[0] = 0
    Z[n - 1] = 0

    # Arreglos de los coeficientes para cada polinomio.
    a = np.zeros(n - 1, dtype=float)
    b = np.zeros(n - 1, dtype=float)
    c = np.zeros(n - 1, dtype=float)
    d = np.zeros(n - 1, dtype=float)
    for j in range(0, n - 1, 1):
        a[j] = (Z[j + 1] - Z[j]) / (6 * h[j])
        b[j] = Z[j] / 2
        c[j] = (yi[j + 1] - yi[j]) / h[j] - \
            (2 * h[j] * Z[j] + h[j] * Z[j + 1]) / 6
        d[j] = yi[j]

    # GUardado
    px_tabla = []
    for j in range(0, n - 1, 1):
        # Construyo elpolinomio
        new_pol = a[j] * (x - xi[j]) ** 3 + b[j] * \
            (x - xi[j]) ** 2 + c[j] * (x - xi[j]) + d[j]
        new_pol = new_pol.expand()
        px_tabla.append(new_pol)

    return px_tabla



print("-------------------- Calculos de LTI --------------------")

dominio = np.array(tiempo)
volt = np.array(voltaje)

for j in range(1, 2):
    corrimiento = j * 10*np.ones_like(dominio)
    volt = np.concatenate((volt, volt))
    newconcate = np.concatenate((dominio, dominio+corrimiento))
    dominio = newconcate
    lineabase = 70 * np.ones_like(volt)
    volt += lineabase


def F(t, x, v):
    result = np.array([1, 1])
    result[0] = float(x[1])
    result[1] = float(v[int(t/dx)]-3*x[1]-2*x[0]) #FALTA EL 4.9*
    return result


def milneCorrector(Yi, y, t, h,yd):
    yf =Yi[int(t/h)]+h*(1*y[int(t/h)]+81/100*y[int((t+h/2)/h)]+ 4*y[int((t+3/2*h)/h)]+1*y[int((t+5/2*h)/h)]+(7/50)*yd)
    return yf, t+3*h

def milnePredictor( y1, y2d, t, h):
    ysig=y1[int(t/h)]+h*(2*y2d[int((t+1*h)/h)]+1*y2d[int((t+2*h)/h)]+2*y2d[int((t+3*h)/h)])
    return ysig,t+h



n = len(dominio)
polinomios = pol3natural(dominio, volt)


print('Polinomios: ')
for tramo in range(1, n, 1):
    print('Polinomio: '+str(tramo)+' x = [' + str(dominio[tramo - 1])
          + ',' + str(dominio[tramo]) + ']')
    print(str(polinomios[tramo - 1]))

x_v = np.array([])
y_v = np.array([])
y_deriv = np.array([])
for i in range(0, n):
    # Limitres del tramo
    a = dominio[i - 1]
    b = dominio[i]
    x_aux = np.arange(a, b, dx)
    # Polinomio y su derivada
    pol_x = polinomios[i - 1]

    pol = lambdify('x', pol_x)
    
    y_aux = pol(x_aux)
    # Concateno
    x_v = np.concatenate((x_v, x_aux))
    y_v = np.concatenate((y_v, y_aux))



h = dx

vp = [0.0]
vp_p = [0.0]
tp = [0]
vn = [0, 0]
t0 = 0

vn2 = [0, 0]
t02 = 0
v2 = [0.0]
v2_p = [0]
v1_p = [0]
t2 = [0]
for j in range(0, int(10/dx)):
    vn, t1, e12 = Runge_Kutta_2(vn, t0, h, y_v)
    t0 = t1
    vp.append((y_v[int(t1/h)]+3*vn[1]+2*vn[0])/2) #Falta el 4.9*
    tp.append(t1)
    
    vn2, tn2, e_14 = Runge_Kutta_4(vn2, t02, h, y_v)
    t02 = tn2
    v2.append((y_v[int(t1/h)]+3*vn2[1]+2*vn2[0])/2) #Falta el 4.9*
    t2.append(tn2)

    vn2, tn2, e_14 = Runge_Kutta_4(vn2, t02, h, y_v)
    t02 = tn2
    v2.append((y_v[int(t1/h)]+3*vn2[1]+2*vn2[0])/2) #Falta el 4.9*
    v2_p.append(vn2[1])
    v1_p.append(vn2[0])
    t2.append(tn2)


# aplicar proyector corrector
# Condiciones iniciales


y_pi = vn2[0]
y_ppi = vn2[1]
y_y = [0]
tr =[0]
t0 = 0

print(vn2[0])
print('hola')
print(vn2[1])


for j in range(0,int(3/dx)):
    y_pp1, ts = milnePredictor(v1_p,v2_p,t0, dx)
    y_c1 , t1 = milneCorrector(v1_p,v2_p,t0,dx,y_pp1)
    tr.append(t1)
    y_3h = (y_v[int(t1/h)]+3*y_pp1+2*y_c1)/2  #FALTA EL 4.9*
    t0=t1
    y_y.append(y_3h)



plt.xlabel('t [ms]')
plt.ylabel('Potencial de accion [mV]')
plt.plot(dominio,volt, label='original sin interpolar')
plt.plot(x_v, y_v, label='original interpolada Spline3')
plt.plot(tr, y_y,label='predictor-corrector ')
plt.plot(tp,vp,label='rungeKutta 2')
plt.plot(t2, v2,label='rungeKutta 4')



plt.legend()
plt.xlim(0,7.5)
plt.show()