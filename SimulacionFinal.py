#Codigo hecho por: Luis Cruz
import numpy as np
import matplotlib.pyplot as plt

#Constantes cambiadas a SI

q = 1.6*10**(-19) #Coulomb positron
m = 9.1*10**(-31) #Kg
E = 10500 #N/C (Efecto Stark)
B = 0.1 #Tesla (Efecto Zeeman)
c = 1 #(Para deshacerse de las unidades gaussianas)
omega_0 = 2*np.pi*10**(13) #Hz (Vibracion de un atomo)

#Valores iniciales
x_0 = 0 #m
y_0 = 0 #m
z_0 = 0 #m
v_0x = 0 #m/s
v_0y = 2200*10**3 #m/s #Entre mayor sea, puede crear un overflow
v_0z = 2200*10**3 #m/s


#Nuevas constantes definidas
omega_c = q*B/(m*c)
R_1 = -v_0y/(np.sqrt(omega_c**2 + 4*omega_0**2)) 
+ q*E*(1+omega_c/(np.sqrt(omega_c**2 + 4*omega_0**2)))/(2*m*omega_0**2)
-q*E/(m*omega_0**2)
R_2 = v_0y/(np.sqrt(omega_c**2 + 4*omega_0**2)) 
- q*E*(1+omega_c/(np.sqrt(omega_c**2 + 4*omega_0**2)))/(2*m*omega_0**2)

#Rango de la simulacion y pasos
a = 0
b = 6*10**(-13)
n= 500
t = np.linspace(a, b,n+1)
#Calculo analitico

def x_a(t):
    return R_1*np.cos(0.5*(np.sqrt(omega_c**2 + 4*omega_0**2)+omega_c)*t) + \
    R_2*np.cos(0.5*(np.sqrt(omega_c**2 + \
    4*omega_0**2)-omega_c)*t) + q*E/(m*omega_0**2)
def y_a(t):
    return -R_1*np.sin(0.5*(np.sqrt(omega_c**2 + \
    4*omega_0**2)+omega_c)*t) + \
    R_2*np.sin(0.5*(np.sqrt(omega_c**2 + 4*omega_0**2)-omega_c)*t)

def z_a(t):
    return v_0z*np.sin(omega_0*t)/omega_0

#Listas de posiciones de resultados analiticos

def Real_x(a,b,n,x_0): #Solucion analitica
  h=(b-a)/n
  T = []
  T.append(a)
  Y = []
  Y.append(x_0)
  
  for i in range(1,n+1):
    T.append(T[i-1]+h)
    Y.append(x_a(T[i]))
    
  return Y

def Real_y(a,b,n,y_0): #Solucion analitica
  h=(b-a)/n
  T = []
  T.append(a)
  Y = []
  Y.append(y_0)
  
  for i in range(1,n+1):
    T.append(T[i-1]+h)
    Y.append(y_a(T[i]))
    
  return Y

def Real_z(a,b,n,z_0): #Solucion analitica
  h=(b-a)/n
  T = []
  T.append(a)
  Y = []
  Y.append(z_0)
  
  for i in range(1,n+1):
    T.append(T[i-1]+h)
    Y.append(z_a(T[i]))
    
  return Y

#Calculo numerico
def V(u,t):
  mu,x,nu,y = u
  return np.array([omega_c*nu-x*omega_0**2 +q*E/m,mu,-omega_c*mu-y*omega_0**2,nu])

def Numerico(f,u0,a,b,n):
    t = np.linspace(a, b, n+1)
    u = np.array((n+1)*[u0],dtype='f')
    h = t[1]-t[0]
    for i in range(n):
        k1 = h * f(u[i], t[i])    
        k2 = h * f(u[i] + 0.5 * k1, t[i] + 0.5*h)
        k3 = h * f(u[i] + 0.5 * k2, t[i] + 0.5*h)
        k4 = h * f(u[i] + k3, t[i] + h)
        u[i+1] = u[i] + (k1 + 2*(k2 + k3 ) + k4) / 6
    return u, t

u, t  = Numerico(V, np.array([v_0x, x_0, v_0y , y_0]) , a , b , n)
mu,x,nu,y  = u.T

#Creacion de un.txt
np.savetxt('Congreso.txt',np.c_[x,y,Real_z(a,b,n,z_0)],fmt='%1.15f')

#Creacion de graficas

plt.plot(x_a(t), y_a(t), color='purple') #Analitica x-y

#Parametrizacion1
ax = plt.figure().add_subplot(projection='3d')
ax.plot(x_a(t), y_a(t), z_a(t), label='parametric curve')
ax.legend()

#Parametrizacion2
bx = plt.figure().add_subplot(projection='3d')
bx.plot(x_a(t), y_a(t), label='parametric curve x-y')
bx.legend()

#Parametrizacion3
cx = plt.figure().add_subplot(projection='3d')
cx.plot(x,y,Real_z(a,b,n,z_0),label='parametric curve 2')
cx.legend()

plt.rcParams["figure.figsize"] = [10, 10]
plt.rcParams["figure.autolayout"] = True
plt.show()