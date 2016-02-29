import numpy as np
from matplotlib import pyplot as plt
from scipy import special, optimize

"variables de gas 1"
rhol=1.0
vl=0.0
pl =1.0
"variables de gas 2"
rhor=0.125
vr=0.0
pr=0.1
"parametros relavantes"
t=1.0
gamma=1.4
m2=(gamma-1)/(gamma+1)
cleft=np.sqrt(gamma*pl/rhol)
"longitud del tubo"
L=4
"posicion de la membrana"
x0=L/2

N=5000

"funcionar a la que toca encontrarle el cero"
def P_post(p):
    return 2*np.sqrt(gamma)/(gamma-1)*(1-pow(p,(gamma-1)/(2*gamma)))-(p-pr)*(1-m2)/np.sqrt(rhor*(p+m2*pr))
p_post=optimize.newton(P_post,0)
def V_post(p_post):
    return 2*np.sqrt(gamma)/(gamma-1)*(1-pow(p_post,(gamma-1)/(gamma*2)))

def V_shock (rho_post):
    return v_post*(rho_post/rhor)/((rho_post/rhor)-1)
def Rho_post(p_post):
    return(rhor)*((p_post/pr)+m2)/(1+m2*(p_post/pr))

p_post=optimize.newton(P_post,0)
v_post=V_post(p_post)
rho_post=Rho_post(p_post)
v_shock=V_shock(rho_post)
r_mitad=rhol*pow((p_post/pl),1/gamma)


def c_sound(x):
    return m2*(x0-x)/t + ( 1 - m2)*cleft
def u(x):
    return ( 1 - m2 )*(- (x0-x)/t + cleft) ;
def r(x):
    return rhol*pow(( c_sound(x)/cleft ),2/(gamma-1))
def p(x):
    return pl*pow(( r(x)/rhol ),gamma)


x1=x0-cleft*t
x3=x0+v_post*t
x4=x0+v_post*t
x2 =t*v_post/(1-m2)-t*cleft+x0;
x=np.linspace(0,L,N)
e=np.linspace(0,L,N)
dominio=np.zeros((3,N))

for i in range (0, len(x)):
    a=x[i]
    if((a<x1)):
        dominio[0][i]=pl
        dominio[1][i]=rhol
        dominio[2][i]=vl
    
    elif ((a<x2)&(a>=x1)):
        dominio[0][i]=p(a)
        dominio[1][i]=r(a)
        dominio[2][i]=u(a)
    
    elif ((a<x3)&(a>=x2)):
        dominio[0][i]=p_post
        dominio[1][i]=r_mitad
        dominio[2][i]=v_post
    
    elif ((a<x4)&(a>=x3)):
        dominio[0][i]=p_post
        dominio[1][i]=rho_post
        dominio[2][i]=v_post
    
    elif ((a>=x4)):
        dominio[0][i]=pr
        dominio[1][i]=rhor
        dominio[2][i]=vr

for i in range (0,len(x)):
    e[i]=dominio[0][i]/((gamma-1)*dominio[1][i])



fig = plt.figure()

axis1 = fig.add_subplot(221)
axis1.set_title("presion")
axis1.plot(x,dominio[0],'r-')
axis1.set_ylim(0,1.5)

axis2 = fig.add_subplot(222)
axis2.set_title("densidad")
axis2.plot(x,dominio[1],'b-')
axis2.set_ylim(0,1.5)

axis3 = fig.add_subplot(223)
axis3.set_title("velocidad")
axis3.plot(x,dominio[2],'g-')
axis3.set_ylim(0,2.5)

axis4 = fig.add_subplot(224)
axis4.set_title("energia")
axis4.plot(x,e,'r--')





fig.savefig('comparacion.pdf')