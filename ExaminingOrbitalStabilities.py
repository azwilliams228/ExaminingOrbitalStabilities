from numpy import *
import matplotlib.pyplot as plt
%matplotlib notebook
%matplotlib inline
from numpy import linalg as LA
def f(X,t): 
    x,y = X
    return array([-y,x])

def eulerstep(f,X,t,h):        # arguments: function, variable, time, step size
    return X + f(X,t)*h

def improved_euler_step(f,X,t,h):
    f1 = f(X,t)
    X1 = X + f1*h              #provisional Euler Step
    f2 = f(X1,t+h)
    fav = (f1+f2)/2
    return X + h*fav
    
X = array([1.,0.])             # putting . after the number changes it to float
T = 2*pi                       # final time we want to get to
n = 50                         # number of steps
h = T/n                        # time-step ("delta-t")

print(0,X)
plt.subplot(111,aspect=1)

# draw exact solution
t = linspace(0.1*pi,300)
plt.plot(cos(t),sin(t),'k', alpha = 0.5)
for i in range(n):
    newX = improved_euler_step(f,X,t,h)     # Euler step
    plt.plot([X[0],newX[0]],[X[1],newX[1]],'bo-',ms =  2, alpha = 0.3)
    X = newX
    t = (i+1)*h
    #print(t,X,newX)
plt.show()
from IPython.display import Image
Image('improvedeulererror.png')
def f(X,t): 
    x1 = R1 * cos(t)
    y1 = R1 * sin(t)
    x2 = -R2 * cos(t)
    y2 = -R2 * sin(t)
    x,y,u,v = X
    r1 = sqrt(abs(x-x1)**2 + abs(y-y1)**2)
    r2 = sqrt(abs(x-x2)**2 + abs(y-y2)**2)
    r132 = r1**3
    r232 = r2**3
    A = array([u,v,-G*(x-x1)/r132-G*(x-x2)/r232, -G*(y-y1)/r132-G*(y-y2)/r232])            # there are 2 suns
    return A
def rk4(f,X,t,h):
    k1 = h*f(X,t)
    k2 = h*f(X+k1/2,t+h/2)
    k3 = h*f(X+k2/2,t+h/2)
    k4 = h*f(X+k3, t+h)
    return X + 1/6*(k1+2*k2+2*k3+k4)
G = 1                      # The Gravitational constant, 1 for simplicity
M1 = 0.3                   # The mass of Sun # 1
M2 = 1-M1                  # The mass of Sun #2
R2 = M1                    # Ratio of mass is inverse the ratio of the radii of the center of gravity of the two suns
R1 = M2
M = M1 + M2
X = array([1.,0.])                      # putting . after the number changes it to float
T = 5000                                # final time we want to get to                           
c = 0.05

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)

fig = plt.figure(figsize = (4,4));
ax = plt.subplot(111,aspect = 1);
ax.plot(x1,y1,'b');
ax.plot(x2,y2,'r');

for u0 in [1.3]:
    X = array([0,2,u0,0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] )          # step size
        newX = rk4(f,X,t,h)                  # Euler step
        plt.plot([X[0],newX[0]],[X[1],newX[1]],'c-',ms =  2, alpha = 0.3);
        X = newX
        t += h
plt.show();

T = 100000 

fig = plt.figure(figsize = (4,4));
ax = plt.subplot(111,aspect = 1);
ax.plot(x1,y1,'b');
ax.plot(x2,y2,'r');

for u0 in [1.3]:
    X = array([0,2,u0,0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] )          # step size
        newX = rk4(f,X,t,h)                  # Euler step
        plt.plot([X[0],newX[0]],[X[1],newX[1]],'c-',ms =  2, alpha = 0.3);
        X = newX
        t += h
plt.show();
T = 2000
c = 0.025
M = M1 + M2
colors = 'gkcmr'
y0s = [5.0]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))
ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors): 
    X = array([0,y0,sqrt(G*M/y0),0])
    t = 0
    while t < T:
       
        h = c/LA.norm( f(X,t)[2:] ) 
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h

plt.show();
T = 10000
c = 0.05

colors = 'gkcmr'
y0s = [5.0]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))
ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors): 
    X = array([0,y0,sqrt(G*M/y0),0])
    t = 0
    while t < T:
       
        h = c/LA.norm( f(X,t)[2:] ) 
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h

plt.show();
T = 10000
c = 0.05
M = M1 + M2
colors = 'gkcmr'
y0s = [10.0,25.0,50.0]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))
ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors): 
    X = array([0,y0,sqrt(G*M/y0),0])
    t = 0
    while t < T:
       
        h = c/LA.norm( f(X,t)[2:] ) 
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h

plt.show();
T = 10000                         
c = 0.05

colors = 'gkcmr'
u0s = [sqrt(G*M/y0)]
y0s = [10.0, 7.5, 6.5]

for y0 in y0s:
    t = linspace(0,2*pi, 100, endpoint = False)
    x1 = R1 * cos(t)
    y1 = R1 * sin(t)
    x2 = -R2 * cos(t)
    y2 = -R2 * sin(t)
    fig = plt.figure(figsize = (4,4))
    ax = plt.subplot(111,aspect = 1)
    ax.plot(x1,y1,'b')
    ax.plot(x2,y2,'r');
    
    for u0,color in zip(u0s, colors): 
        X = array([0,y0,sqrt(G*M/y0),0])
        t = 0
        while t < T:

            h = c/LA.norm( f(X,t)[2:] ) 
            newX = rk4(f,X,t,h)
            plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
            X = newX
            t += h

plt.show();
T = 10000
c = 0.05

y0s = [1.0,3.0,5.0,6.0,7.0,10.0]
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(18, 18)) 
gs = gridspec.GridSpec(2, 3)

ax1 = plt.subplot(gs[0,0])
ax1.title.set_text('Plot 1, Inititial position at y = 1.0')
ax2 = plt.subplot(gs[1,0])
ax2.title.set_text('Plot 4, Inititial position at y = 6.0')
ax3 = plt.subplot(gs[0,1])
ax3.title.set_text('Plot 2, Inititial position at y = 3.0')
ax4 = plt.subplot(gs[1,1])
ax4.title.set_text('Plot 5, Inititial position at y = 7.0')
ax5 = plt.subplot(gs[0,2])
ax5.title.set_text('Plot 3, Inititial position at y = 5.0')
ax6 = plt.subplot(gs[1,2])
ax6.title.set_text('Plot 6, Inititial position at y = 10.0')

for y0 in y0s:
    t = linspace(0,2*pi, 100, endpoint = False)
    x1 = R1 * cos(t)
    y1 = R1 * sin(t)
    x2 = -R2 * cos(t)
    y2 = -R2 * sin(t)
    #fig = plt.figure(figsize = (4,4))
    #ax = plt.subplot(111,aspect = 1)
    ax1.plot(x1,y1,'b')
    ax1.plot(x2,y2,'r');

    ax2.plot(x1,y1,'b')
    ax2.plot(x2,y2,'r');
    ax3.plot(x1,y1,'b')
    ax3.plot(x2,y2,'r');
    ax4.plot(x1,y1,'b')
    ax4.plot(x2,y2,'r');
    ax5.plot(x1,y1,'b')
    ax5.plot(x2,y2,'r');
    ax6.plot(x1,y1,'b')
    ax6.plot(x2,y2,'r');
    
    X = array([0,y0,sqrt(G*M/y0),0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] ) 
        newX = rk4(f,X,t,h)
        if y0 == y0s[0]:        
            ax1.plot([X[0],newX[0]],[X[1],newX[1]], 'g',ms =  2, alpha = 0.3);    
            ax1.set_xlim(-5,5)
            ax1.set_ylim(-5,5)
        if y0 == y0s[1]:        
            ax3.plot([X[0],newX[0]],[X[1],newX[1]], 'g',ms =  2, alpha = 0.3);
            ax2.set_xlim(-10,10)
            ax2.set_ylim(-10,10)
        if y0 == y0s[2]:        
            ax5.plot([X[0],newX[0]],[X[1],newX[1]], 'g',ms =  2, alpha = 0.3);
            ax3.set_xlim(-10,10)
            ax3.set_ylim(-10,10)
        if y0 == y0s[3]:        
            ax2.plot([X[0],newX[0]],[X[1],newX[1]], 'g',ms =  2, alpha = 0.3);
            ax4.set_xlim(-10,10)
            ax4.set_ylim(-10,10)
        if y0 == y0s[4]:        
            ax4.plot([X[0],newX[0]],[X[1],newX[1]], 'g',ms =  2, alpha = 0.3);
            ax5.set_xlim(-10,10)
            ax5.set_ylim(-10,10)
        if y0 == y0s[5]:        
            ax6.plot([X[0],newX[0]],[X[1],newX[1]], 'g',ms =  2, alpha = 0.3);
            ax6.set_xlim(-12,12)
            ax6.set_ylim(-12,12)
        X = newX
        t += h

plt.show();
T = 1                                # final time we want to get to                           
c = 0.1

colors = 'gkcmr'
u0s = [0.2,0.5,1.]
y0s = [0.05]

for u0 in u0s:
    t = linspace(0,2*pi, 100, endpoint = False)
    x1 = R1 * cos(t)
    y1 = R1 * sin(t)
    x2 = -R2 * cos(t)
    y2 = -R2 * sin(t)
    fig = plt.figure(figsize = (4,4))
    ax = plt.subplot(111,aspect = 1)
    ax.plot(x1,y1,'b')
    ax.plot(x2,y2,'r');
    #print(x1,x2,y1,y2)
    
    for y0,color in zip(y0s, colors):
        X = array([0.8,y0,u0,0])
        t = 0
        while t < T:
            h = c/LA.norm( f(X,t)[2:] )          # time-step ("delta-t")
            newX = rk4(f,X,t,h)
            plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
            X = newX
            t += h
plt.show();
T = 20
c = 1

colors = 'gkcmr'
u0s = [0.2,0.5,1.]
y0s = [0.05]

for u0 in u0s:
    t = linspace(0,2*pi, 100, endpoint = False)
    x1 = R1 * cos(t)
    y1 = R1 * sin(t)
    x2 = -R2 * cos(t)
    y2 = -R2 * sin(t)
    fig = plt.figure(figsize = (4,4))
    ax = plt.subplot(111,aspect = 1)
    ax.plot(x1,y1,'b')
    ax.plot(x2,y2,'r');

    for y0,color in zip(y0s, colors):
        X = array([0.8,y0,u0,0])
        t = 0
        while t < T:
            h = c/LA.norm( f(X,t)[2:] )          # time-step ("delta-t")
            newX = rk4(f,X,t,h)
            plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
            X = newX
            t += h
plt.show();
T = 20
c = 0.5

colors = 'gkcmr'
y0s = [0.1]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))
ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors):
    X = array([0.7,y0,sqrt(G*M1/y0),0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] )          # time-step ("delta-t")
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h
plt.show();
T = 200
c = 0.5

colors = 'gkcmr'
y0s = [0.1]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))
ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors):
    X = array([0.7,y0,sqrt(G*M/y0),0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] )          # time-step ("delta-t")
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h
plt.show();
T = 10000
c = 0.05

colors = 'gkcmr'
y0s = [-10.0,-20.0]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))

ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors):
    X = array([0,y0,sqrt(G*M/(abs(y0))),0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] )
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.show();
T = 100000
c = 0.05

colors = 'gkcmr'
y0s = [5.0,10.0,20.0]

t = linspace(0,2*pi, 100, endpoint = False)
x1 = R1 * cos(t)
y1 = R1 * sin(t)
x2 = -R2 * cos(t)
y2 = -R2 * sin(t)
fig = plt.figure(figsize = (4,4))

ax = plt.subplot(111,aspect = 1)
ax.plot(x1,y1,'b')
ax.plot(x2,y2,'r');

for y0,color in zip(y0s, colors):
    X = array([0,y0,-sqrt(G*M/(abs(y0))),0])
    t = 0
    while t < T:
        h = c/LA.norm( f(X,t)[2:] )
        newX = rk4(f,X,t,h)
        plt.plot([X[0],newX[0]],[X[1],newX[1]], color,ms =  2, alpha = 0.3);
        X = newX
        t += h
plt.xlim(-30,30)
plt.ylim(-30,30)
plt.show();
