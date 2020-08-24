"""
Created on Thu Mar 19 22:19:34 2020

@author: Luis Alexander Juarez Mercedes
"""

from Robot6_Modeling import *
from Robot6_Modeling import Dynamics, Tm
import numpy as np
from numpy import *
from numpy import pi, sqrt, array, linspace
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
params = {'text.latex.preamble' : [r'\usepackage{xfrac}']}
plt.rcParams.update(params)
params = {'xtick.labelsize': 15, 'ytick.labelsize': 15} # 30 en grafAcum
plt.rcParams.update(params)

"""
Graficar Variables articulares
"""

def graf(mat,letra):
    if letra == 'pos':
        # lb = 'q'
        # tl = 'Posición'
        var = ('q','[°]','Posición')
    elif letra == 'vel':
        # lb = 'qp'
        # tl = 'Velocidad'
        var = ('qp','[°/s]','Velocidad')
    elif letra == 'acel':
        # lb = 'qpp'
        # tl = 'Aceleración'
        var = ('q','[°/s^2]','Aceleración')
    elif letra == 'tau':
        # lb = 'T'
        # tl = 'Torque'
        var = ('T','[Nm]','Torque')
        
    q = np.degrees(mat) # Transformar a degress
    size = np.shape(q)  # Dimensiones de matriz
    t = np.arange(size[1])  # Vector X segun puntos
    pos = 0
    # GRAFICAS
    fig, ax = plt.subplots(3,2, figsize = (12,7), sharex = True)
    for i in range(3):
        for j in range(2):
            # Cálculos necesarios
            Q = q[pos,:]
            pos = pos + 1
            d = max(Q)-min(Q)   # Tamaño de grafica en Y
            if d == 0:  d = 10  # Si no hay variacion Tamaño en Y = 10
            xmin = t[0]
            xmax = t[-1]
            ymin = min(Q)-0.1*d
            ymax = max(Q)+0.1*d
            # Gráfica
            ax[i,j].plot(t,Q, linewidth = 3)
            if i == 2:  ax[i,j].set_xlabel('t[s]')  # Xlabel solo en ultima fila
            ax[i,j].set_ylabel(var[0] + str(pos) + var[1])
            ax[i,j].set_title(var[2] + ' de Art ' + str(pos), fontSize=15)
            ax[i,j].axis([xmin,xmax,ymin,ymax])
            ax[i,j].grid(linestyle='--')
    plt.show()

def grafArt(mat_q):   # Graficar Matrices de datos
    q = mat_q[0,:,:]
    qp = mat_q[1,:,:]
    qpp = mat_q[2,:,:]
    graf(q,'pos')
    graf(qp,'vel')
    graf(qpp,'acel')
    # Dinámica
    matT = zeros((6,q.shape[1]))
    for i in range(q.shape[1]):
        # A06 = cinematicaDirecta(q[:,i])
        # R = A06[0:3,0:3]
        F0 = [0,0,0.5]    # [kg]  Fuerza espresado en el sistema Base
        M0 = [0,0,0]    # [kg.cm]  Momento expresado en el sistema Base
        # Fe = R.T@F0*9.81        # fuerza en el cuerpo
        # Me = R.T@M0*9.81/100    # Momento en el cuerpo
        matT[:,i] = Dynamics(q[:,i],qp[:,i],qpp[:,i],Fe,Me)[6]*pi/180
    graf(matT,'tau')
    return q, qp, qpp, matT

def grafAcum(mat,letra):
    if letra == 'pos':
        var = (r'Posición [$^\circ$]','Posición',r'${q_{')
    elif letra == 'vel':
        var = (r'Velocidad [$\sfrac{^\circ}{s}$]','Velocidad',r'$\dot{q_{')
    elif letra == 'acel':
        var = (r'Aceleración [$\sfrac{^\circ}{s^2}$]','Aceleración',r'$\ddot{q_{')
    elif letra == 'tau':
        var = ('Torque [kg$\cdot$cm]','Torque',r'${\tau_{')
        
    Q = np.degrees(mat) # Transformar a degress
    # if letra == 'pos':
    #     Q[1,:] += 90
    #     Q[2,:] -= 90
    size = np.shape(Q)  # Dimensiones de matriz
    t = linspace(0,Tm*size[1],size[1])  # Vector X segun puntos
    pos = 0
    # GRAFICAS
    fig, ax = plt.subplots(1,1, figsize = (12,9), sharex = True)
    # Cálculos necesarios
    pos = pos + 1
    d = np.amax(Q) - np.amin(Q)   # Tamaño de grafica en Y
    if d == 0:  d = 10  # Si no hay variacion Tamaño en Y = 10
    xmin = t[0]
    xmax = t[-1]
    ymin = np.amin(Q) -0.1*d
    ymax = np.amax(Q) + 0.1*d
    # Gráfica
    for i in range(6):
        q = Q[i,:]
        ax.plot(t,q, linewidth = 2, label = var[2]+str(i+1)+'}}$')
    ax.set_xlabel('Tiempo [s]', fontSize=30)  # Xlabel solo en ultima fila
    ax.set_ylabel(var[0], fontSize=30)
    # ax.set_title(var[1] + ' vs Tiempo', fontSize=30)
    ax.axis([xmin,xmax,ymin,ymax])
    ax.grid(linestyle='--')
    ax.legend(prop={'size': 20})
    plt.show()
    plt.savefig("CIRC-"+var[1]+".pdf")
    return t

def grafArtAcum(mat_q,modo='plot'):   # Graficar Matrices de datos
    q = mat_q[0,:,:]
    qp = mat_q[1,:,:]
    qpp = mat_q[2,:,:]
    t = grafAcum(q,'pos')
    grafAcum(qp,'vel')
    grafAcum(qpp,'acel')
    # Dinámica
    matT = zeros((6,q.shape[1]))
    for i in range(q.shape[1]):
        # A06 = cinematicaDirecta(q[:,i])
        # R = A06[0:3,0:3]
        F0 = [0,0,0.5]    # [kg]  Fuerza espresado en el sistema Base
        M0 = [0,0,0]    # [kg.cm]  Momento expresado en el sistema Base
        # Fe = R.T@F0*9.81        # fuerza en el cuerpo
        # Me = R.T@M0*9.81/100    # Momento en el cuerpo
        matT[:,i] = Dynamics6(q[:,i],qp[:,i],qpp[:,i],F0,M0)[6]*pi/180*100/9.81
    grafAcum(matT,'tau')
    return t, q, qp, qpp, matT


def grafFull(mat):
    q = mat[0,:,:]
    qp = mat[1,:,:]
    qpp = mat[2,:,:]
    matT = zeros((1,6,q.shape[1]))
    for i in range(mat.shape[2]):
        # A06 = cinematicaDirecta(q[:,i])
        # R = A06[0:3,0:3]
        F0 = [0,0,0.5]    # [kg]  Fuerza espresado en el sistema Base
        M0 = [0,0,0]    # [kg.cm]  Momento expresado en el sistema Base
        # Fe = R.T@F0*9.81        # fuerza en el cuerpo
        # Me = R.T@M0*9.81/100    # Momento en el cuerpo
        matT[0,:,i] = Dynamics(q[:,i],qp[:,i],qpp[:,i],F0,M0)[6]*pi/180*100/9.81
    mat = concatenate((mat,matT),axis=0)
    #----------------------------------------------------------
    # GRAFICAR
    size = np.shape(mat[0,:,:])  # Dimensiones de matriz
    t = linspace(0,Tm*size[1],size[1])  # Vector X segun puntos
    # GRAFICAS
    fig, ax = plt.subplots(4,1, figsize = (12,15), sharex = False)
    # Cálculos necesarios
    for i in range(4):
        for j in range(1):
            if i==0:
                var = (r'Position [$^\circ$]','Posición',r'${q_{')
            elif i==1:
                var = (r'Velocity [$\sfrac{^\circ}{s}$]','Velocidad',r'$\dot{q_{')
            elif i==2:
                var = (r'Acceleration [$\sfrac{^\circ}{s^2}$]','Aceleración',r'$\ddot{q_{')
            elif i==3:
                var = ('Torque [kg$\cdot$cm]','Torque',r'${\tau_{')
            
            Q = np.degrees(mat[i,:,:]) # Transformar a degress
            # if i ==0:
            #     Q[1,:] += 90
            #     Q[2,:] -= 90
            d = np.amax(Q) - np.amin(Q)   # Tamaño de grafica en Y
            if d == 0:  d = 10  # Si no hay variacion Tamaño en Y = 10
            xmin = t[0]
            xmax = t[-1]
            ymin = np.amin(Q) - 0.1*d
            ymax = np.amax(Q) + 0.1*d
            # Gráfica
            for k in range(6):
                q = Q[k,:]
                ax[i].plot(t,q, linewidth = 1, label = var[2]+str(k+1)+'}}$')
                # tiempo = [int(3/Tm),int(6/Tm)]
                # points = q[tiempo]
                # ax[i].scatter([3,6],points, s=15)
            ax[3].set_xlabel('Time [s]', fontSize=15)  # Xlabel solo en ultima fila
            ax[i].set_ylabel(var[0], fontSize=15)
            ax[i].axis([xmin,xmax,ymin,ymax])
            ax[i].grid(linestyle='--')
            ax[0].set_title('Joint variables vs Time', fontSize=20)
            ax[i].legend(prop={'size': 11})
    plt.show()
    plt.savefig("GRAF-TRAYECTO_ANDESCON.pdf")
    return t, mat[0,:,:], mat[1,:,:], mat[2,:,:], mat[3,:,:]


def grafFullINTER(mat):
    q = mat[0,:,:]
    qp = mat[1,:,:]
    qpp = mat[2,:,:]
    matT = zeros((1,6,q.shape[1]))
    for i in range(mat.shape[2]):
        # A06 = cinematicaDirecta(q[:,i])
        # R = A06[0:3,0:3]
        F0 = [0,0,0.5]    # [kg]  Fuerza espresado en el sistema Base
        M0 = [0,0,0]    # [kg.cm]  Momento expresado en el sistema Base
        # Fe = R.T@F0*9.81        # fuerza en el cuerpo
        # Me = R.T@M0*9.81/100    # Momento en el cuerpo
        matT[0,:,i] = Dynamics6(q[:,i],qp[:,i],qpp[:,i],F0,M0)[6]*pi/180*100/9.81
    mat = concatenate((mat,matT),axis=0)
    #----------------------------------------------------------
    # GRAFICAR
    size = np.shape(mat[0,:,:])  # Dimensiones de matriz
    t = linspace(0,Tm*size[1],size[1])  # Vector X segun puntos
    # GRAFICAS
    fig, ax = plt.subplots(3,1, figsize = (7,8), sharex = False)
    # Cálculos necesarios
    for i in range(3):
        if i == 2:
            K = 3
        else:
            K = i
        for j in range(1):
            if K==0:
                var = (r'Position [$^\circ$]','Posición',r'${q_{')
            elif K==1:
                var = (r'Velocity [$\sfrac{^\circ}{s}$]','Velocidad',r'$\dot{q_{')
            elif K==2:
                var = (r'Acceleration [$\sfrac{^\circ}{s^2}$]','Aceleración',r'$\ddot{q_{')
            elif K==3:
                var = ('Torque [kg$\cdot$cm]','Torque',r'${\tau_{')
            
            Q = np.degrees(mat[K,:,:]) # Transformar a degress
            # if i ==0:
            #     Q[1,:] += 90
            #     Q[2,:] -= 90
            d = np.amax(Q) - np.amin(Q)   # Tamaño de grafica en Y
            if d == 0:  d = 10  # Si no hay variacion Tamaño en Y = 10
            xmin = t[0]
            xmax = t[-1]
            ymin = np.amin(Q) - 0.1*d
            ymax = np.amax(Q) + 0.1*d
            # Gráfica
            for k in range(6):
                q = Q[k,:]
                ax[i].plot(t,q, linewidth = 1, label = var[2]+str(k+1)+'}}$')
            ax[2].set_xlabel('Time [s]', fontSize=15)  # Xlabel solo en ultima fila
            ax[i].set_ylabel(var[0], fontSize=15)
            ax[i].axis([xmin,xmax,ymin,ymax])
            ax[i].grid(linestyle='--')
            ax[0].set_title('Joint variables vs Time', fontSize=20)
            ax[i].legend(prop={'size': 10})
    plt.show()
    plt.savefig("GRAF-TRAYECTO_INTERCON.pdf")
    return t, mat[0,:,:], mat[1,:,:], mat[2,:,:], mat[3,:,:]















def plot3(a,b,c,mark="o",col="r"):
  from matplotlib import pyplot
  import pylab
  from mpl_toolkits.mplot3d import Axes3D
  pylab.ion()
  fig = pylab.figure()
  ax = Axes3D(fig)
  ax.scatter(a, b, c,marker=mark,color=col)

    
def f(mat_Q):
    mat_q = mat_Q[0,:,:]
    n = shape(mat_q)[1]
    
    # Longitudes del Robot6
    L1 = 200    # [mm]
    L2 = 162    # [mm]
    L3 = 160    # [mm]
    L4 =  86    # [mm]
    # Parámetros Denavit-Hartenberg del Robot6
    d        = [   L1,    0,         0,    L3,    0,   L4]
    a        = [    0,   L2,         0,     0,    0,    0]
    alpha    = [ pi/2,    0,      pi/2, -pi/2, pi/2,    0]
    
    # # Longitudes del Robot
    # L1 =   329  # 199    # [mm]
    # L2 =    50  #   0    # [mm]
    # L3 =   330  # 162    # [mm]
    # L4 = 35.03  #   0    # [mm]
    # L5 = 333.4  # 160    # [mm]
    # L6 =    80  # 120    # [mm]
    # # Parámetros Denavit-Hartenberg del robot
    # d       = [   L1,         0,    0,    L5,    0,   L6]
    # a       = [   L2,        L3,   L4,     0,    0,    0]
    # alpha   = [ pi/2,         0, pi/2, -pi/2, pi/2,    0]
    # Vector de posicion (x, y, z) del sistema de coordenadas de referencia
    x0 = 0; y0 = 0; z0 = 0
    
    
    lim_x = L2 + L3 + L4 + L5 + L6 + 50
    lim_y = lim_x
    lim_z = L1 + L3 + L4 + L5 + L6 + 50
    
    x = zeros((n,8));   X = zeros(n)
    y = zeros((n,8));   Y = zeros(n)
    z = zeros((n,8));   Z = zeros(n)
    
    MTHEF = eye(4);
    MTH1 = eye(4);
    
    for i in range(n):
        # Variables articulares del brazo robot
        teta1 = mat_q[0,i]
        teta2 = mat_q[1,i] + pi/2
        teta3 = mat_q[2,i]
        teta4 = mat_q[3,i]
        teta5 = mat_q[4,i]
        teta6 = mat_q[5,i]
        # Matrices de transformación homogénea entre sistemas de coordenadasconsecutivos
        A01 = denavit(teta1, d[0], a[0], alpha[0])
        A0X = denavit(teta1, d[0], 0, 0)
        A12 = denavit(teta2, d[1], a[1], alpha[1])
        A23 = denavit(teta3, d[2], a[2], alpha[2])
        A34 = denavit(teta4, d[3], a[3], alpha[3])
        A45 = denavit(teta5, d[4], a[4], alpha[4])
        A56 = denavit(teta6, d[5], a[5], alpha[5])
        # Matrices de transformación del primer sistema al correspondiente
        A02 = A01 @ A12
        A03 = A02 @ A23
        A04 = A03 @ A34
        A05 = A04 @ A45
        A06 = A05 @ A56
        AUX = A06 @ MTH1
        AEF = A06 @ MTHEF
        # Vectores direccion de cada sistea de coordenadas
        nx = AEF[0,0]; ny = AEF[1,0]; nz = AEF[2,0]
        ox = AEF[0,1]; oy = AEF[1,1]; oz = AEF[2,1]
        ax = AEF[0,2]; ay = AEF[1,2]; az = AEF[2,2]
        # Vector de posicion (x, y, z) de cada sistema de coordenadas
        xX = A0X[0,3]; yX = A0X[1,3]; zX = A0X[2,3]
        x1 = A01[0,3]; y1 = A01[1,3]; z1 = A01[2,3]
        x2 = A02[0,3]; y2 = A02[1,3]; z2 = A02[2,3]
        x3 = A03[0,3]; y3 = A03[1,3]; z3 = A03[2,3]
        x4 = A04[0,3]; y4 = A04[1,3]; z4 = A04[2,3]
        x5 = A05[0,3]; y5 = A05[1,3]; z5 = A05[2,3]
        x6 = A06[0,3]; y6 = A06[1,3]; z6 = A06[2,3]
        xEF = AEF[0,3]; yEF = AEF[1,3]; zEF = AEF[2,3]
        # Crear curva de trayectoria
        X[i] = xEF
        Y[i] = yEF
        Z[i] = zEF
        # Se dibuja el robot
        x[i,:] = array([x0, xX, x1, x2, x3, x4, x5, x6])
        y[i,:] = array([y0, yX, y1, y2, y3, y4, y5, y6])
        z[i,:] = array([z0, zX, z1, z2, z3, z4, z5, z6])
    return x,y,z,n




def Animate(mat_Q):
    X, Y, Z, n = f(mat_Q)
    
    import plotly.graph_objects as go
    import numpy as np
    np.random.seed(1)
    
    # Longitudes del Robot
    L1 =   329  # 199    # [mm]
    L2 =    50  #   0    # [mm]
    L3 =   330  # 162    # [mm]
    L4 = 35.03  #   0    # [mm]
    L5 = 333.4  # 160    # [mm]
    L6 =    80  # 120    # [mm]
    
    lim_x = L2 + L3 + L4 + L5 + L6 + 50
    lim_y = lim_x
    lim_z = L1 + L3 + L4 + L5 + L6 + 50
    
    x0 = 0;         y0 = 0;     z0 = 0
    x1 = 0;         y1 = 0;     z1 = L1
    x2 = L2;        y2 = 0;     z2 = z1
    x3 = x2;        y3 = 0;     z3 = z2 + L3
    x4 = x3;        y4 = 0;     z4 = z3 + L4
    x5 = x4 + L5;   y5 = 0;     z5 = z4
    x6 = x5 + L6;   y6 = 0;     z6 = z5
    
    Robot = go.Scatter3d(
        x = X[0,:],
        y = Y[0,:],
        z = Z[0,:],
        mode = 'lines+markers',
        marker = dict(
            size = 5,
            color = 'red',                # set color to an array/list of desired values
            colorscale = 'Viridis',   # choose a colorscale
            opacity = 0.8
        ),
        line = dict(
            width = 5,
            color = 'blue',
        ),
    )
    
    fig = go.Figure(
        frames=[go.Frame(
            data=go.Scatter3d(
                x = X[k,:],
                y = Y[k,:],
                z = Z[k,:],
                mode = 'lines+markers',
                marker = dict(
                    size = 5,
                    color = 'red',                # set color to an array/list of desired values
                    colorscale = 'Viridis',   # choose a colorscale
                    opacity = 0.8
                ),
                line = dict(
                    width = 5,
                    color = 'blue',
                ),
            )) for k in range(n)]
    )
    
    fig.add_trace(Robot)
    fig.layout.title = 'Robot 6'
    
    fig.update_layout(
        updatemenus=[dict(type="buttons",
                              buttons=[dict(label="Play",
                                            method="animate",
                                            args=[None])])],
        scene = dict(
            xaxis_title='X AXIS TITLE',
            yaxis_title='Y AXIS TITLE',
            zaxis_title='Z AXIS TITLE',
            xaxis = dict(nticks=4, range=[-lim_x,lim_x],),
            yaxis = dict(nticks=4, range=[-lim_y,lim_y],),
            zaxis = dict(nticks=4, range=[     0,lim_z],),
            camera=dict(
                center = dict(x=0, y=0, z=0),
                eye = dict(x=1, y=1.2, z=0.6),
                up = dict(x=0, y=0, z=1),
            ),
            annotations=[dict(
                x = X[0,0],
                y = Y[0,0],
                z = Z[0,0],
                ax = 50,
                ay = 0,
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Base",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="left",
                yanchor="bottom"
            ),dict(
                x = X[0,2],
                y = Y[0,2],
                z = Z[0,2],
                ax = 48,
                ay = -sqrt(50**2-48**2),
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Joint1",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="left",
                yanchor="bottom"
            ),dict(
                x = X[0,3],
                y = Y[0,3],
                z = Z[0,3],
                ax = -45,
                ay = sqrt(50**2-45**2),
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Joint2",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="right",
                yanchor="bottom"
            ),dict(
                x = X[0,4],
                y = Y[0,4],
                z = Z[0,4],
                ax = 45,
                ay = -sqrt(50**2-45**2),
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Joint3",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="left",
                yanchor="bottom"
            ),dict(
                x = X[0,5],
                y = Y[0,5],
                z = Z[0,5],
                ax = 35,
                ay = -sqrt(50**2-35**2),
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Joint4",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="left",
                yanchor="bottom"
            ),dict(
                x = X[0,6],
                y = Y[0,6],
                z = Z[0,6],
                ax = -20,
                ay = -sqrt(50**2-20**2),
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Joint3{S3}",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="left",
                yanchor="bottom"
            ),dict(
                x = X[0,7],
                y = Y[0,7],
                z = Z[0,7],
                ax = -45,
                ay = -sqrt(50**2-45**2),
                font = dict(
                    color = "black",
                    size = 12
                ),
                text = "Joint3{S3}",
                arrowcolor = "black",
                arrowsize = 1.4,
                arrowwidth = 1.5,
                arrowhead = 1,
                xanchor="right",
                yanchor="bottom"
            ),],
                
        ),
        scene_aspectmode='manual',
        scene_aspectratio = dict(x = 1, y = 1, z = lim_x/lim_z),
        title_text = 'Robot6',
    )
         
    fig.show()

  




