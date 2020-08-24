"""
@authors:   Luis Juarez Mercedes
            Pedro Vidal Arias
"""

import Actuators
from numpy import pi, shape
from numpy import array, concatenate, zeros, degrees, savetxt
from Robot6_Modeling import P2P, line, circle, Tm
from Graficas import grafArt, Animate, grafArtAcum, grafFull, grafFullINTER
from time import sleep

global matQ     # Total Data Matrix
matQ = zeros((3,6,0))
global HOME     # HOME Position
x0 = 407.999;   r = 0
y0 = 0;         p = -90
z0 = 200;       y = 180
HOME = array([x0,y0,z0,r,p,y])
global actual
actual = array(HOME)

"""
Set angles in each Robot Actuator
"""

def send(mat_q, Tm, op):
    # Number of points
    n = shape(mat_q)[1]
    for i in range(n):    
        q = mat_q[:,i]*180/pi
        q[0] = (q[0] + 132)
        q[1] = (q[1] + 41)
        q[2] = (q[2] + 130)
        q[3] = (q[3] + 85)
        q[4] = (q[4] + 135)
        q[5] = (q[5] + 85)
        try:
            motor.servo[1].angle = float(q[0])
            motor.servo[2].angle = float(q[1])
            motor.servo[3].angle = float(q[2])
            motor.servo[4].angle = float(q[3])
            motor.servo[5].angle = float(q[4])
            motor.servo[6].angle = float(q[5])
        except:
            print('Out of range')
        # sampling Time
        sleep(Tm)

"""
Send PTP Trajectory
"""

def PTP(P,tiempo,v):   # PTP Motion
    global actual
    P0    = actual[0:3]
    rpy0  = actual[3:6]
    conf0 = [1,1,1]
    P1    = P[0:3]
    rpy1  = P[3:6]
    conf1 = [1,1,1]
    # -------------- Trajectory Planning ----------------
    mat_q = P2P(P0,rpy0,conf0,P1,rpy1,conf1,tiempo,v)
    # ------------------ Data Matrix --------------------
    actual[0:3] = P1
    actual[3:6] = rpy1
    global matQ
    matQ = concatenate((matQ,mat_q),axis=2)
    
"""
Enviar trayectorias en espacio de la tarea
"""

def LIN(P,tiempo,v):   # LIN Motion
    global actual
    P0    = actual[0:3]
    rpy   = actual[3:6]
    P1    = P[0:3]
    conf  = [1,1,1]
    # -------------- Trajectory Planning ----------------
    mat_q = line(P0,P1,rpy,conf,tiempo,v)
    # ------------------ Data Matrix --------------------
    actual[0:3] = P1
    actual[3:6] = P[3:6]
    global matQ
    matQ = concatenate((matQ,mat_q),axis=2)

def CIRC(Pa,P,tiempo,v):   # CIRC Motion
    global actual
    P0    = actual[0:3]
    rpy   = actual[3:6]
    conf  = [1,1,1]
    Pa    = Pa[0:3]
    P1    = P[0:3]
    # -------------- Trajectory Planning ----------------
    mat_q = circle(P0,Pa,P1,rpy,conf,tiempo,v)
    # ------------------ Data Matrix --------------------
    actual[0:3] = P1
    actual[3:6] = P[3:6]
    global matQ
    matQ = concatenate((matQ,mat_q),axis=2)

# Export Data for Solidworks Motion Analysis
def exportar(t,matq):
    for i in range(6):
        pos = degrees(matq[i,:])
        savetxt('q'+str(i+1)+'.txt',array([t,pos]).T,delimiter=',')




# -------------------------------------------------
# -------------- Programming Console --------------
# -------------------------------------------------
P1 = [120, 100, 300, 0, -90, 180]
P2 = [350,   0, 100, 0, -90, 180]
Pa = [205,-162, 182, 0, -90, 180]
P3 = [312,-154,  58, 0, -90, 180]
PTP(P1,3,0.625)
LIN(P2,3,0.625)
CIRC(Pa,P3,4,0.625)

# t,q,qp,qpp,tau = grafFull(matQ)
# exportar(t,q)














