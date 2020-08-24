"""
@authors:   Luis Juarez Mercedes
            Pedro Vidal Arias
"""

from numpy import *
from numpy import cos, sin, degrees, radians, array, sqrt, pi
from numpy import eye, dot, cross, concatenate
from numpy import zeros, linspace
from numpy.linalg import inv, norm
from math import atan2, acos
import numpy as np
np.set_printoptions(precision=4,suppress=True)

Tm = 0.01   # Sampling time

# Robot lengths [mm]
L1 = 200;
L2 = 162;
L3 = 160;
L4 =  86;
# Link masses [kg]
m1 =  42.75e-3
m2 = 419.92e-3
m3 =  98.45e-3
m4 = 112.32e-3
m5 =  78.64e-3
m6 =   1.49e-3
m = [m1,m2,m3,m4,m5,m6]
# Centers of mass [m]
c1 = array([  0.00, -17.85, -0.87])*1e-3
c2 = array([-80.77,   0.00, -0.09])*1e-3
c3 = array([  0.00,  -6.08, 42.48])*1e-3
c4 = array([  0.00,  28.21, -0.18])*1e-3
c5 = array([  0.00,  -8.04, 43.37])*1e-3
c6 = array([  0.00,   0.00, -1.14])*1e-3
rc = [c1,c2,c3,c4,c5,c6]
# Inertia tensors [kg.m^2]
I1 = array([[ 38089.06,       0.00,       0.00],
            [     0.00,   28794.35,    -663.45],
            [     0.00,    -663.45,   15228.60]])*1e-6

I2 = array([[123753.47,       0.00,   -1432.62],
            [     0.00, 1468922.83,       0.00],
            [ -1432.62,    -663.45, 1407205.97]])*1e-6

I3 = array([[98943.42,        0.00,       0.00],
            [    0.00,    67176.78,  -11572.80],
            [    0.00,   -11572.80,   40844.44]])*1e-6

I4 = array([[87863.58,        0.02,       0.00],
            [    0.02,    22365.32,    -232.28],
            [    0.00,     -232.28,   72919.07]])*1e-6

I5 = array([[43833.38,        0.00,       0.00],
            [    0.00,    29183.59,   -4439.40],
            [    0.00,    -4439.40,   19702.24]])*1e-6

I6 = array([[   37.32,        0.00,       0.00],
            [    0.00,       37.32,       0.00],
            [    0.00,        0.00,      72.24]])*1e-6
I = [I1, I2, I3, I4, I5, I6]
# Gravity acceleration
g = 9.81
g0  = array([0,0,-g])

"""
Additional Functions
"""    

def zeroEval(num):  # Zero evaluation
    if abs(num) <= 0.0000007:
        return 0
    else:
        return num

def unit(v):    # Unit vector
    modulo = np.linalg.norm(v)
    if modulo == 0:
        return v
    else:
        return v/modulo


"""
Orientation Representation
"""

def rotx(theta):    # Rotation X
    angle = radians(theta)
    return array([[ 1,      0,               0],
                  [ 0, cos(angle), -sin(angle)],
                  [ 0, sin(angle),  cos(angle)]])

def roty(theta):    # Rotation Y
    angle = radians(theta)
    return array([[  cos(angle), 0, sin(angle)],
                  [           0, 1,          0],
                  [ -sin(angle), 0, cos(angle)]])

def rotz(theta):    # Rotation Z
    angle = radians(theta)
    return array([[ cos(angle), -sin(angle), 0],
                  [ sin(angle),  cos(angle), 0],
                  [          0,           0, 1]])

def rpy2mat(rpy):   # Euler Angles XYZ
    return rotz(rpy[2]) @ roty(rpy[1]) @ rotx(rpy[0])

def mat2rpy(R):     # Rotation Matrix to Euler Angles XYZ
    Cp = zeros(4)
    Cp[0] = sqrt((R[2,1])**2+(R[2,2])**2)
    Cp[1] = -Cp[0]
    Cp[2] = sqrt((R[1,0])**2+(R[0,0])**2)
    Cp[3] = -Cp[2]
    Sp = -R[2,0]
    if R[2,0] == 0:
        Sp = 0
    rpy = zeros((4,3))
    for i in range(4):
        p = atan2(Sp,Cp[i])
        r = atan2(R[2,1]/Cp[i],R[2,2]/Cp[i])
        y = atan2(R[1,0]/Cp[i],R[0,0]/Cp[i])
        rpy[i,:] = degrees(array([r,p,y]))
    return rpy

def rot(k,theta):   # Rotation pair
    r = unit(k)
    rx = r[0]
    ry = r[1]
    rz = r[2]
    C = cos(radians(theta))
    S = sin(radians(theta))
    R = array([[    rx**2*(1-C)+C,     rx*ry*(1-C)-rz*S,   rx*rz*(1-C)+ry*S],
               [ rx*ry*(1-C)+rz*S,        ry**2*(1-C)+C,   ry*rz*(1-C)-rx*S],
               [ rx*rz*(1-C)-ry*S,     ry*rz*(1-C)+rx*S,      rz**2*(1-C)+C]])
    return R

def rotAxis(p,k,angle): # Rotation pair
    theta = radians(angle)
    r = p*cos(theta) + cross(k,p)*sin(theta) + k*dot(k,p)*(1-cos(theta))
    return r

"""
Denavit-Hartenberg Algorithm
"""

def dh(theta,d,a,alpha):
    T = array([[ cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta)],
               [ sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
               [          0,             sin(alpha),             cos(alpha),            d],
               [          0,                      0,                      0,            1]])
    return T


"""
Robot Kinematics
"""

def directKinematic(q):
    # Denavit-Hartenberg parameters
    theta    = [ q[0], q[1], q[2]+pi/2,  q[3], q[4], q[5]]
    d        = [   L1,    0,         0,    L3,    0,   L4]
    a        = [    0,   L2,         0,     0,    0,    0]
    alpha    = [ pi/2,    0,      pi/2, -pi/2, pi/2,    0]
    # Homogeneous Transformation Matrices
    A01 = dh(theta[0], d[0], a[0], alpha[0])
    A12 = dh(theta[1], d[1], a[1], alpha[1])
    A23 = dh(theta[2], d[2], a[2], alpha[2])
    A34 = dh(theta[3], d[3], a[3], alpha[3])
    A45 = dh(theta[4], d[4], a[4], alpha[4])
    A56 = dh(theta[5], d[5], a[5], alpha[5])
    # {S6} respect to {S0} HTM
    A06 = A01 @ A12 @ A23 @ A34 @ A45 @ A56
    return A06


def inverseKinematic(T,shoulder,elbow,wrist):
    # ---------------------------------------------
    # ********* First 3 Degrees of Freedom ********
    # ---------------------------------------------
    P = T[0:3,3]     # End Efector position
    z6 = T[0:3,2]    # z6 vector
    # Kinematic decoupling
    Pm = P - L4*z6  # Wrist position
    # Robot Configuration
    shoulder = -1*(-1)**shoulder    # Shoulder
    elbow = (-1)**elbow
    elbow = elbow*shoulder          # Elbow
    # Zero evaluation of wrist position
    x = zeroEval((Pm[0]))
    y = zeroEval((Pm[1]))
    z = zeroEval((Pm[2]))
    # Required computation
    R = sqrt(x**2 + y**2)
    r = sqrt(R**2 + (z-L1)**2)
    C3 = (r**2 - L2**2 - L3**2)/(2*L2*L3)
    S3 = elbow*sqrt(1 - C3**2)
    # Zero evaluation
    S3 = zeroEval(S3)
    C3 = zeroEval(C3)
    # Auxiliar
    fhi = atan2(L3*S3, L2+L3*C3)
    # Arm joints
    q1 = atan2(shoulder*y,shoulder*x)
    q3 = atan2(S3,C3)
    q2 = atan2(z-L1,shoulder*R) - fhi
    q2 = zeroEval(q2)
    # --------------------------------------------
    # ********* Last 3 Degrees of Freedom ********
    # --------------------------------------------
    muneca = -1*(-1)**wrist
    # Required computation
    A01 = dh(      q1, L1,  0, pi/2)
    A12 = dh(      q2,  0, L2,    0)
    A23 = dh( q3+pi/2,  0,  0, pi/2)
    A03 = A01 @ A12 @ A23
    R03 = A03[0:3,0:3]
    R06 = T[0:3,0:3]
    # {S6} respect to {S0} Rotation Matrix
    R36 = (R03.T) @ R06
    # Required computation
    C5 = R36[2,2]
    C5 = zeroEval(C5)
    S5 = muneca*sqrt(1 - C5**2)
    S5 = zeroEval(S5)
    S4 = R36[1,2]/S5
    C4 = R36[0,2]/S5
    S6 = R36[2,1]/S5
    C6 = -R36[2,0]/S5
    # Zero evaluation
    S4 = zeroEval(S4)
    C4 = zeroEval(C4)
    S6 = zeroEval(S6)
    C6 = zeroEval(C6)
    # Wrist joints
    q4 = atan2(S4,C4)
    q5 = atan2(S5,C5)
    q6 = atan2(S6,C6)
    # Joint vector
    q = array([q1, q2, q3, q4, q5, q6])
    return q

"""
Inverse Kinematic P[x,y,z,r,p,y]
"""

def IK(P,rpy,conf):
    T = np.eye(4)
    R = rpy2mat(rpy)
    T[0:3,0:3] = R
    T[0:3,3] = P
    q = inverseKinematic(T,conf[0],conf[1],conf[2])
    return q


"""
Point-to-Point Trajectory Planning
"""

def trajectoryVT(q0,qf,t,v):
    tf = t
    deltaQ = qf - q0
    qpc_max = 2*(deltaQ)/tf
    qpc = v*qpc_max
    tc = 0
    for k in range(6):
        if deltaQ[k] != 0:
            tc = tf - (deltaQ)[k]/qpc[k]
            break
    tb = tf - tc
    t = linspace(0,tf,int(tf/Tm))
    intervals = [t<tc, (tc<=t)&(t<tb), t>=tb]
    mat_q = zeros((3,6,int(tf/Tm)))
    # Solution
    for i in range(6):
        if deltaQ[i] != 0:
            positions = [lambda t: q0[i] + 1/2*qpc[i]/tc*t**2,
                         lambda t: q0[i] + qpc[i]*(t - tc/2),
                         lambda t: qf[i] - 1/2*qpc[i]/tc*(tf-t)**2]
            velocities = [lambda t: qpc[i]/tc*t,
                          lambda t: qpc[i],
                          lambda t: qpc[i]/tc*(tf-t)]
            accelerations = [lambda t:  qpc[i]/tc,
                             lambda t:  0,
                             lambda t: -qpc[i]/tc]
            q = piecewise(t, intervals, positions)
            qp = piecewise(t, intervals, velocities)
            qpp = piecewise(t, intervals, accelerations)
        # Data Matrix
            mat_q[0,i,:] = q
            mat_q[1,i,:] = qp
            mat_q[2,i,:] = qpp
    return mat_q

# POINT 2 POINT: Generates an articular trajectory given an initial and final
#                position in the task space.
#       P0      : Initial Position
#       rpy0    : Initial Orientation
#       conf0   : Initial Configuration
#       P1      : Final Position
#       rpy1    : Final Orientation
#       conf1   : Final Configuration
#       t       : Time
#       v       : Percentage of Maximum Velocity allowed
def P2P(P0,rpy0,conf0,P1,rpy1,conf1,t,v):
    # Initial and final conditions
    qp0 = zeros(6)
    qpf = zeros(6)
    qpp0 = zeros(6)
    qppf = zeros(6)
    # Initial articular configuration
    q0 =  IK(P0,rpy0,conf0)
    # Final articular configuration
    qf =  IK(P1,rpy1,conf1)
    # Trajectory Planning
    mat_q = trajectoryVT(q0,qf,t,v)
    return mat_q


"""
Operational Trajectories
"""

def trapezoidal(L,t,v):
    tf = t
    mat_o = zeros((3,int(tf/Tm)))
    if L != 0:
        Vc_max = 2*L/tf
        Vc = v*Vc_max
        tc = 0
        tc = tf - L/Vc
        tb = tf - tc
        t = linspace(0,tf,int(tf/Tm))
        intervals = [t<tc, (tc<=t)&(t<tb), t>=tb]
        # Solution
        positions = [lambda t: 1/2*Vc/tc*t**2,
                     lambda t: Vc*(t - tc/2),
                     lambda t: L - 1/2*Vc/tc*(tf-t)**2]
        velocities = [lambda t: Vc/tc*t,
                      lambda t: Vc,
                      lambda t: Vc/tc*(tf-t)]
        accelerations = [lambda t:  Vc/tc,
                         lambda t:  0,
                         lambda t: -Vc/tc]
        p = piecewise(t, intervals, positions)
        v = piecewise(t, intervals, velocities)
        a = piecewise(t, intervals, accelerations)
        # Data Matrix
        mat_o[0,:] = p
        mat_o[1,:] = v
        mat_o[2,:] = a
    return mat_o

# RECTILINEAR LINE: Generates a Linear Trajectory given an initial and final
#                   position in the task space.
#       P0      : Initial Position
#       P1      : Final Position
#       rpy     : Constant Orientation
#       conf    : Constant Configuration
#       t       : Time
#       v       : Percentage of Maximum Velocity allowed
def line(P0,P1,rpy,conf,t,v):
    # Generate Linear Path
    if np.all(P0 == P1):    # {S6} doesn't move
        L = 0
        u = array([1,1,1])
    else:       # Else Computation
        dP = array(P1) - P0
        L = sqrt(dP[0]**2+dP[1]**2+dP[2]**2)
        u = dP/L    # Motion direction
    # ------------------ Interpolation -----------------
    s = trapezoidal(L,t,v)
    # ------------------ Path Sampling -----------------
    puntos = int(t/Tm)
    mat_X = zeros((3,6,puntos))
    for i in range(puntos):
        p = P0+(s[0,i])*u   # Cartesian Position
        Q = IK(p,rpy,conf)  # Joint Space
        mat_X[0,:,i] = Q    # Data Matrix

    for j in range(6):
        mat_X[1,j,1:] = np.diff(mat_X[0,j,:])
        mat_X[1,j,:] = mat_X[1,j,:]/Tm
        mat_X[2,j,1:] = np.diff(mat_X[1,j,:])
        mat_X[2,j,:] = mat_X[2,j,:]/Tm
    return mat_X

# CIRCLE: Generates a Circular Trajectory given an initial, auxiliar and final
#         position in the task space.
#       P0      : Posicion inicial
#       PA      : Posicion auxiliar
#       P1      : Posicion final
#       rpy     : Orientación inicial
#       conf    : Configuración de Robot
#       t       : Time
#       v       : Percentage of Maximum Velocity allowed
def circle(P0,PA,P1,rpy,conf,t,v):
    # Generate Circular Path
    P = array(P0)
    Q = array(PA)
    R = array(P1)
    PQ = Q - P
    PR = R - P
    n = cross(PQ,PR)
    k = unit(n)     # Rotation Axis
    # Coeficients Matrix
    A = array([[Q[0]-P[0], Q[1]-P[1], Q[2]-P[2]],
               [R[0]-P[0], R[1]-P[1], R[2]-P[2]],
               [     n[0],      n[1],      n[2]]])
    # Column Vector
    B = array([(norm(Q)**2-norm(P)**2)/2,
               (norm(R)**2-norm(P)**2)/2,
               n[0]*P[0]+n[1]*P[1]+n[2]*P[2]])

    C = (inv(A) @ B)    # Center of Circle
    r = P - C           # Radio Vector
    ri = P - C          # Initial Radio Vector
    ra = Q - C          # auxiliar Radio Vector
    rf = R - C          # Final Radio Vector
    # Validates angles over 180 degrees (pass througth auxiliar point)
    alpha1 = acos(dot(ri,ra)/(norm(ri)*norm(ra)))
    alpha2 = acos(dot(ra,rf)/(norm(ra)*norm(rf)))
    alpha3 = acos(dot(ri,rf)/(norm(ri)*norm(rf)))
    if round(alpha1+alpha2,3) == alpha3:
        alpha = alpha3
    else:
        alpha = alpha1 + alpha2
    alpha = degrees(alpha)
    # ------------------ Interpolation -----------------
    angle = trapezoidal(alpha,t,v)
    # ------------------ Path Sampling -----------------
    puntos = int(t/Tm)
    mat_X = zeros((3,6,puntos))
    for i in range(puntos):
        p = C + rotAxis(r,k,angle[0,i]) # Cartesian Position
        Q = IK(p,rpy,conf)              # Joint Space
        mat_X[0,:,i] = Q                # Data Matrix
        
    for j in range(6):
        mat_X[1,j,1:] = np.diff(mat_X[0,j,:])
        mat_X[1,j,:] = mat_X[1,j,:]/Tm
        mat_X[2,j,1:] = np.diff(mat_X[1,j,:])
        mat_X[2,j,:] = mat_X[2,j,:]/Tm
    return mat_X
        
"""
Dynamics
"""

def Dynamics(q,qp,qpp,F0,M0):
    A06 = directKinematic(q)    # End Effector Position
    R06 = A06[0:3,0:3]          # End Effector Rotation
    Fe = R06.T@F0*9.81          # External Force {S6}
    Me = R06.T@M0*9.81/100      # External Momentum {S6}
    #-------------------------------------------------------------------------
    n = len(q)
    F = [0]*(n+1);  F[-1] = Fe
    M = [0]*(n+1);  M[-1] = Me
    tau = [0]*2*n
    #-------------------------------------------------------------------------
    # Robot lengths [m]
    L1 = 200e-3
    L2 = 162e-3
    L3 = 160e-3
    L4 =  86e-3
    # Denavit-Hartenberg parameters
    theta    = [ q[0], q[1], q[2]+pi/2,  q[3], q[4], q[5]]
    d        = [   L1,    0,         0,    L3,    0,   L4]
    a        = [    0,   L2,         0,     0,    0,    0]
    alpha    = [ pi/2,    0,      pi/2, -pi/2, pi/2,    0]
    #-------------------------------------------------------------------------
    R = [0]*n+[eye(3)];     r = [0]*n;  ac = [0]*n
    w  = [zeros(3)]+[0]*n;  wp = [zeros(3)]+[0]*n;     vp = [-g0]+[0]*n
    tauF = [0]*n;   tauM = [0]*n;   M0 = [0]*n
    z0 = array([0,0,1])
    # Forward Recursion
    for i in range(n):
        A = dh(theta[i],d[i],a[i],alpha[i])     # Matrix (i-1)A(i)
        R[i]  = A[0:3,0:3]                      # Matrix (i-1)R(i)
        r[i]  = array([a[i], d[i]*sin(alpha[i]), d[i]*cos(alpha[i])])
        # Velocidades y aceleraciones
        w[i+1]   = R[i].T@(w[i] + qp[i]*z0)
        wp[i+1]  = R[i].T@(wp[i] + qpp[i]*z0) + cross(w[i], array([0,0,qp[i]]))
        vp[i+1]  = R[i].T@vp[i] + cross(wp[i+1], r[i]) + cross(w[i+1], cross(w[i+1], r[i]))
    
    for i in range(n):
        ac[i] = vp[i+1] + cross(wp[i+1], rc[i]) + cross(w[i+1], cross(w[i+1], rc[i]))
    
    # Backward Recursion
    for j in range(n):
        i = n - j - 1
        F[i] = R[i+1]@F[i+1] + m[i]*ac[i]
        M[i] = R[i+1]@(M[i+1] + cross(R[i+1].T@r[i],F[i+1])) + cross(r[i]+rc[i],m[i]*ac[i]) + I[i]@wp[i+1] + cross(w[i+1],I[i]@w[i+1])
        z = R[i][2,0:3]
        tauF[i] = F[i].T@z
        tauM[i] = M[i].T@z
        M0[i] = R[i]@M[i]   # Momentum at Base System {S0}
    Tau = array(tauM)
    return w,wp,vp,ac,M,M0,Tau

















