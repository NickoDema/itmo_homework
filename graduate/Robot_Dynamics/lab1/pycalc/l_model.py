#!/usr/bin/python3
"""
This document conatin math model of 3 joints plain manipulator
"""

import sympy as sp
import numpy as np
from sympy import cos, sin
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sp.init_printing()

# define symbols
q_1, q_2, q_3 = sp.symbols('q_1, q_2, q_3')
dq_1, dq_2, dq_3 = sp.symbols('\dot{q}_1, \dot{q}_2, \dot{q}_3')
ddq_1, ddq_2, ddq_3 = sp.symbols('\ddot{q}_1, \ddot{q}_2, \ddot{q}_3')
l_1, l_2, l_3 = sp.symbols('l_1, l_2, l_3')
l_1c, l_2c, l_3c = sp.symbols('l_1c, l_2c, l_3c')
m_1, m_2, m_3, g = sp.symbols('m_1, m_2, m_3, g')
I_1c, I_2c, I_3c = sp.symbols('I_1c, I_2c, I_3c')

# I[1] Kinetic energy
# See Robot Modelling and Control (Spong)
J_c1 = sp.Matrix([[-l_1c*sin(q_1), 0, 0], [l_1c*cos(q_1), 0, 0], [0, 0, 0]])
J_c2 = sp.Matrix([[-l_1*sin(q_1) - l_2c*sin(q_1 + q_2), -l_2c*sin(q_1 + q_2), 0],
    [l_1*cos(q_1) + l_2c*cos(q_1 + q_2), l_2c*cos(q_1 + q_2), 0],
    [0, 0, 0]])

jc3_11 = -l_1*sin(q_1) - l_2*sin(q_1 + q_2) - l_3c*sin(q_1 + q_2 + q_3)
jc3_12 = -l_2*sin(q_1 + q_2) - l_3c*sin(q_1 + q_2 + q_3)
jc3_13 = -l_3c*sin(q_1 + q_2 + q_3)
jc3_21 = l_1*cos(q_1) + l_2*cos(q_1 + q_2) + l_3c*cos(q_1 + q_2 + q_3)
jc3_22 = l_2*cos(q_1 + q_2) + l_3c*cos(q_1 + q_2 + q_3)
jc3_23 = l_3c*cos(q_1 + q_2 + q_3)
J_c3 = sp.Matrix([[jc3_11, jc3_12, jc3_13], [jc3_21, jc3_22, jc3_23], [0, 0, 0]])

# Calculate D(q) -- dq^T*D(q)*dq, D:
I1 = sp.Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])*I_1c
I2 = sp.Matrix([[1, 1, 0], [1, 1, 0], [0, 0, 0]])*I_2c
I3 = sp.ones(3)*I_3c
D = sp.simplify(m_1*J_c1.T*J_c1 + m_2*J_c2.T*J_c2 + m_3*J_c3.T*J_c3 + I1 + I2 + I3)

# Derivative of Kinetic energy -- Koriolis Matrix
# c_kij         kj/i                         ki/j                    ij/k
c111 =          1/2*(sp.diff(D[0, 0], q_1) + sp.diff(D[0, 0], q_1) - sp.diff(D[0, 0], q_1))
c112 = c121 =   1/2*(sp.diff(D[0, 1], q_1) + sp.diff(D[0, 0], q_2) - sp.diff(D[0, 1], q_1))
c113 = c131 =   1/2*(sp.diff(D[0, 2], q_1) + sp.diff(D[0, 0], q_3) - sp.diff(D[0, 2], q_1))
c122 =          1/2*(sp.diff(D[0, 1], q_2) + sp.diff(D[0, 1], q_2) - sp.diff(D[1, 1], q_1))
c123 = c132 =   1/2*(sp.diff(D[0, 2], q_2) + sp.diff(D[0, 1], q_3) - sp.diff(D[1, 2], q_1))
c133 =          1/2*(sp.diff(D[0, 2], q_3) + sp.diff(D[0, 2], q_3) - sp.diff(D[2, 2], q_1))

c211 =          1/2*(sp.diff(D[1, 0], q_1) + sp.diff(D[1, 0], q_1) - sp.diff(D[0, 0], q_2))
c212 = c221 =   1/2*(sp.diff(D[1, 1], q_1) + sp.diff(D[1, 0], q_2) - sp.diff(D[0, 1], q_2))
c213 = c231 =   1/2*(sp.diff(D[1, 2], q_1) + sp.diff(D[1, 0], q_3) - sp.diff(D[0, 2], q_2))
c222 =          1/2*(sp.diff(D[1, 1], q_2) + sp.diff(D[1, 1], q_2) - sp.diff(D[1, 1], q_2))
c223 = c232 =   1/2*(sp.diff(D[1, 2], q_2) + sp.diff(D[1, 1], q_3) - sp.diff(D[1, 2], q_2))
c233 =          1/2*(sp.diff(D[1, 2], q_3) + sp.diff(D[1, 2], q_3) - sp.diff(D[2, 2], q_2))

c311 =          1/2*(sp.diff(D[2, 0], q_1) + sp.diff(D[2, 0], q_1) - sp.diff(D[0, 0], q_3))
c312 = c321 =   1/2*(sp.diff(D[2, 1], q_1) + sp.diff(D[2, 0], q_2) - sp.diff(D[0, 1], q_3))
c313 = c331 =   1/2*(sp.diff(D[2, 2], q_1) + sp.diff(D[2, 0], q_3) - sp.diff(D[0, 2], q_3))
c322 =          1/2*(sp.diff(D[2, 1], q_2) + sp.diff(D[2, 1], q_2) - sp.diff(D[1, 1], q_3))
c323 = c332 =   1/2*(sp.diff(D[2, 2], q_2) + sp.diff(D[2, 1], q_3) - sp.diff(D[1, 2], q_3))
c333 =          1/2*(sp.diff(D[2, 2], q_3) + sp.diff(D[2, 2], q_3) - sp.diff(D[2, 2], q_3))

c11 = c111*dq_1 + c112*dq_2 + c113*dq_3
c12 = c121*dq_1 + c122*dq_2 + c123*dq_3
c13 = c131*dq_1 + c132*dq_2 + c133*dq_3

c21 = c211*dq_1 + c212*dq_2 + c213*dq_3
c22 = c221*dq_1 + c222*dq_2 + c223*dq_3
c23 = c231*dq_1 + c232*dq_2 + c233*dq_3

c31 = c311*dq_1 + c312*dq_2 + c313*dq_3
c32 = c321*dq_1 + c322*dq_2 + c323*dq_3
c33 = c331*dq_1 + c332*dq_2 + c333*dq_3

C = sp.Matrix([[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]])

# I[2] Potentional energy
P1 = m_1*g*l_1c*sin(q_1)
P2 = m_2*g*(l_1*sin(q_1) + l_2c*sin(q_1 + q_2))
P3 = m_3*g*(l_1*sin(q_1) + l_2*sin(q_1 + q_2) + l_3c*sin(q_1 + q_2 + q_3))
P = P1 + P2 + P3

# Derivative of Potentional energy
phi1 = sp.diff(P, q_1)
phi2 = sp.diff(P, q_2)
phi3 = sp.diff(P, q_3)

G = sp.Matrix([[phi1], [phi2], [phi3]])

def readFile(filename):
    file = open(filename, "r")
    data = np.array([range(29)])
    i = 1
    for line in file:
        if i == 8:
            j = 0
            for a in line.split('\t'):
                if len(a) > 1:
                    # print(str(j) + ": " + a)
                    j = j+1
        if (i > 8):
            words = np.fromstring(line, dtype = np.float, sep = '\t')
            if len(words) > 1:
                data = np.vstack((data, words))
        i = i + 1

    print("Done!")
    return data

def processData(raw_data):
    titles = np.array([('R-Hip Flex-Extension', 6), ('R-Knee Flex-Extension', 12), ('R-Ankle Dorsi-Plantarflex', 17)])
    # Convert degrees to radians
    q_1 = raw_data[1:, int(titles[0, 1])]*np.pi/180
    q_2 = raw_data[1:, int(titles[1, 1])]*np.pi/180
    q_3 = raw_data[1:, int(titles[2, 1])]*np.pi/180
    x = raw_data[1:, 0]

    # Interpolate data
    nx = np.linspace(0, 100, num=501, endpoint=True)
    tck = interpolate.splrep(x, q_1, s=0)
    nq_1 = interpolate.splev(nx, tck, der=0)

    tck = interpolate.splrep(x, q_2, s=0)
    nq_2 = interpolate.splev(nx, tck, der=0)

    tck = interpolate.splrep(x, q_3, s=0)
    nq_3 = interpolate.splev(nx, tck, der=0)

    # Numecical derivative
    dq_1 = np.diff(nq_1)/np.diff(nx)
    dq_2 = np.diff(nq_2)/np.diff(nx)
    dq_3 = np.diff(nq_3)/np.diff(nx)
    ddq_1 = np.diff(dq_1)/np.diff(nx[:-1])
    ddq_2 = np.diff(dq_2)/np.diff(nx[:-1])
    ddq_3 = np.diff(dq_3)/np.diff(nx[:-1])

    size = len(ddq_1)
    result = [nx[0:size], [nq_1[0:size], nq_2[0:size], nq_3[0:size]],
        [dq_1[0:size], dq_2[0:size], dq_3[0:size]],
        [ddq_1[0:size], ddq_2[0:size], ddq_3[0:size]]]

    print("Done!")
    return result

# l - link length, c - vector of center mass
# I - voctor of Moments of Inertia, m - mass of link
def getNumTorque(l, c, I, m, g):
    l_param = []
    c_param = []
    I_param = []
    m_param = []

    for i in range(len(l)):
        l_param.append(("l_" + str(i + 1), l[i]))
        c_param.append(("l_" + str(i + 1) + "c", c[i][0]))
        I_param.append(("I_" + str(i + 1) + "c", I[i]))
        m_param.append(("m_" + str(i + 1), m[i]))
    params = l_param + c_param + I_param + m_param + [("g", g)]

    numD = D.subs(params)
    numC = C.subs(params)
    numG = G.subs(params)

    dq = sp.Matrix([[dq_1], [dq_2], [dq_3]])
    ddq = sp.Matrix([[ddq_1], [ddq_2], [ddq_3]])
    torque = numD*ddq + numC*dq + numG

    print("Done!")
    return torque


def calcTorque(numTorque, data):
    q = data[1]
    dq = data[2]
    ddq = data[3]
    tau1 = []
    tau2 = []
    tau3 = []
    for i in range(len(data[0])):
        q_param = [("q_" + str(j + 1), q[j][i]) for j in range(3)]
        dq_param = [('\dot{q}_' + str(j + 1), dq[j][i]) for j in range(3)]
        ddq_param = [('\ddot{q}_' + str(j + 1), ddq[j][i]) for j in range(3)]
        params = q_param + dq_param + ddq_param
        if i == 3: print(params);

        tau = numTorque.subs(params)
        tau1.append(tau[0])
        tau2.append(tau[1])
        tau3.append(tau[2])

    return [tau1, tau2, tau3]
    print("Done!")

## Calculation
l = [0.300, 0.300, 0.190]
h = [0.20, 0.20, 0.100]
c = [[0.150, 0], [0.150, 0], [0.90, 0]]
m = [0.1, 0.1, 0.35]
g = 9.8
I = [m[i]*(l[i]**2 + h[i]**2)/12 for i in range(3)]

filename = "/home/senex/temp/ModellingData/Erofeev Mikhail (мужчина)/Erofeev_Mihhail_02_15_06_Angles.emt"
raw_data = readFile(filename)
data = processData(raw_data)
numTorque = getNumTorque(l, c, I, m, g)
res = calcTorque(numTorque, data)

## Plot results
pp = PdfPages('multipage.pdf')
plt.figure(figsize=(20, 15))
for i in range(3):
    ## Angles
    plt.subplot(4, 3, (i + 1))
    plt.plot(data[0], data[1][i], 'r')
    plt.xlabel("x")
    plt.ylabel("$q_" + str(i + 1) + "$" + " [rad]")
    plt.grid(True)
    plt.xlim(0, 100)

    ## Angular velocity
    plt.subplot(4, 3, (i + 4))
    plt.plot(data[0], data[2][i], 'r')
    plt.xlabel("x")
    plt.ylabel("$\dot{q}_" + str(i + 1) + "$" + " [rad/s]")
    plt.grid(True)
    plt.xlim(0, 100)

    ## Angular acceleration
    plt.subplot(4, 3, (i + 7))
    plt.plot(data[0], data[3][i], 'r')
    plt.xlabel("x")
    plt.ylabel("$\ddot{q}_" + str(i + 1) + "$" + " [rad/s$^2$]")
    plt.grid(True)
    plt.xlim(0, 100)

    ## Torque
    plt.subplot(4, 3, (i + 10))
    plt.plot(data[0], res[i], 'r')
    plt.xlabel("x")
    plt.ylabel("$M_" + str(i + 1) + "$ [N]")
    plt.grid(True)
    plt.xlim(0, 100)

plt.tight_layout()
plt.savefig(pp, format='pdf')
plt.show()

pp.close()
