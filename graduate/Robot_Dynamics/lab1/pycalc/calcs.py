#!/usr/bin/env python3

import math as m
import sympy as sp
from sympy import sin, cos
sp.init_printing()

x, y, z, L2, L3, d1, th2, th3, pi, l2, l3, m1, m2, m3, g, r2, r3 = sp.symbols('x y z L2 L3 d1 theta2 theta3 pi l2 l3 m1 m2 m3 g r2 r3')
a, al, d, th, d, dt2, dt3 = sp.symbols('a al d th \dot{d} \dot{theta_2} \dot{theta_3}')

# l,q1,q2,r = sp.symbols('l q1 q2 r')
# W = -l((q1+q2)**2)-((q1+q2)**2)*(-l)
# sp.simplify(W)

H = sp.Matrix([[cos(th), -sin(th)*cos(al), sin(th)*sin(al), a*cos(th)],
     [sin(th), cos(th)*cos(al), -cos(th)*sin(al), a*sin(th)],
     [0, sin(al), cos(al), d],
     [0, 0, 0, 1]])

H1 = H.subs([(a,0),(al,-m.pi/2),(th,0)])
H2 = H.subs([(a,L2),(al,0),(d,0),(th,th2)])
H2c= H.subs([(a,l2),(al,0),(d,0),(th,th2)])
H3 = H.subs([(a,L3),(al,m.pi/2),(d,0),(th,th3+pi/2)])
H3c= H.subs([(a,l3),(al,0),(d,0),(th,th3)])

#sp.latex(sp.N(H3,3))
print(sp.latex(sp.N(H1,2)))
sp.N(H1,2)
print(sp.latex(sp.N(H2,2)))
print(sp.latex(sp.N(H3,2)))
print(sp.latex(H2c))
print(sp.latex(H3c))




#------------------------------------------------------
r1c1 = sp.Matrix([0, 0, 0, 1])
r2c2 = sp.Matrix([0, 0, 0, 1])
r3c3 = sp.Matrix([0, 0, 0, 1])
r0c1=H1*r1c1
r0c1
print(sp.latex(H1*r1c1))
r0c2=(H1*H2c*r2c2).evalf(3)
r0c2
print(sp.latex(H1*H2*r2c2.evalf(3)))
r0c3=(H1*H2*H3c*r3c3).evalf(3)
sp.simplify(r0c3)
print(sp.latex(r0c3))
r0c1.row_del(3)
r0c1
r0c2.row_del(3)
r0c2[1] = 0
sp.simplify(r0c2)
r0c3.row_del(3)
r0c3
r0c3[1] = 0
r0c3


gv = sp.Matrix([-g, 0, 0])
U = -gv.T*(m1*r0c1 + m2*r0c2 + m3*r0c3)
U = sp.simplify(sp.expand(U))
print(sp.latex(U))
print(sp.latex(sp.sympify(U)))

sp.diff(U,q[2])

# Kinetik energy calc
I2 = sp.Matrix([[m2*r2*r2/2,    0,                          0],
                [0,             m2*r2*r2/4 + m2*L2*L2/12,   0],
                [0,             0,                          m2*r2*r2/4 + m2*L2*L2/12]])
I3 = sp.Matrix([[m3*r3*r3/2, 0, 0],
                [0, m3*r3*r3/4 + m3*L3*L3/12, 0],
                [0, 0, m3*r3*r3/4 + m3*L3*L3/12]])
Ju1 = sp.Matrix([[0,0,0],
                 [0,0,0],
                 [1,0,0]])
Ju2 = sp.Matrix([[0,    -l2*sin(th2),   0],
                 [0,    0,              0],
                 [1,    l2*cos(th2),    0]])
Ju3 = sp.Matrix([[0,    -L2*sin(th2)-l3*sin(th2+th3),    -l3*sin(th2+th3)],
                 [0,    0,                              0],
                 [1,    L2*cos(th2)+l3*sin(th2+th3),    l3*cos(th2+th3)]])
Jw1 = sp.Matrix([[0,0,0],[0,0,0],[0,0,0]])
Jw2 = sp.Matrix([[0,0,0],[0,1,0],[0,0,0]])
Jw3 = sp.Matrix([[0,0,0],[0,1,1],[0,0,0]])

R1 = sp.Matrix(H1)
R1.row_del(3)
R1.col_del(3)
R1[1,1] = 0
R1[2,2] = 0
R2 = sp.Matrix(H2c)
R2.row_del(3)
R2.col_del(3)
R3 = sp.Matrix(H3c)
R3.row_del(3)
R3.col_del(3)
#R3[0,1] = 0
#R3[1,1] = 0
#R3[2,2] = 0
R1
R2
R3
#
#R3 = sp.Matrix([[cos(th3), -sin(th3), 0],[sin(th3), cos(th3), 0],[0, 0, 1]])

qp = sp.Matrix([[d],[dt2],[dt3]])

q = sp.Matrix([d1, th2, th3])

R02 = R1*R2
R02
R03 = R02*R3
R03
K = 1/2*qp.T*(m1*Ju1.T*Ju1 + m2*Ju2.T*Ju2 + m3*Ju3.T*Ju3 + Jw2.T*R02*I2*R02.T*Jw2 + Jw3.T*R03*I3*R03.T*Jw3)*qp
sp.simplify(sp.expand(K))
print(sp.latex(sp.simplify(sp.expand(K))))
D = sp.simplify(sp.expand(m1*Ju1.T*Ju1 + m2*Ju2.T*Ju2 + m3*Ju3.T*Ju3 + Jw2.T*R02*I2*R02.T*Jw2 + Jw3.T*R03*I3*R03.T*Jw3))
D = E.evalf(4)
D
print(sp.latex(D))
K
A = sp.simplify(K)
A
print(sp.latex(A))

# Koriolis Matrix
#kij            kj/i                         ki/j                    ij/k
c111 =          1/2*(sp.diff(D[0, 0], q[0]) + sp.diff(D[0, 0], q[0]) - sp.diff(D[0, 0], q[0]))
c112 = c121 =   1/2*(sp.diff(D[0, 1], q[0]) + sp.diff(D[0, 0], q[1]) - sp.diff(D[0, 1], q[0]))
c113 = c131 =   1/2*(sp.diff(D[0, 2], q[0]) + sp.diff(D[0, 0], q[2]) - sp.diff(D[0, 2], q[0]))
c122 =          1/2*(sp.diff(D[0, 1], q[1]) + sp.diff(D[0, 1], q[1]) - sp.diff(D[1, 1], q[0]))
c123 = c132 =   1/2*(sp.diff(D[0, 2], q[1]) + sp.diff(D[0, 1], q[2]) - sp.diff(D[1, 2], q[0]))
c133 =          1/2*(sp.diff(D[0, 2], q[2]) + sp.diff(D[0, 2], q[2]) - sp.diff(D[2, 2], q[0]))

c211 =          1/2*(sp.diff(D[1, 0], q[0]) + sp.diff(D[1, 0], q[0]) - sp.diff(D[0, 0], q[1]))
c212 = c221 =   1/2*(sp.diff(D[1, 1], q[0]) + sp.diff(D[1, 0], q[1]) - sp.diff(D[0, 1], q[1]))
c213 = c231 =   1/2*(sp.diff(D[1, 2], q[0]) + sp.diff(D[1, 0], q[2]) - sp.diff(D[0, 2], q[1]))
c222 =          1/2*(sp.diff(D[1, 1], q[1]) + sp.diff(D[1, 1], q[1]) - sp.diff(D[1, 1], q[1]))
c223 = c232 =   1/2*(sp.diff(D[1, 2], q[1]) + sp.diff(D[1, 1], q[2]) - sp.diff(D[1, 2], q[1]))
c233 =          1/2*(sp.diff(D[1, 2], q[2]) + sp.diff(D[1, 2], q[2]) - sp.diff(D[2, 2], q[1]))

c311 =          1/2*(sp.diff(D[2, 0], q[0]) + sp.diff(D[2, 0], q[0]) - sp.diff(D[0, 0], q[2]))
c312 = c321 =   1/2*(sp.diff(D[2, 1], q[0]) + sp.diff(D[2, 0], q[1]) - sp.diff(D[0, 1], q[2]))
c313 = c331 =   1/2*(sp.diff(D[2, 2], q[0]) + sp.diff(D[2, 0], q[2]) - sp.diff(D[0, 2], q[2]))
c322 =          1/2*(sp.diff(D[2, 1], q[1]) + sp.diff(D[2, 1], q[1]) - sp.diff(D[1, 1], q[2]))
c323 = c332 =   1/2*(sp.diff(D[2, 2], q[1]) + sp.diff(D[2, 1], q[2]) - sp.diff(D[1, 2], q[2]))
c333 =          1/2*(sp.diff(D[2, 2], q[2]) + sp.diff(D[2, 2], q[2]) - sp.diff(D[2, 2], q[2]))
#ijk                 sp.diff(D[k, j], q[i]) + sp.diff(D[k, i], q[j]) - sp.diff(D[i, j], q[k])

c122
c123
c133
c213
c222
c223
c233
c312
c322


#ki?
#ji   ckij         ckij         ckij
c11 = c111*qp[0] + c112*qp[1] + c113*qp[2]
c12 = c121*qp[0] + c122*qp[1] + c123*qp[2]
c13 = c131*qp[0] + c132*qp[1] + c133*qp[2]

c21 = c211*qp[0] + c212*qp[1] + c213*qp[2]
c22 = c221*qp[0] + c222*qp[1] + c223*qp[2]
c23 = c231*qp[0] + c232*qp[1] + c233*qp[2]

c31 = c311*qp[0] + c312*qp[1] + c313*qp[2]
c32 = c321*qp[0] + c322*qp[1] + c323*qp[2]
c33 = c331*qp[0] + c332*qp[1] + c333*qp[2]

#                                               kij
print(sp.latex(sp.simplify(sp.expand(c111)))) #c111 =
print(sp.latex(sp.simplify(sp.expand(c112)))) #c112 = c121
print(sp.latex(sp.simplify(sp.expand(c113)))) #c113 = c131
print(sp.latex(sp.simplify(sp.expand(c122)))) #c122 =
print(sp.latex(sp.simplify(sp.expand(c123)))) #c123 = c132
print(sp.latex(sp.simplify(sp.expand(c133)))) #c133 =
################################################################################
print(sp.latex(sp.simplify(sp.expand(c211)))) #c211 =
print(sp.latex(sp.simplify(sp.expand(c212)))) #c212 = c221
print(sp.latex(sp.simplify(sp.expand(c213)))) #c213 = c231
print(sp.latex(sp.simplify(sp.expand(c222)))) #c222 =
print(sp.latex(sp.simplify(sp.expand(c223)))) #c223 = c232
print(sp.latex(sp.simplify(sp.expand(c233)))) #c233 =
#####################################################################################
print(sp.latex(sp.simplify(sp.expand(c311)))) #c311 =
print(sp.latex(sp.simplify(sp.expand(c312)))) #c312 = c321
print(sp.latex(sp.simplify(sp.expand(c313)))) #c313 = c331
print(sp.latex(sp.simplify(sp.expand(c322)))) #c322 =
print(sp.latex(sp.simplify(sp.expand(c323)))) #c323 = c332
print(sp.latex(sp.simplify(sp.expand(c333)))) #c333 =


print(sp.latex(sp.simplify(sp.expand(c11))))
print(sp.latex(sp.simplify(sp.expand(c12))))
print(sp.latex(sp.simplify(sp.expand(c13))))
print(sp.latex(sp.simplify(sp.expand(c21))))
print(sp.latex(sp.simplify(sp.expand(c22))))
print(sp.latex(sp.simplify(sp.expand(c23))))
print(sp.latex(sp.simplify(sp.expand(c31))))
print(sp.latex(sp.simplify(sp.expand(c32))))
print(sp.latex(sp.simplify(sp.expand(c33))))


C = sp.Matrix([[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]])
C

# (sp.diff(D[k,j], q[i]) + sp.diff(D[k,i], q[j]) - sp.diff(D[i,j], q[k]))/2
# (sp.diff(D[2,2], q[1]) + sp.diff(D[2,1], q[2]) - sp.diff(D[1,2], q[2]))/2
#
# D.shape[0]-1
#
# kq = 3
# Ck = sp.MutableDenseNDimArray([0] * 27 ,(3,3,3))
# Ck
# for k in range(0, 2):
#     for i in range(0, D.shape[0]):
#         for j in range(0, D.shape[1]):
#             Ck[k,i,j] = (sp.diff(D[k,j], q[i]) + sp.diff(D[k,i], q[j]) - sp.diff(D[i,j], q[k]))/2
# #Ck[1,1,1] = (sp.diff(E[1,1], q[1]))
# #sp.simplify(Ck[1,1,1])
# Ck
#
# C = sp.MutableDenseNDimArray([0]*(3*3), (3,3))
# for k in range(0, kq-1):
#     for j in range(0, D.shape[0]-1):
#         for i in range(0, D.shape[1]-1):
#             C[k,j] = C[k,j] + Ck[k,j,i]*qp[i]
# C







R = sp.Matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
p = sp.Matrix([0, 0, 1])
R*p

u1, u2, u3 = sp.symbols('u1:4')
v1, v2, v3 = sp.symbols('v1:4')
P = (sp.Matrix([ u1 ,u2 , u3 ])).cross (sp.Matrix ([ v1 ,v2 , v3 ]))
P.subs([(u1,0),(u2,-1),(u3,0),(v1,5),(v2,0),(v3,2),])
