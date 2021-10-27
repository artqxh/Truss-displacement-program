import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.linalg import matrix_power

elementnumbercheck = 1

elements = np.loadtxt('points.txt')

maxvalue = 0
for element in elements:
    if maxvalue < element[5]:
        maxvalue = element[5]

    elif maxvalue < element[6]:
        maxvalue = element[6]

#print(maxvalue)

lelements = len(elements)


'''
p1= (1,0,0)
p2= (2,0,2)
p3= (3,2,2)
p4= (4,2,0)
p5= (5,4,0)
'''


e1 = (0, 0, 0, 2, 2)
e2 = (0, 0, 2, 0, 2)
e3 = (0, 2, 2, 0, 2.82)
e4 = (0, 2, 2, 2, 2)
e5 = (2, 2, 2, 0, 2)
e6 = (2, 0, 4, 0, 2)
e7 = (2, 2, 4, 0, 2.82)


def ssin(x):
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    ssin = (x[4]-x[2])/l
    return ssin

def ccos(x):
    l = math.sqrt((x[3] - x[1]) ** 2 + (x[4] - x[2]) ** 2)
    ccos = (x[3]-x[1])/l

    return ccos

def DC(ccos, ssin):
    DC = np.array([[ccos, ssin, 0, 0], [(-ssin), ccos, 0, 0], [0, 0, ccos, ssin], [0, 0, (-ssin), ccos]])

    return DC

DC1 = DC(ccos(elements[elementnumbercheck]), ssin(elements[elementnumbercheck]))
#print(DC1)
#print(DC1.T)
'''
def div(matrixDC):
    m00 = matrixDC[0][0]
    m01 = matrixDC[0][1]
    m10 = matrixDC[1][0]
    m11 = matrixDC[1][1]

    ARR = np.array([[m00,m01],[m10,m11]])
    return ARR
'''
'''
for i in range(lelements):
    print(ssin(elements[i]))
    print(ccos(elements[i]))
i=0
'''

#print(div(DC1))


for i in range(lelements):
    blank = []
    plt.plot([elements[i][1], elements[i][3]], [elements[i][2], elements[i][4]], 'b')
    i = i+1

plt.title('Kratownica')
plt.xlabel('Dlugosc [m]')
plt.ylabel('Dlugosc [m]')
plt.show()


def matrixk():
    matrixi = np.array([[1, 0, -1, 0],
                        [0, 0, 0, 0],
                        [-1, 0, 1, 0],
                        [0, 0, 0, 0]])
    k = matrixi
    return k

#print(matrixk(elements[elementnumbercheck]))

def matrixkprim(DC, k, x):
    E = 2 * (10 ** 11)
    A = 3.14 * (10 ** (-4))
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    matrixkprim = ((A*E)/l)*(np.matmul(np.matmul(DC.T, k), DC))
    return matrixkprim

k1 = matrixk()

#print(matrixkprim(DC1, k1, elements[elementnumbercheck]))


def spr(x, ssin, ccos):
    E = 2 * (10 ** 11)
    A = 3.14 * (10 ** (-4))
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    spr = ((A*E)/l)*np.array([[ccos*ccos, ssin*ccos, -ccos*ccos, -ssin*ccos],
                               [ssin*ccos, ssin*ssin, -ssin*ccos, -ssin*ssin],
                               [-ccos*ccos, -ssin*ccos, ccos*ccos, ssin*ccos],
                               [-ssin*ccos, -ssin*ssin, ssin*ccos, ssin*ssin]])
    return spr

#print('\n')
#print(spr(elements[elementnumbercheck], ssin(elements[elementnumbercheck]), ccos(elements[elementnumbercheck])))

Z1 = np.zeros((2*int(maxvalue), 2*int(maxvalue)))

for element in elements:
    Z = np.zeros((2*int(maxvalue), 2*int(maxvalue)))
    M = matrixkprim((DC(ccos(element), ssin(element))), k1, element)
    #print(M)

    m = int(element[5])
    n = int(element[6])
    #print(m, n)

    m = 2 * m - 1
    n = 2 * n - 1
    Z[m - 1][m - 1] = M[0][0]
    Z[m - 1][m] = M[1][0]
    Z[m][m - 1] = M[0][1]
    Z[m][m] = M[1][1]

    Z[n - 1][n - 1] = M[2][2]
    Z[n - 1][n] = M[3][2]
    Z[n][n - 1] = M[2][3]
    Z[n][n] = M[3][3]

    Z[m - 1][n - 1] = M[0][2]
    Z[m - 1][n] = M[1][2]
    Z[m][n - 1] = M[0][3]
    Z[m][n] = M[1][3]

    Z[n - 1][m - 1] = M[2][0]
    Z[n - 1][m] = M[3][0]
    Z[n][m - 1] = M[2][1]
    Z[n][m] = M[3][1]
    Z1 = Z1 + Z
print(Z1)

F = np.zeros((10, 1))
F[8][0] = 100000
F[9][0] = -1000
print(F)
Z1[:][0] = 0
Z1[:][1] = 0
Z1[:][2] = 0
#g = 1
#h = 2

#g = element[g]
#h = element[h]
#print(g, h)
for zzz in range(len(Z1)):
    Z1[zzz][0] = 0
    Z1[zzz][1] = 0
    Z1[zzz][2] = 0
print(Z1)
Z1[0][0] = 1
Z1[1][1] = 1
Z1[2][2] = 1

#np.set_printoptions(formatter={'float': "{0:0.0e}".format})

print(Z1)
U = np.linalg.pinv(Z1, -1)
#U2 = np.linalg.matrix_power(Z1, -1)
U = np.matmul(np.linalg.pinv(Z1, -1), F)
print(U)

i = 0


for element in elements:
    j = int(element[5])
    j = 2*j -1
    elements[i][1] = elements[i][1] + U[j - 1][0]
    print(elements[i][1], U[j - 1][0])
    elements[i][2] = elements[i][2] + U[j][0]
    k = int(element[6])
    k = 2*k -1
    elements[i][3] = elements[i][1] + U[k - 1][0]
    elements[i][4] = elements[i][2] + U[k][0]
    i = i +1

print(elements)


for i in range(lelements):
    blank = []
    plt.plot([elements[i][1], elements[i][3]], [elements[i][2], elements[i][4]], 'b')
    i = i+1

plt.title('Kratownica')
plt.xlabel('Dlugosc [m]')
plt.ylabel('Dlugosc [m]')
plt.show()
