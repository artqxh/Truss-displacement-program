import numpy as np
import math
import matplotlib.pyplot as plt

#Set print options

np.set_printoptions(formatter={'float': "{0:0.0e}".format})

#Import from file

elements = np.loadtxt('points.txt')

#Looking for points number

maxvalue = 0
for element in elements:
    if maxvalue < element[5]:
        maxvalue = element[5]

    elif maxvalue < element[6]:
        maxvalue = element[6]

#Looking for elements number

lelements = len(elements)

#Cos, sin calculation

def ssin(x):
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    ssin = (x[4]-x[2])/l
    return ssin

def ccos(x):
    l = math.sqrt((x[3] - x[1]) ** 2 + (x[4] - x[2]) ** 2)
    ccos = (x[3]-x[1])/l
    return ccos

#Cosinus direction matrix

def DC(ccos, ssin):
    DC = np.array([[ccos, ssin, 0, 0], [(-ssin), ccos, 0, 0], [0, 0, ccos, ssin], [0, 0, (-ssin), ccos]])
    return DC

#Drawing the main truss

for i in range(lelements):
    blank = []
    plt.plot([elements[i][1], elements[i][3]], [elements[i][2], elements[i][4]], 'b')
    i = i+1

#Stiffness matrix

def matrixk():
    k = np.array([[1, 0, -1, 0],
                  [0, 0, 0, 0],
                  [-1, 0, 1, 0],
                  [0, 0, 0, 0]])
    return k


def matrixkprim(DC, k, x):
    E = 2 * (10 ** 11)
    A = 3.14 * (10 ** (-4))
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    matrixkprim = ((A*E)/l)*(np.matmul(np.matmul(DC.T, k), DC))
    return matrixkprim

k1 = matrixk()

#Stiffness matrix verification

def spr(x, ssin, ccos):
    E = 2 * (10 ** 11)
    A = 3.14 * (10 ** (-4))
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    spr = ((A*E)/l)*np.array([[ccos*ccos, ssin*ccos, -ccos*ccos, -ssin*ccos],
                               [ssin*ccos, ssin*ssin, -ssin*ccos, -ssin*ssin],
                               [-ccos*ccos, -ssin*ccos, ccos*ccos, ssin*ccos],
                               [-ssin*ccos, -ssin*ssin, ssin*ccos, ssin*ssin]])
    return spr

#Stiffness matrix after agregation

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

#Force vector - length have to be 2*number of points

F = np.zeros((10, 1))
F[8][0] = 1000000
F[9][0] = -1000000

#Dirichlet's boundary conditions

Z1[:][0] = 0
Z1[:][1] = 0
Z1[:][2] = 0

i = 0

for i in range(len(Z1)):
    Z1[i][0] = 0
    Z1[i][1] = 0
    Z1[i][2] = 0

#print(Z1)
Z1[0][0] = 1
Z1[1][1] = 1
Z1[2][2] = 1

#Calculating points displacements

U = np.linalg.pinv(Z1, -1)
#U = np.linalg.matrix_power(Z1, -1)
U = np.matmul(np.linalg.pinv(Z1, -1), F)

#Creating new points after displacements

i = 0
for element in elements:
    j = int(element[5])
    j = 2*j -1
    elements[i][1] = elements[i][1] + U[j - 1][0]
    #print(elements[i][1], U[j - 1][0])
    elements[i][2] = elements[i][2] + U[j][0]
    k = int(element[6])
    k = 2*k -1
    elements[i][3] = elements[i][3] + U[k - 1][0]
    elements[i][4] = elements[i][4] + U[k][0]
    i = i +1

#Drawing truss after displacements

i = 0
for i in range(lelements):
    plt.plot([elements[i][1], elements[i][3]], [elements[i][2], elements[i][4]], 'r')
    i = i+1

plt.title('Truss')
plt.xlabel('Length [m]')
plt.ylabel('Length [m]')
plt.show()
