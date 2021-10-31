import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.linalg as la

#Material parameters
global E, A, ro
A = 8 * (10 ** (-6))
E = 2 * (10 ** 11)
ro = 7860

#Files
pointsfile = 'points2.txt'
forcefile = 'force2.txt'
dirichletsfile = 'xyz2.txt'

#Set print options

np.set_printoptions(formatter={'float': "{0:0.0e}".format})

#Import from file

elements = np.loadtxt(pointsfile)

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
    plt.annotate(elements[i][5], (elements[i][1], elements[i][2]))
    plt.annotate(elements[i][6], (elements[i][3], elements[i][4]))
    i = i+1

#Stiffness matrix

#Statics

def matrixk():
    k = np.array([[1, 0, -1, 0],
                  [0, 0, 0, 0],
                  [-1, 0, 1, 0],
                  [0, 0, 0, 0]])
    return k


def matrixkprim(DC, k, x):
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    matrixkprim = ((A*E)/l)*(np.matmul(np.matmul(DC.T, k), DC))
    return matrixkprim

k1 = matrixk()

#Dynamics

def matrixm():
    m = np.array([[2,0,1,0],
                 [0,2,0,1],
                 [1,0,2,0],
                 [0,1,0,2]])
    return m

def matrixmprim(DC, m, x):
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    matrixmprim = ((ro*A*l)/6)*(np.matmul(np.matmul(DC.T, m), DC))
    return matrixmprim

m1 = matrixm()

#Stiffness matrix verification

def spr(x, ssin, ccos):
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    spr = ((A*E)/l)*np.array([[ccos*ccos, ssin*ccos, -ccos*ccos, -ssin*ccos],
                               [ssin*ccos, ssin*ssin, -ssin*ccos, -ssin*ssin],
                               [-ccos*ccos, -ssin*ccos, ccos*ccos, ssin*ccos],
                               [-ssin*ccos, -ssin*ssin, ssin*ccos, ssin*ssin]])
    return spr

#Inertia matrix after agregation
#Z1 - inertia matrix - statics
#M_ini - inertia matrix - dynamics

Z1 = np.zeros((2*int(maxvalue), 2*int(maxvalue)))
M_ini = np.zeros((2*int(maxvalue), 2*int(maxvalue)))

for element in elements:
    Z = np.zeros((2*int(maxvalue), 2*int(maxvalue)))
    M = matrixkprim((DC(ccos(element), ssin(element))), k1, element)


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

    # Dynamics ---------------------------------

    Z_P2 = np.zeros((2 * int(maxvalue), 2 * int(maxvalue)))
    M_P2 = matrixmprim((DC(ccos(element), ssin(element))), m1, element)

    Z_P2[m - 1][m - 1] = M_P2[0][0]
    Z_P2[m - 1][m] = M_P2[1][0]
    Z_P2[m][m - 1] = M_P2[0][1]
    Z_P2[m][m] = M_P2[1][1]

    Z_P2[n - 1][n - 1] = M_P2[2][2]
    Z_P2[n - 1][n] = M_P2[3][2]
    Z_P2[n][n - 1] = M_P2[2][3]
    Z_P2[n][n] = M_P2[3][3]

    Z_P2[m - 1][n - 1] = M_P2[0][2]
    Z_P2[m - 1][n] = M_P2[1][2]
    Z_P2[m][n - 1] = M_P2[0][3]
    Z_P2[m][n] = M_P2[1][3]

    Z_P2[n - 1][m - 1] = M_P2[2][0]
    Z_P2[n - 1][m] = M_P2[3][0]
    Z_P2[n][m - 1] = M_P2[2][1]
    Z_P2[n][m] = M_P2[3][1]
    M_ini = M_ini + Z_P2


#Force vector - length have to be 2*number of points

#Force manually
'''
F = np.zeros((10, 1))
F[8][0] = 1000000
F[9][0] = -1000000
'''

#Import force from file

F = np.loadtxt(forcefile)
length_F = len(F)
F = F.reshape(length_F, 1)

#Dirichlet's boundary conditions from file


P = np.loadtxt(dirichletsfile, dtype='str')

i = 0
L = []

for i in range(len(P)):
    #print(P[i][0])
    if P[i][1] == 'n':
        L.append(2*int(P[i][0])-2)
        L.append(2*int(P[i][0])-1)
    elif P[i][1] == 'p' and P[i][2] == 'x':
        L.append(2 * int(P[i][0]) - 2)
    elif P[i][1] == 'p' and P[i][2] == 'y':
        L.append(2 * int(P[i][0]) - 1)


i = 0
f = 0

for i in range(len(Z1)):
    for f in L:
        Z1[:][f] = 0
        Z1[i][f] = 0
        Z1[f][f] = 1

        #Dynamics

        M_ini[:][f] = 0
        M_ini[i][f] = 0
        M_ini[f][f] = 1


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
#plt.show()

#Creating eigenvector and eigenmatrix

w, u = la.eig(Z1, M_ini)

#Deleting values connected with dirichlet's conditions

w = np.delete(w, np.where(w == 1))

#Vibration frequency

f = 1/(2*np.pi)*np.sqrt(w)
f = np.sort(f)

print(f)


