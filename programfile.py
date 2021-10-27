import numpy as np
import math
import matplotlib.pyplot as plt

elementnumbercheck = 2

file = open('points.txt', 'r')
allfile = file.readlines()


elements = []
#print(matrixi)


for line in allfile:
    preelements = []
    line = str(line)
    line = line.replace(" ", "")
    line = line.replace("\n", "")
    line = line[:9]
    for i in line:
        i = float(i)
        preelements.append(i)
    elements.append(preelements)

#print(elements[1][0])


lelements = len(elements)


'''
p1= (1,0,0)
p2= (2,0,2)
p3= (3,2,2)
p4= (4,2,0)
p5= (5,4,0)
'''


e1 = (0,0,0,2,2)
e2 = (0,0,2,0,2)
e3 = (0,2,2,0,2.82)
e4 = (0,2,2,2,2)
e5 = (2,2,2,0,2)
e6 = (2,0,4,0,2)
e7 = (2,2,4,0,2.82)


def ssin(x):
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    ssin = (x[4]-x[2])/l
    return ssin

def ccos(x):
    l = math.sqrt((x[3] - x[1]) ** 2 + (x[4] - x[2]) ** 2)
    ccos = (x[3]-x[1])/l

    return ccos

def DC(ccos,ssin):
    DC = np.array([[ccos,ssin,0,0],[(-ssin),ccos,0,0],[0,0,ccos,ssin],[0,0,(-ssin),ccos]])

    return DC

DC1 = DC(ccos(elements[elementnumbercheck]),ssin(elements[elementnumbercheck]))
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

'''
for i in range(7):
    plt.plot([elements[i][1], elements[i][3]],[elements[i][2], elements[i][4]],'b')
    i = i+1
plt.title('Kratownica')
plt.xlabel('Dlugosc [m]')
plt.ylabel('Dlugosc [m]')
plt.show()
'''

def matrixk():
    matrixi = np.array([[1, 0, -1, 0],
                        [0, 0, 0, 0],
                        [-1, 0, 1, 0],
                        [0, 0, 0, 0]])
    k = matrixi
    return k

#print(matrixk(elements[elementnumbercheck]))

def matrixkprim(DC,k,x):
    E = 2 * (10 ** (11))
    A = 3.14 * (10 ** (-4))
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    matrixkprim = ((A*E)/l)*(np.matmul(np.matmul(np.transpose(DC),k), DC))
    return matrixkprim

k1 = matrixk()

print(matrixkprim(DC1,k1,elements[elementnumbercheck]))


def spr(x,ssin,ccos):
    E = 2 * (10 ** (11))
    A = 3.14 * (10 ** (-4))
    l = math.sqrt((x[3]-x[1])**2+(x[4]-x[2])**2)
    spr = ((A*E)/l)*np.array([[ccos*ccos,ssin*ccos,-ccos*ccos,-ssin*ccos],
                               [ssin*ccos,ssin*ssin,-ssin*ccos,-ssin*ssin],
                               [-ccos*ccos,-ssin*ccos,ccos*ccos,ssin*ccos],
                               [-ssin*ccos,-ssin*ssin,ssin*ccos,ssin*ssin]])
    return spr

print('\n')
print(spr(elements[elementnumbercheck],ssin(elements[elementnumbercheck]),ccos(elements[elementnumbercheck])))

#k = np.zeros([lelements,lelements])

#print(k)
