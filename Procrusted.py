# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:16:45 2023

@author: Romain Chauvaux
"""

"""
Compare two different solutions of the MDGP problem, plot them ang give their Frobenius error. 
Based on the work of Amaral et al.
"""

import numpy as np
import matplotlib.pyplot as plt

# Solution obtained by the algorithm
file1 = "Sol_1ptq.txt"

# Solution given by the Protein Data Bank
file2 = "1ptq.pdb"

# Type of the protein
Type= "Atom"
#Type= "HETATM"

fileCoo = open(file1,"r")


lignes = fileCoo.readlines()
ligne = lignes[0].split(" ")
N = int(ligne[0])
Xk = np.zeros((3,N));
for i in range(N):
    ligne = lignes[i+1].split(" ")
    Xk[0,i]= float(ligne[0]) 
    Xk[1,i]= float(ligne[1])
    Xk[2,i]= float(ligne[2])
fileCoo.close()


XStar = np.zeros((3,N))

file = open("file2","r")
lignes = file.readlines()
index = 0
while(True):
    ligne = lignes[index].split(" ")
    a = ligne[0]
    if ligne[0] == Type:
        break
    index+=1
for i in range(N):
    ligne = lignes[index+i].split(" ")
    ligne2 = [ele for ele in ligne if ele != '']
    XStar[0,i]= float(ligne2[6])
    XStar[1,i]= float(ligne2[7])
    XStar[2,i]= float(ligne2[8])
file.close()




J = np.identity(N)
J = J - 1/N*np.ones((N,N))

XkJ = Xk@J
XStarJ = XStar@J



C = XkJ@XStarJ.T

U, s, V = np.linalg.svd(C)
Q = V.T@U.T

Err_Frobenius = np.linalg.norm(Q@XkJ - XStarJ)**2


X = Q@XkJ
fig = plt.figure()
ax = fig.add_subplot(projection='3d')



one = ax.scatter(X[0,:],X[1,:],X[2,:], c='r', marker='o',s=10,label = "Xk")
two = ax.scatter(XStarJ[0,:],XStarJ[1,:],XStarJ[2,:], c='b', marker='v',s=10,label = "X*")


ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title("Comparison between X* and X for " + file2 + "protein", fontsize = 15)
plt.legend(["X*",'X'],loc='lower left',numpoints=1, ncol=3, fontsize=10, bbox_to_anchor=(0, 0))

plt.show()


Exi= np.zeros(N)
for i in range(N):
    Sub = X - XStarJ
    Norme = max(abs(Sub[:,i]))
    Exi[i] = Norme/max(1,max(abs(XStarJ[:,i])))
Ex = max(Exi)
    
# The error 
print("The Frobenius error is : {}".format(Err_Frobenius))