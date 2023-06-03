# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 12:10:22 2023

@author: Romain Chauvaux
"""

"""
Create the data file describing the Moré and Wu problems 
"""
import numpy as np

# The size of the problems
s = 4
# The name of the output file
name = "Moré_and_Wu_s4.txt"



Coordinates = np.zeros((s**3,3))
Distance_Matrix = np.zeros((s**3,s**3))

for i in range(s):
    for j in range(s):
        for k in range(s):
            index = i + j*s + k*s**2
            Coordinates[index][0] = i
            Coordinates[index][1] = j
            Coordinates[index][2] = k

for i in range(s**3):
    for j in range(s**3):
        if i != j:
            Distance_Matrix[i][j] = np.linalg.norm(Coordinates[i] - Coordinates[j])
            
Practical_Distance = np.zeros((s**3*s**3,3))
count = 0
for i in range(s**3):
    for j in range(i+1,s**3,1):
        if Distance_Matrix[i][j] <= 2:
            Practical_Distance[count][0] = i
            Practical_Distance[count][1] = j
            Practical_Distance[count][2] = Distance_Matrix[i][j]
            count += 1


fichier = open(name, "w")
str_n = "{}\n".format(s**3)
str_count = "{}\n".format(count)
fichier.write(str_n)
fichier.write(str_count)
fichier.write("{}\n".format(0))
fichier.write("{}\n".format(2))
for i in range(count):
    txt= "{} {} {}\n".format(Practical_Distance[i][0]+1,Practical_Distance[i][1]+1,Practical_Distance[i][2])
    fichier.write(txt)
fichier.close()
            
            
        