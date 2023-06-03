# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 10:47:55 2023

@author: HP
"""



fileCoo = open("Sol_1ao2.txt","r")
lignes = fileCoo.readlines()
ligne = lignes[0].split(" ")
N = int(ligne[0])
#New_Coo = np.zeros((N,3));
#for i in range(N):
#    ligne = lignes[i].split(" ")
#    New_Coo[i,0]= float(ligne[0])
#    New_Coo[i,1]= float(ligne[1])
#    New_Coo[i,2]= float(ligne[2])

fileBefore = open("Brut/1ao2.pdb","r")

fileAfter= open("1ao2_changed.pdb","w")

lignes_B = fileBefore.readlines()
index = 0
while(True):
    ligne_b = lignes_B[index].split(" ")
    #if ligne_b[0] == 'ATOM':
    if ligne_b[0] == 'HETATM':
        break
    else:
        fileAfter.write(lignes_B[index])
    index+=1
    
for i in range(N):
    ligne_b = lignes_B[index+i].split(" ")
    ligne = lignes[i+1].split(" ")
    #ligne2 = [ele for ele in ligne if ele != '']
    count = 0
    x1 = float(ligne[0])
    x1_str = "{:.3f}".format(x1)
    x2 = float(ligne[1])
    x2_str = "{:.3f}".format(x2)
    x3 = float(ligne[2].strip())
    x3_str = "{:.3f}".format(x3)
    
    offset = 0
    if i == 8:
        i
    for j in range(len(ligne_b)):
        if ligne_b[j + offset] != '':
            count+=1
            if len(ligne_b[j + offset]) < 6 and count >= 7 and count <= 9:
                del(ligne_b[j + offset -1])
                offset-=1
            if len(ligne_b[j + offset]) > 6 and count >= 7 and count <= 9:
                ligne_b.insert(j+offset, '')
                offset+=1
        if count == 7 and x1_str != '':
            if len(x1_str) < 6:               
                ligne_b[j + offset] = x1_str
                ligne_b.insert(j+offset, '')
                offset +=1
            elif len(x1_str) > 6:
                ligne_b[j + offset] = x1_str
                del(ligne_b[j + offset -1])
                offset-=1
            else:
                ligne_b[j + offset] = x1_str
            x1_str = ''
        elif count == 8 and x2_str != '':
            if len(x2_str) < 6:               
                ligne_b[j + offset] = x2_str
                ligne_b.insert(j+offset, '')
                offset +=1
            elif len(x2_str) > 6:
                ligne_b[j + offset] = x2_str
                del(ligne_b[j + offset -1])
                offset-=1
            else:
                ligne_b[j + offset] = x2_str
            x2_str = ''
        elif count == 9 and x3_str != '':
            if len(x3_str) < 6:               
                ligne_b[j + offset] = x3_str
                ligne_b.insert(j+offset, '')
                offset +=1
            elif len(x3_str) > 6:
                ligne_b[j + offset] = x3_str
                del(ligne_b[j + offset -1])
                offset-=1
            else:
                ligne_b[j + offset] = x3_str
            x3_str = ''
        elif count > 9:
            break
    separator = " "
    ligne_a = separator.join(ligne_b)
    fileAfter.write(ligne_a)

S = len(lignes_B)
ind = index +N
while ind <S:
    fileAfter.write(lignes_B[ind])
    ind +=1
fileCoo.close()
fileBefore.close()
fileAfter.close()

    
    
    
    
    
    