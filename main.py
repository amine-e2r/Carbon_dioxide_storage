import numpy as np


#Variables globales

#Constantes ******************************
alpha = 1
beta = 1
gamma = 1
delta = 1
K = 1
#Constantes ******************************

#Autres variables ************************
max_iter = 10000
eps = 1e-6
Tf = 1000 #temps final en secondes


    #Conditions initiales
CA0 = 0.1
CT0 = 0.1
CS0 = 0.1
C0 = [CA0, CT0, CS0]
#Plus tard, le tableau de tableau C contiendra à chaque ligne n le tableau Cn (n ième itération de la methode de newton)
t0 = 0 #temps initial en secondes
h = 0.1 # il y a 0.1 secondes entre tn et tn+1 (on pourra le réduire ou l'augmenter en fonction de nos besoins)



#Autres variables ************************



#fonctions********************************

def f(n, C): 
    #Ici C est un tableau de tableaux de taille 3 chacun(numpy array ou bien juste un array)
    #n correspond à l'indice de tn
    #renvoie un np array
    y1 = C[n][1] * (-alpha * (1 - (C[n][1])/K) + beta) + delta * C[n][2]
    y2 = C[n][1] * (alpha * (1 - (C[n][1])/K) -beta - delta - gamma)
    y3 = C[n][1] * (gamma + delta) - delta * C[n][2]
    return np.array([y1,y2,y3])

def F(n, C):
    #l'argument C peut etre un np array ou un array (les deux marchent)
    C = np.array(C)
    return C[n-1] + f(n,C)

#fonctions********************************



#methodes numeriques**********************
"""
def point_fixe(n, C, _F, _eps=eps, _max_iter = max_iter):
    C = np.array(C)
    Ck = C[0]
    for i in range(_max_iter):
        if(C[])
"""
#methodes numeriques**********************
