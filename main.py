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


"""
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
def point_fixe(C0, F, eps=eps, max_iter = max_iter):
    C = np.array(C0)
    k = 1
    while True:
        C = np.append(C,F(k,C))
        if(k >= max_iter or np.linalg.norm(C[k] - C[k-1],2) < eps):
            break
        k+=1
    return C[k]
#methodes numeriques**********************
"""
#fonctions********************************

def f(Cn): 
    #Ici Cn est un tableau de taille 3(numpy array ou bien juste un array)
    #n correspond à l'indice de tn
    #renvoie un np array
    y1 = Cn[1] * (-alpha * (1 - (Cn[1])/K) + beta) + delta * Cn[2]
    y2 = Cn[1] * (alpha * (1 - (Cn[1])/K) - beta - delta - gamma)
    y3 = Cn[1] * (gamma + delta) - delta * Cn[2]
    return np.array([y1,y2,y3])

def F(Cnk_suiv, Cnk_prec):
    #l'argument C peut etre un np array ou un array (les deux marchent)
    Cnk_prec = np.array(Cnk_prec)
    Cnk_suiv = np.array(Cnk_suiv)
    return Cnk_prec + h*f(Cnk_suiv)

#fonctions********************************



#methodes numeriques**********************

def point_fixe(X0, _F, _eps=eps, _max_iter = max_iter):
    #résout X  = _F(X) par methode du point fixe
    C = np.array(C)
    Xk = X0 #
    Xk_1 = X0 #Cn,k+1
    for i in range(_max_iter):
        Xk_1 = _F(Xk, X0)
        if(np.norm(Xk, Xk_1) < eps):
            break
        Xk = Xk_1
    return Xk_1

def eulerImplicite(C0,_F):
    N = Tf//h
    C = np.zeros((3,N))
    C[:,0] = C0
    for k in range(1,N):
        C[:,k] = point_fixe(C[:,k-1],_F)
    return C
        
#methodes numeriques**********************
