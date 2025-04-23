import numpy as np
import matplotlib.pyplot as plt

#Variables globales

#Constantes ******************************
alpha = 0.1
beta = 0.02
gamma = 0.03
delta = 0.01
K = 100
#Valeurs trouvées grâce à chat gpt et internet
#Constantes ******************************

#Autres variables ************************
max_iter = 10000
eps = 1e-6
Tf = 1000 #temps final en secondes


    #Conditions initiales
CA0 = 80
CT0 = 20
CS0 = 10
#Valeurs trouvées grâce à chat gpt et internet
C0 = [CA0, CT0, CS0]
#Plus tard, le tableau de tableau C contiendra à chaque ligne n le tableau Cn (n ième itération de la methode de newton)
t0 = 0 #temps initial en secondes
h = 0.1 # il y a 0.1 secondes entre tn et tn+1 (on pourra le réduire ou l'augmenter en fonction de nos besoins)



#Autres variables ************************

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
    X0 = np.array(X0)
    Xk = X0 #
    Xk_1 = X0 #Cn,k+1
    for i in range(_max_iter):
        Xk_1 = _F(Xk, X0)
        if(np.linalg.norm(Xk - Xk_1) < eps):
            break
        Xk = Xk_1
    return Xk_1

def eulerImplicite(_C0,_F):
    #On calcul tout les Cn grâce à la méthode du point fixe
    #C0 appartient à R3
    t = t0
    C = np.zeros((3,1))
    T = [t0]
    C[:,0] = _C0
    k = 1
    while(t < Tf):
        t = t+h
        C = np.append(C, np.transpose([point_fixe(C[:,k-1],_F)]), axis=1) #ajoute comme on veut le nouveau Cn+1 au tableau C
        T.append(t)
        k = k+1
    return C, T


def eulerExplicite(_C0):
    #On calcul directement les Cn
    #C0 appartient à R3
    t = t0
    C = np.zeros((3,1))
    T = [t0]
    C[:,0] = _C0
    k = 1
    while(t < Tf):
        t = t+h
        C = np.append(C, np.transpose( [C[:,k-1] + h*f(C[:,k-1])] ), axis = 1) #ajoute comme on veut le nouveau Cn+1 au tableau C
        T.append(t)
        k = k+1
    return C, T
#methodes numeriques**********************


C1, T1 = eulerImplicite(C0, F1)
C2, T2 = eulerExplicite(C0)

plt.subplot(2,1,1)
plt.plot(T1, C1[0], label='CA(t)')
plt.plot(T1, C1[1], label='CT(t)')
plt.plot(T1, C1[2], label='Cs(t)')
plt.legend()
plt.title("Point fixe à partir d'Euler Implicite")

plt.subplot(2,1,2)
plt.plot(T2, C2[0], label='CA(t)')
plt.plot(T2, C2[1], label='CT(t)')
plt.plot(T2, C2[2], label='CS(t)')
plt.legend()
plt.xlabel("t")
plt.title("Euler explicite")

plt.show()
