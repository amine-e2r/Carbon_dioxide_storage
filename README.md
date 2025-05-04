# Modélisation du Stockage du Dioxyde de Carbone
<br><br/>
_Le cycle du carbone est un processus naturel qui régule la répartition du carbone entre l’atmosphère, la biosphère (forêts, sols), les océans et la lithosphère. Les forêts jouent un rôle essentiel dans ce cycle en absorbant le dioxyde de carbone (CO₂) atmosphérique par photosynthèse, le transformant en matière organique (glucides, cellulose) et en le stockant dans les arbres et les sols. La réaction chimique principale, la photosynthèse, est :_

$$
6 CO_2 + 6 H_2O \rightarrow C_6H_{12}O_6 + 6O_2
$$

<br><br/>
<br><br/>
***
### **Introduction**
***
Ce projet a pour objectif de modéliser les échanges de carbone entre les différents compartiments de l’écosystème à l’aide d’un système d’équations différentielles ordinaires (EDO). Bien sûr, nous ne prendrons pas en compte l’ensemble des compartiments écologiques, car cela mènerait à un système plus complexe.
<br><br/>
Pour aller plus loin dans la modélisation du cycle du carbone et obtenir une description plus réaliste des échanges, il est nécessaire d’utiliser des équations aux dérivées partielles (EDP). Contrairement au modèle présenté ici, basé sur des EDO qui ne tiennent pas compte de la dimension spatiale, les EDP permettent de modéliser la répartition du carbone dans l’espace — par exemple, sa diffusion dans le sol ou son transport dans l’atmosphère.
<br><br/>
Ces équations sont essentielles pour représenter la variabilité des sols, du climat ou de la végétation à grande échelle, ainsi que pour étudier des phénomènes tels que la décomposition de la matière organique en profondeur. Ces modèles, plus complexes, sont utilisés dans les simulateurs climatiques et les modèles globaux de végétation, et nécessitent des méthodes numériques avancées pour être résolus.
<br><br/>
<br><br/>
***
### **Formulation Mathématique**
***
Dans ce projet, nous étudions l'évolution de trois variables en fonction du temps :

- $C_A(t)$ : Quantité de carbone dans l'atmosphère  
- $C_T(t)$ : Quantité de carbone dans les arbres  
- $C_S(t)$ : Quantité de carbone dans les sols  

Le système d’équations gouvernant les échanges de carbone dans ce projet peut être modélisé par:

$$
\begin{align}
(1) \quad \frac{dC_A}{dt} &= -S(C_T) + \beta C_T + \delta C_S \\
(2) \quad \frac{dC_T}{dt} &= S(C_T) - \beta C_T - \delta C_T - \gamma C_T\\
(3) \quad \frac{dC_S}{dt} &= \gamma C_T + \delta C_T - \delta C_S
\end{align}
$$

avec:
  - $S(C_T) = \alpha C_T \left(1 - \dfrac{C_T}{K} \right)$: taux de séquestration du carbone dans les arbres,
  - $\beta C_T$: respiration des arbres vers l'atmosphère,
  - $\delta C_T$: transfert de carbone des arbres vers les sols,
  - $\delta C_S$: respiration des sols vers l'atmosphère,
  - $\gamma C_T$: litière végétale (feuilles mortes, débris organiques) transférée des arbres vers les sols.
<br><br/>
<br><br/>
***
### **Méthodes et Implémentation**
***
On pose:

$$
C(t) = \begin{pmatrix}
    C_A(t) \\
    C_T(t) \\
    C_S(t)
\end{pmatrix}
$$

et

$$
f(t,C(t)) = \begin{pmatrix}
    -\alpha C_T(t) \left(1 - \frac{C_T(t)}{K}\right) + \beta C_T(t) + \delta C_S(t) \\
    \alpha C_T(t) \left(1 - \frac{C_T(t)}{K}\right) - \beta C_T(t) - \delta C_T(t) - \gamma C_T(t) \\
    \gamma C_T(t) - \delta C_S(t) + \delta C_T(t)
\end{pmatrix}
$$

Ainsi le système revient à:

$$
\frac{dC(t)}{dt} = f(t,C(t))
$$

Si on discrétise uniformément l'intervalle de temps en N intervalles $[t_n, t_{n+1}]$, avec $n = \{1, 2, ..., N-1\}$, tel que $t_{n+1} = t_n + h$ et $h = \frac{1}{N}$.<br><br/>
On peut appliquer une methode d'Euler:

$$
\begin{align}
  C_n &= C(t_n)\\
  C_{A,n} &= C_A(t_n)\\
  C_{T,n} &= C_T(t_n)\\
  C_{S,n} &= C_S(t_n)
\end{align}
$$

$$
\begin{equation*}
    \int_{t_n}^{t_{n+1}} \frac{dC(t)}{dt} \ dt
    \ = \int_{t_n}^{t_{n+1}} f(t,C(t)) \ dt
\end{equation*}
$$

#### _Euler Implicite_

$$
\begin{equation*}
    C_{n+1} - C_n \approx (t_{n+1} - t_n)f(t_{n+1},C_{n+1})
\end{equation*}
$$

$$
\begin{equation*}
    C_{n+1} = C_n + h.f(t_{n+1},C_{n+1})
\end{equation*}
$$

Cette équation étants non-linéaire, il est mieux de se ramener à un problème de point fixe :

$$
F(C_{n+1}) = C_n + h \cdot f(t_n, C_{n+1})
$$

On considère alors $\forall n \ge 1$ la suite $(C_{n,k})_{k \in \mathbb{N}}$ telle que:

$$
\begin{equation}
\begin{cases}
  C_{n,k+1} &= \quad F(C_{n,k}) \\
  C_{n,0} &= \quad C_{n-1}
\end{cases}
\end{equation}
$$

Algorithme:

```
Fonction eulerImplicite(C0, F)
    Initialiser t ← t0
    Créer un tableau C de taille 3×1 rempli de zéros
    Créer une liste T contenant t0
    Mettre la première colonne de C égale à C0
    Initialiser k ← 1

    Tant que t < Tf faire
        t ← t + h
        Cn_plus_1 ← point_fixe(C[:, k-1], F)
        Ajouter Cn_plus_1 comme nouvelle colonne à C
        Ajouter t à la liste T
        k ← k + 1
    Fin Tant que

    Retourner C et T
Fin Fonction
```
```
Fonction point_fixe(X0, F, eps, max_iter)
    Convertir X0 en vecteur (si nécessaire)
    Initialiser Xk ← X0
    Initialiser Xk_1 ← X0

    Pour i allant de 0 à max_iter - 1 faire
        Xk_1 ← F(Xk, X0)
        Si la norme de (Xk - Xk_1) est inférieure à eps alors
            Sortir de la boucle
        Fin Si
        Xk ← Xk_1
    Fin Pour

    Retourner Xk_1
Fin Fonction
```

<br><br/>
#### _Euler Explicite_

$$
\begin{equation*}
    C_{n+1} - C_n \approx (t_{n+1} - t_n)f(t_{n},C_{n})
\end{equation*}
$$

$$
\begin{equation*}
    C_{n+1} = C_n + h.f(t_{n},C_{n})
\end{equation*}
$$

Algorithme:

```
Fonction eulerExplicite(C0)
    Initialiser t ← t0
    Créer un tableau C de taille 3×1 rempli de zéros
    Créer une liste T contenant t0
    Mettre la première colonne de C égale à C0
    Initialiser k ← 1

    Tant que t < Tf faire
        t ← t + h
        Cn_plus_1 ← C[:, k-1] + h × f(C[:, k-1])
        Ajouter Cn_plus_1 comme nouvelle colonne à C
        Ajouter t à la liste T
        k ← k + 1
    Fin Tant que

    Retourner C et T
Fin Fonction
```
<br><br/>
#### _Méthode de Trapèze_

$$
\begin{equation*}
        C_{n+1} = C_n + \frac{h}{2}[f(t_{n+1},C_{n+1})+f(t_n, C_n)]
    \end{equation*}
$$

Tout comme la méthode d'Euler Implicite, nous avons une équation non-linéaire, cette fois si nous allons utiliser une méthode de Newton qui va nous garantir une convergence quadratique et qui est donc plus performant que la méthode de point fixe, on pose alors:

$$
\begin{equation*}
        F(C_{n+1}) = C_n + \frac{h}{2}[f(t_{n+1},C_{n+1})+f(t_n, C_n)]
    \end{equation*}
$$

et

$$
f_{newton}(C_{n+1}) = F(C_{n+1}) - C_{n+1}
$$

Cela revient à résoudre $f_{newton}(C_{n+1}) = 0$

$$
\begin{align}
    \frac{df_{newton}(C_{n+1})}{d(C_{n+1})}
     &= \frac{d}{dC_{n+1}}F(C_{n+1}) - \frac{d}{dC_{n+1}}C_{n+1} \\
     &= \frac{d}{dC_{n+1}}(C_n +\frac{h}{2}(f(t_{n+1},C_{n+1})+f(t_n, C_n)) - I_3 \\
     &=\frac{h}{2}\frac{d}{dC_{n+1}}f(t_{n+1},C_{n+1}) - I_3 \\
\end{align}
$$

ainsi:

$$
\begin{equation*}
    \frac{df_{newton}(C_{n+1})}{d(C_{n+1})} = \frac{h}{2} 
    \begin{pmatrix}
        0 & -\alpha + \frac{2C_{T,n+1}}{K}+\beta & \delta \\
        0 & \alpha - \frac{2C_{T,n+1}}{K} - \beta - \delta - \gamma  & 0 \\ 
        0 & \gamma + \delta & - \delta \\
    \end{pmatrix} - I_3
\end{equation*}
$$

On à alors que chaque $C_{n}$ est solution du problème $f_{newton}(X) = 0$ telle que:

$$
\begin{equation}
\begin{cases}
    X_{n,0} &= C_{n-1} \\
    X_{n,k+1} &= X_{n,k} - (\frac{df_{newton}(X_{n,k-1})}{X_{n,k-1}})^{-1} \cdot f_{newton(X_{n,k-1})}
\end{cases}
\end{equation}
$$


Algorithme:

```
Fonction trapeze_newton(C0)
    Initialiser t ← t0
    Créer un tableau C de taille 3×1 rempli de zéros
    Créer une liste T contenant t0
    Mettre la première colonne de C égale à C0
    Initialiser k ← 1

    Tant que t < Tf faire
        t ← t + h
        Cn_plus_1 ← newton(C[:, k-1], f_newton, df_newton)
        Ajouter Cn_plus_1 comme nouvelle colonne à C
        Ajouter t à la liste T
        k ← k + 1
    Fin Tant que

    Retourner C et T
Fin Fonction
```

```
Fonction newton(X0, f, df, eps, max_iter)
    Convertir X0 en vecteur (si nécessaire)
    Initialiser Xk ← X0

    Pour i allant de 0 à max_iter - 1 faire
        Si la norme de f(Xk, X0) est inférieure à eps alors
            Sortir de la boucle
        Fin Si

        Xk ← Xk - inverse(df(Xk)) × f(Xk, X0)
    Fin Pour

    Retourner Xk
Fin Fonction
```

<br><br/>
<br><br/>
***
### **Resultas**
***
<br><br/>
<br><br/>

***
### **Impact des Paramètres sur la Séquestration du Carbone**
***

#### _Le paramètre 𝛼_
Le paramètre 𝛼 influence le taux de séquestration du carbone par les arbres. Il apparaît dans la fonction non linéaire :

$$
S(C_T) = \alpha C_T \left(1 - \frac{C_T}{K} \right)
$$

Cette fonction décrit la capture du CO₂ par les arbres. Lorsque 𝛼 augmente (pour un $C_T$ et un $K$ fixés), la valeur de $S(C_T)$ augmente également. Cela signifie que 𝛼 influence directement la rapidité d’absorption du carbone, ce qui conduit à un stockage plus rapide dans la biomasse végétale.

<br><br/>

#### _Le paramètre 𝐾_
Le paramètre $K$ représente la quantité maximale de carbone que les arbres peuvent stocker. En reprenant l'équation :

$$
S(C_T) = \alpha C_T \left(1 - \frac{C_T}{K} \right)
$$

on observe que faire varier $K$ modifie le taux de séquestration. Une valeur élevée de $K$ permet aux arbres d’absorber du carbone pendant une période plus longue avant d’atteindre leur capacité maximale.

<br><br/>

#### _Le paramètre β_
Le paramètre $\beta$, qui intervient dans le terme $\beta C_T$, représente la fraction de carbone que les arbres restituent à l’atmosphère par respiration. Plus $\beta$ est élevé, plus le retour de CO₂ dans l’air est important, ce qui réduit l’efficacité globale de la séquestration du carbone.

<br><br/>

#### _Le paramètre γ_
Le paramètre $\gamma$ contrôle le flux de matière organique morte (litière) des arbres vers les sols, modélisé par le terme $\gamma C_T$. Il ne prend pas en compte la respiration des arbres vers le sol, qui est décrite par le terme $\delta C_T$. Une valeur élevée de $\gamma$ favorise l’enrichissement des sols en matière organique, ce qui améliore leur fertilité et leur capacité à stocker du carbone sur le long terme.

<br><br/>

#### _Le paramètre δ_
Le paramètre $\delta$ intervient à deux niveaux : dans la respiration des arbres vers les sols ($\delta C_T$), et dans la respiration des sols vers l’atmosphère ($\delta C_S$). Il régule donc les pertes de carbone par respiration. Une valeur élevée de $\delta$ accélère le cycle du carbone, augmentant la quantité de CO₂ retournant à l’atmosphère. À l’inverse, un $\delta$ faible réduit ces pertes, ce qui favorise le stockage du carbone, notamment dans les sols. Cependant, une respiration trop faible peut ralentir le recyclage des nutriments essentiels à la croissance des plantes.
<br><br/>
<br><br/>

***
### **Analyse de l'impact des paramètres sur la modélisation**
***

#### Impact de $K$

La capacité de stockage maximale des arbres est représentée par le paramètre $K$. Si la quantité de carbone initiale dans les arbres $C_{T0}$ est supérieure à $K$, alors la fonction de séquestration du carbone, $S(C_T)$ est négative. Donc du carbone présent dans les arbres est rejeté dans le sol et dans l'atmosphère pour que $C_T$ passe sous $K$.

![Courbe avec \(C_T0 > K\) et \(h = 10\)](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/CT0%20sup%20K%20h%20%3D%2010.png)

![Courbe avec \(C_T0 > K\) et \(h = 0.1\)](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/CT0%20sup%20K%20h%20%3D%200.1.png)

On remarque qu’en effet au début le carbone dans l’air augmente et celui dans le sol a aussi une augmentation plus forte.  
Il y a une autre remarque : c'est le fait que la modélisation pour $h = 10$ et Euler Implicite ne marche pas. On peut imaginer que si $C_{T0} > K$ alors $C_T$ évolue rapidement pour passer sous $K$, tandis que $C_S$ et $C_A$ sont plus lents. Il faut donc un pas faible pour pouvoir capter correctement ces variations.


#### Absorption de carbone par les arbres uniquement

Dans ce cas, $\gamma = \beta = \delta = 0$ et $\alpha = 0.3$, les arbres absorbent du carbone sans jamais le relâcher, ni dans l’atmosphère, ni dans le sol. Dans ce scénario, le carbone absorbé par les arbres augmente jusqu'à atteindre la capacité maximale $K$, tandis que le carbone dans l'atmosphère diminue proportionnellement, jusqu'à ce que le niveau de $C_A$ atteigne $C_{A0} - (K - C_{T0})$.

![Courbe pour $\alpha = 0.3$, $\gamma = \beta = \delta = 0$](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/tout%20%3D%200%2C%20alpha%20!%3D.png)  
*Paramètre :* $K = 200$, $C_{A0} = 80$, $C_{T0} = 1$, $C_{S0} = 10$


#### Transfert de carbone vers le sol

Si on permet le transfert de carbone vers le sol grâce à la litière (feuilles mortes) en posant $\gamma = 0.2$  
Le carbone absorbé par les arbres est transféré dans le sol, ce qui permet une séquestration plus efficace.

![Courbe avec $\gamma = 0.2$](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/uniquement%20gamma%20alpha.png)  
Cette représentation nous montre qu’ici le transfert de carbone entre les arbres et le sol est assez important pour avoir une séquestration très forte et très rapide.

Si on pose $\gamma = 0.001$,  
![Courbe avec $\gamma = 0.001$](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/uniquement%20gamma%20alpha%20%2B%20faible.png)

On remarque que la séquestration est plus faible. **Mais quelque soit la valeur des paramètres $\alpha$ et $\gamma$ si $\beta = \delta = 0$ on n'atteindra pas d'équilibre** puisque aucun carbone n'est rejeté dans l'atmosphère.


#### Transfert arbres, sol et atmosphère

Lorsque le carbone est autorisé à être rejeté dans le sol par les arbres et dans l'atmosphère par le sol grâce à la respiration, avec des valeurs spécifiques pour $\delta = 0.01$, $\alpha = 0.1$, et $\gamma = 0.02$ et toujours $\beta = 0$, un équilibre est atteint. Dans cet équilibre, la quantité de carbone absorbée par les arbres est égale à celle rejetée dans l'atmosphère et le sol.

![Courbe avec $\delta = 0.01$, $\alpha = 0.1$, $\gamma = 0.02$](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/equilibre.png)  
*Paramètres :* $C_{A0} = 80$, $C_{T0} = 30$, $C_{S0} = 10$


#### Transfert entre arbres et atmosphère

On pose $\beta = 0.01$, on autorise la respiration directe entre les arbres et l'atmosphère.

![Courbe avec $\beta = 0.01$](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/equilibre%202.png)

On atteint l'équilibre à une quantité de séquestration plus faible car le système perd plus de carbone dans l'atmosphère en raison de l'augmentation de $\beta$, et la capacité d'absorption des arbres est limitée. On pourrait compenser cette perte en augmentant $\alpha$ et on pourra retrouver le même graphique que précédemment.

Si on augmente encore $\beta$ avec $\beta = 0.08$ on observe une augmentation du carbone dans l'atmosphère. La quantité de carbone rejetée est plus élevée que celle absorbée.

![Courbe avec $\beta = 0.08$](https://github.com/amine-e2r/Carbon_dioxide_storage/blob/main/courbe/augmentation%20ca.png)
### **Amélioration du Modèle**
***
Dans notre modèle, on autorise $C_A$ à devenir négatif, ce qui n'est pas réaliste d'un point de vue physique. Pour améliorer le modèle, on pourrait imposer une borne inférieure à $C_A$.
On peut faire cela en remplaçant une partie du code dans la boucle par:
```
Cn = #calcul du nouveau vecteur
Cn[0] = max(0, Cn[0]) #Cn[0] correspond au terme CA
C = np.append(C, np.transpose([Cn]))
```
Les océans sont les plus grand puit de carbone de la planète, il est donc légitime de les prendre en compge dans notre modélisation. 
En présence d'un océean, on peut rajouter une autre variable $C_O(t)$ la quantité de carbone dans l'eau. On a une interaction entre l'eau et l'atmosphère.

$\frac{dC_A(t)}{dt} = -S(C_T(t)) + \beta C_T(t) + \delta C_S(t) - \epsilon C_A(t) + \omega C_O(t)$
$\frac{dC_T(t)}{dt} = S(C_T(t)) - \beta C_T(t) - \delta C_T(t) - \gamma C_T(t)$
$\frac{dC_S(t)}{dt} = \gamma C_T(t) - \delta C_S(t) + \delta C_T(t)$
$\frac{dC_O(t)}{dt} = \epsilon C_A(t) - \omega C_O(t)$

avec $\epsilon C_A(t)$ qui représente le taux de séquestration du carbone dans l'océean et $\omega C_O(t)$ décrit l’effet de respiration de l'océan vers l'atmosphère.

***
### **Conclusion**
***
<br><br/>
<br><br/>
***
### **Annexe**
***




















    
