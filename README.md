# Modeling of carbon dioxide storage
<br><br/>
<br><br/>
<br><br/>
***
### **Introduction**
***
_Ce projet a pour objectif de modéliser les échanges de carbone entre les différents compartiments de l’écosystème à l’aide d’un système d’équations différentielles ordinaires (EDO). Bien sûr, nous ne prendrons pas en compte l’ensemble des compartiments écologiques, car cela mènerait à un système plus complexe.
<br><br/>
Pour aller plus loin dans la modélisation du cycle du carbone et obtenir une description plus réaliste des échanges, il est nécessaire d’utiliser des équations aux dérivées partielles (EDP). Contrairement au modèle présenté ici, basé sur des EDO qui ne tiennent pas compte de la dimension spatiale, les EDP permettent de modéliser la répartition du carbone dans l’espace — par exemple, sa diffusion dans le sol ou son transport dans l’atmosphère.
<br><br/>
Ces équations sont essentielles pour représenter la variabilité des sols, du climat ou de la végétation à grande échelle, ainsi que pour étudier des phénomènes tels que la décomposition de la matière organique en profondeur. Ces modèles, plus complexes, sont utilisés dans les simulateurs climatiques et les modèles globaux de végétation, et nécessitent des méthodes numériques avancées pour être résolus._
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

