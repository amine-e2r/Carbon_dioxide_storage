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
Dans ce projet Nous étudions l'évolution de trois variables en fonction du temps :

- **CA(t)** : Quantité de carbone dans l'atmosphère  
- **CT(t)** : Quantité de carbone dans les arbres  
- **CS(t)** : Quantité de carbone dans les sols  

Le système d’équations gouvernant les échanges de carbone dans ce projet peut être modélisé par:

&`
dCA/dt = -S(CT) + β·CT + δ·CS        (1)
dCT/dt = S(CT) - β·CT - δ·CT - γ·CT  (2)
dCS/dt = γ·CT + δ·CT - δ·CS          (3)
`&






A rajouter :
Méthode point fixe mais avec trapèze au lieux de rectangle, essayer avec simpson aussi + autres idées,
Affiner le modèle pour coller plus aux vrais observations (lucah c le pro il es trop fort en CO2)
Voir l'impact des différentes constantes, et voir différentes conditions initiales
