### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ b1f454c0-1be0-11f1-1707-8dc83ac1495c
begin
	using PlutoUI
	using Plots
end

# ╔═╡ 22a2d086-00f9-43fb-ab60-c6c7df65f0c9
html"""
<div style="display: flex; justify-content: center; background-color: white; padding: 10px; border-radius: 8px;">
	<img src="https://raw.githubusercontent.com/GaldinoMagalhaes/geomechanics_classes/98574ff7d01e1b5c722538b70605a89c3fceddfb/images/logo_GRMM.png" width="800">
</div>
"""

# ╔═╡ 27c9ebff-4e97-4a7d-aba0-1c26583d2de8
html"""
<p style="background-color:lightgray" xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">&nbsp&nbsp💻&nbsp<b>Anisotropy (Part I)</b></span> par <span property="cc:attributionName">Magalhães et Cacciari (2026)</span> est licencié sous <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:10px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>
"""

# ╔═╡ 4942766d-0e7c-44fd-9dcd-379e8be02cdb
begin
	md"""
	# Considérations initiales
	
	Le matériel suivant a été développé comme support pour les cours de premier cycle du cours de Mécanique des Roches à Polytechnique Montréal. Ce notebook fait partie d'une série en trois parties sur le comportement anisotrope des roches, composée de :
	
	* Partie I - Déformation des roches anisotropes
	* Partie II - Critères de rupture anisotropes
	* Partie III - Tenseurs de structure

	Le contenu de ce notebook — _Partie I_ — commence par la relation constitutive entre contrainte et déformation en utilisant une forme généralisée de la loi de Hooke. Initialement, la condition d'anisotropie maximale est présentée, évoluant vers des cas spécifiques de symétrie jusqu'à atteindre la condition théorique d'isotropie complète. Bien que le matériel s'appuie sur des concepts introductifs de physique, il est conçu pour fournir une contextualisation géologique et géotechnique des connaissances. 
	
	Les matériaux sont élaborés par le groupe de recherche Geomechanics & Rock Mass Modeling et peuvent être utilisés à des fins éducatives à condition que la référence appropriée soit donnée. Pour des suggestions, corrections ou pour tout contact général, veuillez utiliser l'adresse courriel [gabrielgaldinodm@gmail.com](mailto:gabrielgaldinodm@gmail.com). 
	
	💡 **Le contenu est principalement destiné à la révision et à la pratique des concepts, et les étudiants sont encouragés à continuer à consulter la littérature.** 

	"""
end

# ╔═╡ 9c03c6b6-63d8-46a3-8c84-8a0e9f424055
md"""

# Loi de Hooke généralisée

Du point de vue de la **mécanique des milieux continus**, la relation entre contrainte et déformation pour un matériau linéaire et élastique peut être décrite par une forme tensorielle généralisée de la **loi de Hooke**, donnée par la relation constitutive suivante :

$\sigma_{ij} =C_{ijkl}\epsilon_{kl}\quad\quad (Eq.01)$

où :

* **$\epsilon$** : Tenseur de déformation du second ordre.  
* **$\sigma$** : Tenseur des contraintes du second ordre.  
* **$C$** : Tenseur de rigidité élastique (quatrième ordre) qui décrit la résistance du matériau à la déformation. Il est l’inverse du tenseur de compliance ($S_{ijkl}$).  

L’Eq. 01 utilise la **notation indicielle** (convention de sommation d’Einstein), où la répétition d’indices implique une sommation sur toutes les directions spatiales. Il s’agit de la forme la plus générale de la loi de Hooke, en supposant les principes physiques d’élasticité et de linéarité entre contrainte et déformation.

## Relations de symétrie

En partant d’une condition générale, nous pouvons réfléchir aux relations possibles de symétrie. En commençant par le cas de **symétrie maximale**, nous obtenons les **matériaux isotropes**, dans lesquels aucune direction préférentielle n’est attribuée au matériau.

L’équilibre des moments dans un élément infinitésimal exige que le tenseur des contraintes soit symétrique, c’est-à-dire $\sigma_{ij} = \sigma_{ji}$, et, par extension, la même symétrie s’applique au tenseur des déformations, $\epsilon_{ij} = \epsilon_{ji}$. Ainsi, les deux tenseurs du second ordre possèdent seulement **six composantes indépendantes**.

En revenant à l’équation constitutive (Eq.01), si $\epsilon_{kl}$ est symétrique, seules les combinaisons avec $k \leq l$ sont physiquement indépendantes. De plus, le tenseur des contraintes résultant ($\sigma_{ij}$) doit également être symétrique. Cela implique les **symétries mineures** du tenseur élastique :

* Si la déformation est symétrique, changer l’ordre de $kl$ ne modifie pas l’égalité  
$C_{ij\color{#F15A22}kl} = C_{ij\color{#F15A22}lk} \quad\quad (Eq.02)$

* Si le tenseur résultant est symétrique, l’ordre $ij$ ne modifie pas l’égalité  
$C_{\color{#F15A22}ji\color{black}kl} = C_{\color{#F15A22}ij\color{black}kl} \quad\quad (\text{Eq.03})$

Une **troisième symétrie** (symétrie majeure) peut également être déduite à partir de l’énergie de déformation élastique, ce qui garantit que :

$C_{ijkl} = C_{klij}$

Pour plus de détails, voir **Exercice 1**.

Une fois ces symétries établies, il est possible d’adopter la **notation de Voigt** (Eq.04), un mappage algébrique dans lequel les tenseurs symétriques du second ordre sont représentés par des **vecteurs de dimension six**, tandis que le tenseur du quatrième ordre devient une **matrice 6×6**.

$[\epsilon] = [C][\sigma]\quad\quad (Eq.04)$

Ou sous forme matricielle complète :

$$\begin{bmatrix} 
\sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} 
\end{bmatrix} 
=
\begin{bmatrix}
C_{11} & C_{12} & C_{13} & C_{14} & C_{15} & C_{16} \\
\color{#F15A22}{C_{21}} & C_{22} & C_{23} & C_{24} & C_{25} & C_{26} \\
\color{#F15A22}{C_{31}} & \color{#F15A22}{C_{32}} & C_{33} & C_{34} & C_{35} & C_{36} \\
\color{#F15A22}{C_{41}} & \color{#F15A22}{C_{42}} & \color{#F15A22}{C_{43}} & C_{44} & C_{45} & C_{46} \\
\color{#F15A22}{C_{51}} & \color{#F15A22}{C_{52}} & \color{#F15A22}{C_{53}} & \color{#F15A22}{C_{54}} & C_{55} & C_{56} \\
\color{#F15A22}{C_{61}} & \color{#F15A22}{C_{62}} & \color{#F15A22}{C_{63}} & \color{#F15A22}{C_{64}} & \color{#F15A22}{C_{65}} & C_{66}
\end{bmatrix}
\begin{bmatrix} 
\epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ 2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12} 
\end{bmatrix}$$

Dans cette notation, la matrice $[C]$ est appelée **matrice constitutive**. Les relations de symétrie permettent ici de réduire le nombre de constantes indépendantes de **36 à 21**, caractérisant l’état **anisotrope général**, communément appelé **triclinique**. Les éléments dépendants sont indiqués en $\color{#F15A22}orange$.

> **Remarque :** La multiplication par 2 des déformations de cisaillement garantit que le travail interne calculé à partir du vecteur de Voigt est identique à celui calculé à partir du tenseur original. Dans les étapes suivantes, la substitution $\gamma_{ij} = 2\epsilon_{ij}$ sera adoptée.

"""

# ╔═╡ d7887b31-b1ed-447a-9887-2237574c093f
md"""
## Exercice 01
À partir du concept d’énergie spécifique de déformation élastique, $$U = \frac{1}{2} \sigma_{ij} \epsilon_{ij}$$, expliquez pourquoi le tenseur élastique doit satisfaire la symétrie majeure ($$C_{ijkl} = C_{klij}$$).
"""

# ╔═╡ 95d91406-842e-442c-9a07-0109252942cc
md"""
## Exercice 02
Quelles sont les limitations pratiques (en laboratoire) pour la détermination de toutes les composantes de la matrice constitutive d’une roche ?
"""

# ╔═╡ ad7b4db9-1116-4e26-a861-4e8b42342c28
md"""
Dans les étapes suivantes, nous allons détailler les principaux modèles:
"""

# ╔═╡ 5718c387-4455-4a49-933d-63ffef23fffa
html"""
<div style="display: flex; justify-content: center; background-color: white; padding: 10px; border-radius: 8px;">
	<img src="https://raw.githubusercontent.com/GaldinoMagalhaes/geomechanics_classes/b7ccd1847c908bfa1ce519ff07467c0326341b05/images/models_fr.png" width="800">
</div>
"""

# ╔═╡ 0bcd566e-b7dd-40f8-8bb6-ae67db90bc41
md"""
# Modèle Orthotrope
Compte tenu de la complexité des modèles complètement anisotropes, certains modèles avec symétrie élastique sont introduits pour expliquer des scénarios géologiques spécifiques. Le modèle orthotrope en particulier est caractérisé par trois plans de symétrie élastique. Ce type de symétrie peut être observé, par exemple, dans des roches avec trois familles de fractures dans des directions perpendiculaires. Au niveau minéralogique, le modèle orthotrope décrit des cristaux de symétrie orthorhombique.  

Pour le contexte de l’orthotropie, l’exigence d’invariance de la loi constitutive sous les opérations de symétrie impose des restrictions qui forcent certains coefficients constitutifs à être nuls. Cette condition modifie la matrice de rigidité (ou de compliance), la rendant dépendante d’un nombre plus faible de constantes indépendantes (moins de degrés de liberté). Cette relation sera démontrée dans les sections suivantes.

## Relations de symétrie
Dans le contexte du modèle orthotrope, nous avons trois plans de symétrie $x_2x_3$, $x_1x_3$ et $x_1x_2$, avec les relations suivantes :

* Réflexion dans le plan $x_2x_3$ (échange $x_1\rightarrow -x_1$)
* Réflexion dans le plan $x_1x_3$ (échange $x_2\rightarrow -x_2$)
* Réflexion dans le plan $x_1x_2$ (échange $x_3\rightarrow -x_3$)

En développant le cas de réflexion dans $x_2x_3$, nous avons que la transformation devient $(x'_1, x'_2, x'_3) = (-x_1,x_2, x_3)$. Comme la description des composantes de contrainte et de déformation change de signe, dans ce cas, nous pouvons généraliser que les composantes avec un nombre impair d’indices 1 changent de signe, tandis que celles avec un nombre pair (ou aucun) restent invariantes.

| Composante | Contrainte ($\sigma$) | Déformation ($\varepsilon$) | Statut |
| :---: | :---: | :---: | :--- |
| $11 \to 1$ | $\sigma_1' = \sigma_1$ | $\varepsilon_1' = \varepsilon_1$ | Invariante |
| $22 \to 2$ | $\sigma_2' = \sigma_2$ | $\varepsilon_2' = \varepsilon_2$ | Invariante |
| $33 \to 3$ | $\sigma_3' = \sigma_3$ | $\varepsilon_3' = \varepsilon_3$ | Invariante |
| $23 \to 4$ | $\sigma_4' = \sigma_4$ | $\varepsilon_4' = \varepsilon_4$ | Invariante |
| $13 \to 5$ | $\sigma_5' = -\sigma_5$ | $\varepsilon_5' = -\varepsilon_5$ | **Inverse** |
| $12 \to 6$ | $\sigma_6' = -\sigma_6$ | $\varepsilon_6' = -\varepsilon_6$ | **Inverse** |


Si l’invariance est valide, la relation constitutive doit se comporter de la même manière avant et après la réflexion, c’est-à-dire que les équations Eq.05 et Eq.06 doivent être valides respectivement avant et après la réflexion. 

$\epsilon = [C]\sigma\quad\quad (Eq.05)$ 
$\epsilon' = [C]\sigma'\quad\quad(Eq.06)$

En développant les équations précédentes pour la composante $\epsilon_1$ pour le scénario initial et après réflexion, nous avons :

$\varepsilon _1=C_{11}\sigma _1+C_{12}\sigma _2+C_{13}\sigma _3+C_{14}\sigma _4+C_{15}\sigma _5+C_{16}\sigma _6 \quad\quad (Eq.07)$

$\varepsilon _1'=C_{11}\sigma _1'+C_{12}\sigma _2'+C_{13}\sigma _3'+C_{14}\sigma _4'+C_{15}\sigma _5'+C_{16}\sigma _6' \quad\quad (Eq.08)$

En substituant les signes de l’Eq.07 :

$\varepsilon _1=C_{11}\sigma _1+C_{12}\sigma _2+C_{13}\sigma _3+C_{14}\sigma _4-C_{15}\sigma _5-C_{16}\sigma _6 \quad\quad (Eq.09)$

Pour que l’égalité soit maintenue après la rotation, il est nécessaire que les composantes $C_{15}$ et $C_{16}$ soient nulles. En développant les autres lignes de la multiplication ($\epsilon_2, \epsilon_3$), nous obtenons :

$C_{15}=C_{16}=C_{25}=C_{26}=C_{35}=C_{36}=0$

En réalisant la même procédure pour les autres plans de symétrie, nous obtenons la relation constitutive caractéristique du modèle orthotrope :

$$\left[ \begin{matrix}\varepsilon _1\\ \varepsilon _2\\ \varepsilon _3\\ \varepsilon _4\\ \varepsilon _5\\ \varepsilon _6\end{matrix}\right] =\left[ \begin{matrix}C_{11}&C_{12}&C_{13}&0&0&0\\ C_{12}&C_{22}&C_{23}&0&0&0\\ C_{13}&C_{23}&C_{33}&0&0&0\\ 0&0&0&C_{44}&0&0\\ 0&0&0&0&C_{55}&0\\ 0&0&0&0&0&C_{66}\end{matrix}\right] \left[ \begin{matrix}\sigma _1\\ \sigma _2\\ \sigma _3\\ \sigma _4\\ \sigma _5\\ \sigma _6\end{matrix}\right]\quad\quad (Eq.10)$$

## Constantes d’Ingénierie ($E$ et $\nu$)
L’Eq.10 décrit la relation constitutive du modèle orthotrope basée sur les entrées de la matrice de rigidité. Cependant, nous pouvons réécrire la loi de Hooke en fonction du module de Young ($E$) et de la constante de Poisson ($\nu$) pour une meilleure interprétabilité d’un point de vue ingénierie. 

*Attention !*
Dans ce cas, la manipulation algébrique de la matrice de compliance est plus simple que celle de la matrice de rigidité. Ainsi, réécrivons $[\epsilon] = [S][\sigma]$ comme :

$$\begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ 2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12} \end{bmatrix} = 
\begin{bmatrix} 
1/E_1 & -\nu_{21}/E_2 & -\nu_{31}/E_3 & 0 & 0 & 0 \\
-\nu_{12}/E_1 & 1/E_2 & -\nu_{32}/E_3 & 0 & 0 & 0 \\
-\nu_{13}/E_1 & -\nu_{23}/E_2 & 1/E_3 & 0 & 0 & 0 \\
0 & 0 & 0 & 1/G_{23} & 0 & 0 \\
0 & 0 & 0 & 0 & 1/G_{13} & 0 \\
0 & 0 & 0 & 0 & 0 & 1/G_{12} 
\end{bmatrix}
\begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{bmatrix}\quad\quad (Eq.11)$$

Dans l’Eq.11, nous avons **9 termes indépendants**, à savoir :

* 3 Modules de Young : $E_1, E_2, E_3$ (rigidité selon chaque axe).
* 3 Modules de Cisaillement : $G_{12}, G_{13}, G_{23}$ (rigidité au cisaillement dans les plans).
* 3 Coefficients de Poisson : $\nu_{12}, \nu_{13}, \nu_{23}$ (les 3 autres constantes, $\nu_{21}, \nu_{31}, \nu_{32}$, dépendent de celles-ci en raison de la symétrie de la matrice). 

D’une part, les paramètres d’élasticité $E$ et $\nu$ sont connus et traités ici seulement de manière directionnelle, le paramètre $G_{vh}$ est introduit pour capturer l’effet de glissement entre les couches (cisaillement entre les lamelles).

Bien que davantage de termes de Poisson apparaissent dans la matrice, seuls 3 sont indépendants, étant donné la validité de la relation de symétrie $\frac{\nu_{ij}}{E_i} = \frac{\nu_{ji}}{E_j}$. 
"""

# ╔═╡ 4bce3618-e8a3-4940-a21e-fcdbd61fc1d2
md"""
## Exercice 03
Ceci est un exercice conceptuel. Utilisez la description du problème pour les questions A à E.

**Description :** Considérez un massif rocheux avec 3 ensembles de plans de discontinuités orientés orthogonalement entre eux. Le premier ensemble est le feuilletage métamorphique ($S_n$), pénétratif et très fréquent, orienté selon la direction $x_1$ ; le deuxième ensemble de plans ($F1$) consiste en des fractures peu fréquentes dans la direction $x_2$ ; dans la direction $x_3$, on observe des fractures ($F2$) avec des caractéristiques très similaires à $F1$, mais considérablement plus fréquentes.

A) Quelle est la relation attendue entre le module de Young ($E_i$) et l’intensité de fracturation dans les directions $x_2$ et $x_3$ ? Pourquoi cela se produit-il ?
"""

# ╔═╡ 916a9876-b586-46a3-848f-af54e2e92a83
md"""
B) Est-il attendu que $\nu_{23}$ et $\nu_{32}$ soient égaux ? Justifiez votre réponse en vous basant sur la matrice de compliance.
"""

# ╔═╡ 535aaadf-0402-421f-9035-fdb502e20f7b
md"""
C) Sous un état de contrainte uniaxiale appliquée dans la direction $x_3$, attend-on une déformation axiale plus grande ou plus faible par rapport à la même contrainte appliquée dans la direction $x_2$ ? Prouvez mathématiquement et justifiez en vous basant sur la configuration structurelle du massif.
"""

# ╔═╡ 85f9d2b4-e450-4f8e-bec9-f1ed6f15a2e9
md"""
D) Une équipe responsable du dynamitage de roche doit réaliser des forages pour l’insertion d’explosifs. En considérant uniquement la réponse élastique anisotrope du massif, et en négligeant les effets dynamiques et non linéaires, quelle serait la direction la plus favorable pour un dynamitage plus efficace ? Justifiez.
"""

# ╔═╡ 70cc842d-6b86-47d2-a356-5746cd0f7524
md"""
E) Quelles conditions pourraient justifier l’adoption d’un système encore plus simplifié (moins d’axes d’anisotropie) ?
"""

# ╔═╡ 6e83a91d-bba7-4c84-9e18-429719b1ee0f
md"""
# Modèle Transversal Isotrope

Alors que les matériaux orthotropes sont invariants par rapport aux **réflexions** et possèdent 3 plans de symétrie, les matériaux **transversalement isotropes** se caractérisent par des propriétés constantes le long d’un plan de symétrie, tandis qu’elles varient dans la direction orthogonale – c’est-à-dire que toute **rotation** autour de l’axe de symétrie n’altère pas la réponse mécanique du matériau. Ce type de symétrie est une approximation suffisamment bonne pour différents contextes géologiques, comme les roches sédimentaires à stratification bien définie ou les roches métamorphiques avec un plan de foliation régulier.  

> **Note** : il est important de noter ici que, algébriquement, les concepts de réflexion et de rotation sont différents, ce qui souligne d’ailleurs la différence entre les modèles orthotrope et transversalement isotrope – dans ce sens, un matériau TI est, par définition, un cas particulier de l’orthotrope. 

## Relations de symétrie

Dans l’Eq.10, nous avons présenté la matrice de rigidité pour le modèle orthotrope en supposant l’axe $x_2x_3$ comme référence de symétrie. Considérons maintenant le même plan comme référence pour la symétrie du modèle transversalement isotrope. Pour cela, nous reprenons le modèle précédent, mais en imposant l’invariance à la rotation autour de l’axe $x_1$ pour tout $\theta$ :

$$\begin{cases}
x_1' = x_1 \\
x_2' = x_2 \cos\theta - x_3 \sin\theta \\
x_3' = x_2 \sin\theta + x_3 \cos\theta
\end{cases}$$

L’invariance pour tout $\theta$ impose :

**A)** $C_{22} = C_{33} \quad\text{et}\quad C_{12} = C_{13}$

La rigidité dans les axes $x_2$ et $x_3$ doit être identique pour garantir l’isotropie dans le plan $x_2x_3$. Même chose pour les déformations latérales en $x_1$. 

**B)** $C_{55}=C_{66}$

L’élément $C_{66}$ est lié au cisaillement dans le plan $1-2$, et $C_{55}$ au plan $1-3$. Comme l’axe $x_1$ est le centre de rotation, tout plan qui le contient doit avoir la même résistance au cisaillement.

**C)** Isotropie pour $C_{44}$

La relation entre $E$, $\nu$ et $G$ pour les matériaux isotropes est donnée par :

$G = \frac{E}{2(1 + \nu)}$

Comme un matériau TI se comporte de manière isotrope dans $x_2x_3$, cette relation doit également être valide pour ce plan. Dans la notation de la matrice de rigidité, nous avons :

$$C_{44} = \frac{C_{22} - C_{23}}{2}$$.

Ainsi, nous pouvons écrire la loi de Hooke pour les matériaux transversalement isotropes en fonction de la matrice de rigidité comme :

$$\begin{bmatrix} 
\sigma_1 \\ \sigma_2 \\ \sigma_3 \\ \sigma_4 \\ \sigma_5 \\ \sigma_6 
\end{bmatrix} =
\begin{bmatrix}
C_{11} & C_{12} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{22} & C_{23} & 0 & 0 & 0 \\
C_{12} & C_{23} & C_{22} & 0 & 0 & 0 \\
0 & 0 & 0 & \frac{C_{22}-C_{23}}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & C_{66} & 0 \\
0 & 0 & 0 & 0 & 0 & C_{66}
\end{bmatrix}
\begin{bmatrix} 
\varepsilon_1 \\ \varepsilon_2 \\ \varepsilon_3 \\ \gamma_{23} \\ \gamma_{13} \\ \gamma_{12} 
\end{bmatrix}\quad\quad (Eq.12)$$

## Constantes d’Ingénierie ($E$ et $\nu$)

L’Eq.01 de ce matériau décrit le scénario généralisé de la loi de Hooke où aucun plan de symétrie n’est supposé – ce qui implique 21 paramètres indépendants. Ensuite, nous avons discuté du modèle orthotrope, où les conditions de symétrie réduisent la formalisation à 9 paramètres. Dans le modèle transversalement isotrope, la formalisation d’ingénierie permet d’utiliser seulement les 5 constantes suivantes pour décrire le modèle :

* **$E_h$** : module de Young dans le plan horizontal (stratification) ;
* **$E_v$** : module de Young dans le plan vertical ;
* **$\nu_{hh}$** : effet de Poisson dans le plan horizontal ;
* **$\nu_{vh}$** : effet de Poisson décrivant la déformation horizontale causée par une contrainte verticale ;
* **$G_{vh}$** : module de cisaillement dans les plans verticaux ;

En plus de ces 5 paramètres, nous avons $G_{hh}$ qui est défini en fonction des propriétés élastiques longitudinales du plan – c’est-à-dire que ce n’est pas un paramètre indépendant et donc pas un degré de liberté du modèle :

$G_{hh} = \frac{E_h}{2(1 +\nu_{hh})}\quad\quad (Eq.13)$

Ainsi, l’équation constitutive en fonction des paramètres d’ingénierie devient :

$$\begin{bmatrix} 
\varepsilon_1 \\ \varepsilon_2 \\ \varepsilon_3 \\ \gamma_{23} \\ \gamma_{13} \\ \gamma_{12} 
\end{bmatrix} =
\begin{bmatrix}
\frac{1}{E_v} & -\frac{\nu_{vh}}{E_v} & -\frac{\nu_{vh}}{E_v} & 0 & 0 & 0 \\
-\frac{\nu_{hv}}{E_h} & \frac{1}{E_h} & -\frac{\nu_{hh}}{E_h} & 0 & 0 & 0 \\
-\frac{\nu_{hv}}{E_h} & -\frac{\nu_{hh}}{E_h} & \frac{1}{E_h} & 0 & 0 & 0 \\
0 & 0 & 0 & \frac{1}{G_{hh}} & 0 & 0 \\
0 & 0 & 0 & 0 & \frac{1}{G_{vh}} & 0 \\
0 & 0 & 0 & 0 & 0 & \frac{1}{G_{vh}}
\end{bmatrix}
\begin{bmatrix} 
\sigma_1 \\ \sigma_2 \\ \sigma_3 \\ \sigma_4 \\ \sigma_5 \\ \sigma_6 
\end{bmatrix}\quad\quad (Eq.14)$$
"""

# ╔═╡ 72ba316d-6105-443b-a2b9-6bb8d49c9ece
md"""
## Exercice 04

Considérons un échantillon de roche métamorphique fortement feuilletée, avec une foliation plane-parallèle et régulière. Le tableau ci-dessous décrit des valeurs de référence issues de tests en laboratoire pour les paramètres élastiques de ce type de roche. Considérons que la foliation est orientée perpendiculairement à l’axe $x_3$ (vertical).

| Paramètre | Valeur |
| :--- | :--- |
| Module de Young Horizontal ($E_h$) | $30\text{ GPa}$ |
| Module de Young Vertical ($E_v$) | $10\text{ GPa}$ |
| Coefficient de Poisson Horizontal ($\nu_{hh}$) | $0.25$ |
| Coefficient de Poisson Vertical-Horizontal ($\nu_{vh}$) | $0.20$ |
| Module de Cisaillement Vertical ($G_{vh}$) | $8\text{ GPa}$ |

Cet échantillon sera soumis à de nouveaux essais triaxiaux avec contraintes principales $\sigma_{11} = 10 \text{ MPa}$, $\sigma_{22} = 5 \text{ MPa}$, $\sigma_{33} = 15 \text{ MPa}$ et $\sigma_{12} = \sigma_{23} = \sigma_{31} = 0$. 
"""

# ╔═╡ ea9cf99c-63df-47fe-bc76-1a8d74ce7d24
begin
    # Paramètres élastiques TI
    E_h = 30_000.0     # MPa
    E_v = 10_000.0     # MPa
    nu_hh = 0.25
    nu_vh = 0.25
    G_vh = 8_000.0     # MPa
    G_hh = E_h / (2*(1 + nu_hh))

    # Contraintes appliquées
    sigma_ti = [10.0, 5.0, 15.0, 0.0, 0.0, 0.0]

    # Initialisation de la matrice de compliance
    S_ti = zeros(6, 6)

    # Composantes normales
    S_ti[1,1] = 1/E_h;  S_ti[1,2] = -nu_hh/E_h;  S_ti[1,3] = -nu_vh/E_v
    S_ti[2,1] = -nu_hh/E_h;  S_ti[2,2] = 1/E_h;  S_ti[2,3] = -nu_vh/E_v
    S_ti[3,1] = -nu_vh/E_v;  S_ti[3,2] = -nu_vh/E_v;  S_ti[3,3] = 1/E_v

    # Composantes de cisaillement
    S_ti[4,4] = 1/G_vh
    S_ti[5,5] = 1/G_vh
    S_ti[6,6] = 1/G_hh

    # Calcul des déformations
    epsilon_ti = zeros(6)
    for i in 1:6
        for j in 1:6
            epsilon_ti[i] += S_ti[i,j] * sigma_ti[j]
        end
    end
end

# ╔═╡ 5bb957f1-ba8d-4ce3-98ca-09ea9a423870
md"""
## Exercice 05

En considérant la même roche que dans l'exercice précédent, calculez l'indice d'anisotropie élastique.
"""

# ╔═╡ 34f2592c-55bf-4eb0-ad55-9eaec705c41b
begin
    A_e = E_h/E_v
    nothing
end

# ╔═╡ c079af61-8a8c-46fb-ac7d-5575204bf40c
md"""
# Modèles en orientations génériques
Dans l'approche théorique construite jusqu'à présent, nous avons adopté quelques simplifications découlant de l'orientation des systèmes de symétrie de façon à ce qu'elles coïncident avec les directions principales de contrainte. Cependant, dans les problèmes pratiques, il est courant que la tectonique locale ou régionale entraîne des configurations structurelles plus complexes. 

Dans ce sens, il est fondamental de formaliser le modèle dans des systèmes de coordonnées généralisés. Pour cela, nous effectuons des rotations sur les tenseurs afin de les réorienter selon les besoins du massif étudié.

Pour effectuer la transposition des données de terrain (mesures géologiques), la stratégie consiste à partir des mesures d'attitude des couches, réaliser la description en fonction des cosinus directeurs et recalculer la matrice de compliance dans le nouveau système. 

La formalisation la plus courante dans la littérature, du point de vue de la mécanique des roches, est d'adopter la conversion en fonction de l'angle entre l'orientation des axes dans le système original et dans le système tourné, conformément à l'Éq. 15 :

$$\begin{array}{lll}
l_1 = \cos(x, x') & l_2 = \cos(x, y') & l_3 = \cos(x, z') \\
m_1 = \cos(y, x') & m_2 = \cos(y, y') & m_3 = \cos(y, z') \\
n_1 = \cos(z, x') & n_2 = \cos(z, y') & n_3 = \cos(z, z')
\end{array}\quad\quad (Éq.15)$$

Les cosinus directeurs sont alors organisés dans la matrice de transformation $[T]$, qui sera responsable de la réorientation du système. 

$$[T] = \begin{bmatrix} 
(l_1)^2 & (m_1)^2 & (n_1)^2 & 2m_1n_1 & 2n_1l_1 & 2l_1m_1 \\
(l_2)^2 & (m_2)^2 & (n_2)^2 & 2m_2n_2 & 2n_2l_2 & 2l_2m_2 \\
(l_3)^2 & (m_3)^2 & (n_3)^2 & 2m_3n_3 & 2n_3l_3 & 2l_3m_3 \\
l_2l_3 & m_2m_3 & n_2n_3 & m_2n_3 + m_3n_2 & n_2l_3 + n_3l_2 & l_2m_3 + l_3m_2 \\
l_3l_1 & m_3m_1 & n_3n_1 & m_3n_1 + m_1n_3 & n_3l_1 + n_1l_3 & l_3m_1 + l_1m_3 \\
l_1l_2 & m_1m_2 & n_1n_2 & m_1n_2 + m_2n_1 & n_1l_2 + n_2l_1 & l_1m_2 + l_2m_1 
\end{bmatrix}\quad\quad (Éq.16)$$

Et sur cette base, nous calculons la matrice de compliance dans le système secondaire en utilisant l'Éq.17 :

$[S'] = [T][S][T]^T\quad\quad (Éq.17)$

La solution, cependant, est peu pratique pour les problèmes géologiques en général, étant donné que les mesures sur le terrain sont collectées avec des orientations de boussole. Pour cela, nous allons utiliser la formalisation basée sur les angles de **strike** ($\alpha$) et de **pendage** ($\beta$) et adopter les équations analytiques pour le calcul des cosinus directeurs :

$$\begin{array}{lll}
l_1 = \sin(\alpha) & l_2 = -\cos(\alpha) & l_3 = 0 \\
m_1 = \cos(\beta)\cos(\alpha) & m_2 = \cos(\beta)\sin(\alpha) & m_3 = -\sin(\beta) \\
n_1 = \sin(\beta)\cos(\alpha) & n_2 = \sin(\beta)\sin(\alpha) & n_3 = \cos(\beta)
\end{array}\quad\quad (Éq.18)$$
"""

# ╔═╡ ccd56e88-96b6-4e74-bc9e-0e41b3273088
md"""
## Exercice 06
Considérez une roche provenant du même affleurement que les échantillons de l'**Exercice 04**, et pour les besoins de l'analyse nous supposerons que les paramètres élastiques restent les mêmes. Cependant, dans cette portion de l'affleurement, des variations structurales sont observées et la foliation présente une orientation 30/120 (pendage/direction de pendage). La roche sera soumise à un essai triaxial dont les contraintes dans le système global sont données par $\sigma_{11} = 10 \text{ MPa}$, $\sigma_{22} = 5 \text{ MPa}$, $\sigma_{33} = 15 \text{ MPa}$.

Calculez :
* a) la matrice de compliance du matériau dans le système global, en considérant la rotation associée à l’orientation de la foliation ;
* b) le tenseur des déformations dans le système global ;
"""

# ╔═╡ f9e4626d-f4bd-4916-98ef-956836082321
md"""### Question A """

# ╔═╡ f081361b-d177-48ed-91e9-85d4e40e450f
# Conversions angulaires
begin
    dip_deg = 30.0
    dip_dir_deg = 120.0
    
    strike_deg = dip_dir_deg - 90.0 # conversion en strike
    α = deg2rad(strike_deg) # conversion en radians
    β = deg2rad(dip_deg)
    
    nothing
end

# ╔═╡ dbee0666-71cb-4f8d-981f-679e5529296b
# Calcul des cosinus directeurs
begin
    # Axe x local (Strike)
    l = [sin(α), -cos(α), 0.0] 
    # Axe y local (Ligne de plongement)
    m = [cos(β)*cos(α), cos(β)*sin(α), -sin(β)]
    # Axe z local (Normal à la couche)
    n = [sin(β)*cos(α), sin(β)*sin(α), cos(β)]
    
    # Attribution des composants individuels pour la matrice T
    l1, l2, l3 = l
    m1, m2, m3 = m
    n1, n2, n3 = n
    
    nothing
end

# ╔═╡ eeb79a78-c9dc-4ab6-b652-a91416ea75be
# Calcul de la matrice de transformation et de la compliance globale
begin
    # Matrice S Locale (Système d'axes 3 normal au plan)
    S_local = [
         1/E_h       -nu_hh/E_h   -nu_vh/E_v    0          0          0;
        -nu_hh/E_h    1/E_h       -nu_vh/E_v    0          0          0;
        -nu_vh/E_v   -nu_vh/E_v    1/E_v        0          0          0;
          0            0            0         1/G_vh       0          0;
          0            0            0           0        1/G_vh       0;
          0            0            0           0          0        1/G_hh
    ]
    
    # Matrice de transformation T 
    T = [
        (l1)^2   (m1)^2   (n1)^2   2*m1*n1         2*n1*l1         2*l1*m1;
        (l2)^2   (m2)^2   (n2)^2   2*m2*n2         2*n2*l2         2*l2*m2;
        (l3)^2   (m3)^2   (n3)^2   2*m3*n3         2*n3*l3         2*l3*m3;
        l2*l3    m2*m3    n2*n3    (m2*n3+m3*n2)   (n2*l3+n3*l2)   (l2*m3+l3*m2);
        l3*l1    m3*m1    n3*n1    (m3*n1+m1*n3)   (n3*l1+n1*l3)   (l3*m1+l1*m3);
        l1*l2    m1*m2    n1*n2    (m1*n2+m2*n1)   (n1*l2+n2*l1)   (l1*m2+l2*m1)
    ]
    
    # Matrice de compliance globale (S')
    S_global = T * S_local * T'
    nothing
end

# ╔═╡ 4cbc7429-75c9-4d4d-aa1a-585d6a0f5163
md"""### Question B"""

# ╔═╡ f3078e9d-ba6f-4ad4-a683-34e0cae10d54
begin
    sigma_global = [0.0, 0.0, 10.0, 0.0, 0.0, 0.0] # MPa
    epsilon_global = S_global * sigma_global
    nothing
end

# ╔═╡ f5c0608b-1150-4e6f-bbb1-18dff2f0e4d5
md"""
## Exercice 07

Considérez un échantillon de roche avec un plan de discontinuité bien défini qui permet de l’approximer comme un matériau à comportement transversalement isotrope. Cet échantillon sera soumis à une compression uniaxiale de 10 GPa. 

Des informations préliminaires indiquent que cette roche possède les paramètres élastiques suivants :

- Module vertical ($E_v$) : $15\text{ GPa}$
- Module horizontal ($E_h$) : $25\text{ GPa}$
- Coefficients de Poisson : $\nu_{vh} = 0.30$ et $\nu_{hh} = 0.25$
- Module de cisaillement vertical ($G_{vh}$) : $7\text{ GPa}$

Sur cette base, répondez aux questions suivantes :

**A)** Calculez la matrice de compliance $[S]$ transformée en considérant que la roche possède une orientation des discontinuités 220/30 (pendage/direction de pendage) ;

**B)** Construisez le graphique de la variation du module de Young vertical ($S_{33}$) en fonction de l’angle de la discontinuité ; 

**C)** Réalisez un graphique montrant la variation de la déformation verticale en fonction de l’angle de la discontinuité ;
"""

# ╔═╡ c1cbf5e0-3dd4-449f-af7f-4c48502f840c
md"""### Question A"""

# ╔═╡ d3a30c94-e288-42dc-a14c-9ef0ce26df98
begin
	# 1. Définition des propriétés élastiques avec les nouveaux noms
	# Unités en GPa
	E_v2 = 15.0
	E_h2 = 25.0
	nu_vh2 = 0.30  # Coefficient de Poisson pour une charge verticale
	nu_hh2 = 0.25  # Coefficient de Poisson pour une charge horizontale dans le plan d’isotropie
	G_vh2 = 7.0    # Module de cisaillement dans le plan vertical
	
	# Module de cisaillement dans le plan d’isotropie (G_hh2)
	G_hh2 = E_h2 / (2 * (1 + nu_hh2))
	
	# 2. Construction de la matrice de compliance locale (S_local_nouvelle)
	# L’axe 3 est l’axe de symétrie (normal à la stratification)
	S_local_nouvelle = [
		 1/E_h2       -nu_hh2/E_h2  -nu_vh2/E_v2   0.0          0.0          0.0;
		-nu_hh2/E_h2   1/E_h2       -nu_vh2/E_v2   0.0          0.0          0.0;
		-nu_vh2/E_v2  -nu_vh2/E_v2   1/E_v2        0.0          0.0          0.0;
		 0.0           0.0           0.0           1/G_vh2      0.0          0.0;
		 0.0           0.0           0.0           0.0          1/G_vh2      0.0;
		 0.0           0.0           0.0           0.0          0.0          1/G_hh2
	]
	
	# 3. Définition de l’attitude de la discontinuité (pendage/direction de pendage)
	# Pour un point spécifique (comme dans l’énoncé) :
	pendage_fixe = 30.0            # Pendage
	direction_pendage_fixe = 220.0 # Direction de pendage
	
	beta_rad_roche = deg2rad(pendage_fixe)
	alpha_rad_roche = deg2rad(direction_pendage_fixe)
	
	nothing
end

# ╔═╡ c2d00603-e869-4109-a730-87183db9a9c6
md"""### Question B"""

# ╔═╡ 4e7dd0cb-a2df-460f-a21f-dc1d097e4281
function calculer_reponse_ti(beta_deg, alpha_rad, S_local)
    beta_rad = deg2rad(beta_deg)
    
    # Cosinus directeurs (Attention à la convention des axes !)
    l1, l2, l3 = sin(alpha_rad), -cos(alpha_rad), 0.0
    m1, m2, m3 = cos(beta_rad)*cos(alpha_rad), cos(beta_rad)*sin(alpha_rad), -sin(beta_rad)
    n1, n2, n3 = sin(beta_rad)*cos(alpha_rad), sin(beta_rad)*sin(alpha_rad), cos(beta_rad)
    
    # Matrice de transformation T (Transformation de Bond pour contraintes/compliance)
    # Notez que pour la compliance S' = T' * S * T ou similaire selon la définition de T
    T = [
        (l1)^2   (m1)^2   (n1)^2   2*m1*n1         2*n1*l1         2*l1*m1;
        (l2)^2   (m2)^2   (n2)^2   2*m2*n2         2*n2*l2         2*l2*m2;
        (l3)^2   (m3)^2   (n3)^2   2*m3*n3         2*n3*l3         2*l3*m3;
        l2*l3    m2*m3    n2*n3    (m2*n3+m3*n2)   (n2*l3+n3*l2)   (l2*m3+l3*m2);
        l3*l1    m3*m1    n3*n1    (m3*n1+m1*n3)   (n3*l1+n1*l3)   (l3*m1+l1*m3);
        l1*l2    m1*m2    n1*n2    (m1*n2+m2*n1)   (n1*l2+n2*l1)   (l1*m2+l2*m1)
    ]
    
    # S_global pour la compliance (Transformation de Bond)
    S_global = T * S_local * T'
    
    # Nous retournons l’inverse du terme 3,3 qui représente le module de Young vertical effectif à l’angle beta
    return 1 / S_global[3, 3] 
end

# ╔═╡ 805e087a-a4de-4de4-949b-2afc3db7533e
begin
    plage_beta = 0.0:1.0:90.0
    
    # Notez que nous passons maintenant S_local_nouvelle et alpha_rad_roche
    resultats_courbe = [calculer_reponse_ti(b, alpha_rad_roche, S_local_nouvelle) for b in plage_beta]
    
    plot(plage_beta, resultats_courbe, 
        title="Module de Young vertical avec le pendage",
        label="S_{33} (GPa)", 
        xlabel="Angle de pendage (beta) [degrés]", 
        ylabel="E [GPa]",
        lw=3, 
        grid=true)
end

# ╔═╡ d891b81d-82ea-4139-b228-deb0966c11b8
md"""### Question C"""

# ╔═╡ a980c0ec-10f6-497f-b9f1-3b77580d0eec
begin
	# 1. Contrainte verticale
	sigma_z_valeur = 0.010 
	
	# 2. Calcul de la déformation verticale pour chaque angle
	deformation_verticale = [ (1/r) * sigma_z_valeur for r in resultats_courbe ]
	
	# 3. Tracé direct
	plot(plage_beta, deformation_verticale,
		title="Réponse géomécanique : Déformation vs Pendage",
		label="Déformation verticale (ε₃₃)",
		xlabel="Angle de pendage (β) [degrés]",
		ylabel="Strain [sans dimension]",
		lw=3,
		lc=:darkgreen,
		grid=true)
end

# ╔═╡ 3e84285a-8629-4f31-960d-9918abc9064a
md"""
# Modèle isotrope
Dans la discussion menée jusqu'à présent, nous sommes partis d'un matériau purement anisotrope et des axes de symétrie ont été progressivement ajoutés. Le modèle isotrope se trouve à l'extrême opposé, où toutes les propriétés sont invariantes dans toutes les directions. 

En étendant les interprétations de symétrie du modèle transversalement isotrope, nous pouvons penser, de manière didactique, que dans ce cas, au lieu de la symétrie en $x_3$ observée précédemment, maintenant tous les axes se comportent de cette manière. Autrement dit, nous pouvons supposer :

* Modules normaux égaux dans toutes les directions ($C_{11}=C_{33}$);
* Accouplements normaux égaux ($C_{12}=C_{13}$);
* Cisaillements équivalents ($C_{44}=C_{55}=C_{66}$);

Ainsi, l'imposition d'un nombre infini de plans de symétrie réduit drastiquement la complexité du tenseur de rigidité. Alors que le modèle orthotrope nécessitait 9 constantes indépendantes et le modèle transversalement isotrope 5, le modèle isotrope requiert seulement 2 constantes indépendantes pour décrire complètement le comportement élastique du matériau.

Cela implique que la propriété de cisaillement n'est plus indépendante, mais une fonction directe des modules normaux et des accouplements. Mathématiquement, pour que la matrice soit invariante sous toute rotation de coordonnées, il faut respecter la relation :

$C_{44}=C_{55}=C_{66} = \frac{1}{2}(C_{11} - C_{12})\quad\quad (Eq.19)$

Cette restriction garantit que le matériau répond de la même manière, que la contrainte soit appliquée selon les axes principaux ou dans une direction oblique quelconque.

Ainsi, la matrice de rigidité mise à jour pour le système isotrope devient :

$$[C] = \begin{bmatrix} 
C_{11} & C_{12} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{11} & C_{12} & 0 & 0 & 0 \\
C_{12} & C_{12} & C_{11} & 0 & 0 & 0 \\
0 & 0 & 0 & \frac{C_{11}-C_{12}}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & \frac{C_{11}-C_{12}}{2} & 0 \\
0 & 0 & 0 & 0 & 0 & \frac{C_{11}-C_{12}}{2}
\end{bmatrix}\quad\quad (Eq.20)$$

Ainsi, la relation constitutive de Hooke complète, en fonction de la matrice de rigidité, pour un milieu isotrope, linéaire et élastique utilisant les paramètres d’ingénierie devient :

$$\begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \tau_{23} \\ \tau_{13} \\ \tau_{12} \end{bmatrix} = \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix} 
1-\nu & \nu & \nu & 0 & 0 & 0 \\
\nu & 1-\nu & \nu & 0 & 0 & 0 \\
\nu & \nu & 1-\nu & 0 & 0 & 0 \\
0 & 0 & 0 & \frac{1-2\nu}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & \frac{1-2\nu}{2} & 0 \\
0 & 0 & 0 & 0 & 0 & \frac{1-2\nu}{2}
\end{bmatrix}
\begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ \gamma_{23} \\ \gamma_{13} \\ \gamma_{12} \end{bmatrix}\quad (Eq.21)$$

La relation devient considérablement plus simple si nous utilisons la matrice de compliance :

$$\begin{bmatrix} 
\epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ \gamma_{23} \\ \gamma_{31} \\ \gamma_{12} 
\end{bmatrix} = 
\frac{1}{E} \begin{bmatrix} 
1 & -\nu & -\nu & 0 & 0 & 0 \\ 
-\nu & 1 & -\nu & 0 & 0 & 0 \\ 
-\nu & -\nu & 1 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 2(1+\nu) & 0 & 0 \\ 
0 & 0 & 0 & 0 & 2(1+\nu) & 0 \\ 
0 & 0 & 0 & 0 & 0 & 2(1+\nu) 
\end{bmatrix} 
\begin{bmatrix} 
\sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{31} \\ \sigma_{12} 
\end{bmatrix}\quad\quad (Eq.22)$$
"""

# ╔═╡ 4f236e1e-4857-4c1b-a3a3-e0d25bc02c75
md"""
## Exercice 07

Considérez un massif rocheux en profondeur ayant un comportement linéaire, élastique et isotrope. L'état de contraintes _in situ_ a pour directions principales de contrainte : $\sigma_{ii}=(15, 20, 30)$. Les paramètres élastiques du massif rocheux sont $E= 35 \text{ GPa}$ et $\nu=0.25$. Calculez le tenseur des déformations. 
"""

# ╔═╡ db6ea286-b90a-473b-8689-578e82f59aa1
md"""
**Remarque :** Le langage Julia permet la déclaration de matrices de manière très directe, et cette approche a été utilisée dans les exercices précédents. Ici, cependant, pour la simplicité de la matrice, nous allons développer les calculs en émulant le processus manuel de multiplication de matrices. L'objectif est purement didactique. 
"""

# ╔═╡ d57de5a6-82bc-4713-a06b-e91fc29bb49a
begin
	# Insertion des paramètres
    E = 35000.0
    nu = 0.25
    sigma = [15.0, 20.0, 30.0, 0.0, 0.0, 0.0]

	# Création d'une matrice nulle pour recevoir les valeurs
    S = zeros(6, 6)

	# Insertion des informations des composantes normales
    for i in 1:3
        for j in 1:3
            S[i, j] = (i == j) ? 1/E : -nu/E
        end
    end

	# Insertion des informations des composantes de cisaillement
    G_inv = 2 * (1 + nu) / E
    for i in 4:6
        S[i, i] = G_inv
    end

	# Boucle de multiplication matrice-vecteur
    epsilon = zeros(6)
    for i in 1:6
        for j in 1:6
            epsilon[i] += S[i, j] * sigma[j]
        end
    end
end

# ╔═╡ aef76a33-3cdb-4a6b-965e-ada9f468e898
md"""
## Exercice 08

Considérez un tunnel linéaire étendu, creusé dans un massif rocheux isotrope, linéaire et élastique. Admettez un état de déformations planes avec l’axe orienté le long de la direction "y". L'état de contraintes _in situ_ dans le plan transversal est donné par $(\sigma_{xx}, \sigma_{zz})=(15,20)$, en MPa. Déterminez $\sigma_{yy}$ imposée par la déformation plane. 
"""

# ╔═╡ 9b4bf833-99dd-418c-9db1-455213680e4c
begin
	# Définition des paramètres d'entrée
	σ_xx = 15.0  # MPa
	σ_zz = 20.0  # MPa
	ν = 0.25      # Coefficient de Poisson (exemple : valeur courante pour les roches)
	
	# Application directe de la déduction : σ_yy = ν(σ_xx + σ_zz)
	σ_yy = ν * (σ_xx + σ_zz)
	nothing
end

# ╔═╡ efd223aa-0faf-4d85-9c02-cef7482b5cf9
md"""
## Exercice 09
Une roche présente trois plans de symétrie en raison du stratifié et d’un système de fractures verticales. Les attitudes des familles de discontinuités (dip-direction / dip) dans le système local de coordonnées sont données par :

- Plan 1 : 120° / 60°
- Plan 2 : 30° / 60°
- Plan 3 : 210° / 30°

Les propriétés de la roche sont les suivantes :
* Modules (GPa) : $E_1 = 20$, $E_2 = 25$, $E_3 = 15$
* Poisson : $\nu_{12}=0,20$, $\nu_{13}=0,25$, $\nu_{23}=0,30$
* Cisaillement (GPa) : $G_{12}=10$, $G_{13}=7$, $G_{23}=8$

L’état de contraintes global est simplifié et donné par la matrice suivante en MPa :

$$\sigma_{global} = \begin{bmatrix} 40 & 0 & 0 \\ 0 & 30 & 0 \\ 0 & 0 & 20 \end{bmatrix}$$

À partir de ces informations, calculez le tenseur des déformations en utilisant les trois modèles étudiés (isotrope, transversalement isotrope et orthotrope). Comparez les résultats et commentez l’effet du couplage mécanique.  

_Remarque : pour les modèles TI et isotrope, utilisez les moyennes arithmétiques des propriétés ci-dessus conformément à la symétrie de chaque système._
"""

# ╔═╡ d0893a8f-1d26-4e54-97c1-6d9a81b1c515
begin
	# 1. Paramètres orthotropes
	E1, E2, E3 = 20.0, 25.0, 15.0
	nu12, nu13, nu23 = 0.20, 0.25, 0.30
	G12, G13, G23 = 10.0, 7.0, 8.0
	
	# Coefficients de Poisson complémentaires
	nu21 = nu12 * (E2/E1)
	nu31 = nu13 * (E3/E1)
	nu32 = nu23 * (E3/E2)

	# 2. Matrice de compliance orthotrope locale
	S_ortho_locale = [
		 1/E1      -nu21/E2   -nu31/E3    0.0        0.0        0.0;
		-nu12/E1    1/E2      -nu32/E3    0.0        0.0        0.0;
		-nu13/E1   -nu23/E2    1/E3       0.0        0.0        0.0;
		 0.0        0.0        0.0        1/G23      0.0        0.0;
		 0.0        0.0        0.0        0.0        1/G13      0.0;
		 0.0        0.0        0.0        0.0        0.0        1/G12
	]

	# 3. Orientation des plans (dip-direction / dip)
	plans = [
		(120.0, 60.0),
		(30.0, 60.0),
		(210.0, 30.0)
	]

	# 4. Contraintes globales (MPa → GPa)
	sigma_g = [0.040, 0.030, 0.020, 0.0, 0.0, 0.0]

	nothing
end

# ╔═╡ 67bbd68f-159c-451a-84df-2fb2fbe77cc4
begin
	# Modèle orthotrope
	S_ortho = S_ortho_locale 

	# 4. Modèle transversal isotrope (TI)
	E_h_ti = (E1 + E2) / 2
	E_v_ti = E3

	nu_hh_ti = nu12
	nu_hv_ti = (nu13 + nu23) / 2
	nu_vh_ti = nu_hv_ti * (E_v_ti/E_h_ti)

	G_v_ti = (G13 + G23) / 2
	G_h_ti = E_h_ti / (2 * (1 + nu_hh_ti))

	S_ti_locale = [
		 1/E_h_ti      -nu_hh_ti/E_h_ti   -nu_hv_ti/E_h_ti  0.0  0.0  0.0;
		-nu_hh_ti/E_h_ti  1/E_h_ti        -nu_hv_ti/E_h_ti  0.0  0.0  0.0;
		-nu_vh_ti/E_v_ti -nu_vh_ti/E_v_ti  1/E_v_ti        0.0  0.0  0.0;
		 0.0  0.0  0.0  1/G_v_ti  0.0  0.0;
		 0.0  0.0  0.0  0.0  1/G_v_ti  0.0;
		 0.0  0.0  0.0  0.0  0.0  1/G_h_ti
	]

	# 5. Modèle isotrope
	E_iso = (E1 + E2 + E3) / 3
	nu_iso = (nu12 + nu13 + nu23) / 3
	G_iso = E_iso / (2 * (1 + nu_iso))

	S_iso_locale = [
		 1/E_iso      -nu_iso/E_iso -nu_iso/E_iso  0.0  0.0  0.0;
		-nu_iso/E_iso  1/E_iso      -nu_iso/E_iso  0.0  0.0  0.0;
		-nu_iso/E_iso -nu_iso/E_iso  1/E_iso       0.0  0.0  0.0;
		 0.0  0.0  0.0  1/G_iso  0.0  0.0;
		 0.0  0.0  0.0  0.0  1/G_iso  0.0;
		 0.0  0.0  0.0  0.0  0.0  1/G_iso
	]

	nothing
end

# ╔═╡ f8a43654-5385-47bc-b285-d442e4858635
begin
	function calculer_deformations_modeles(dir1_deg, dip1_deg, dir2_deg, dip2_deg, dir3_deg, dip3_deg, S_loc, sigma_g)

	    α1, β1 = deg2rad(dir1_deg), deg2rad(dip1_deg)
	    α2, β2 = deg2rad(dir2_deg), deg2rad(dip2_deg)
	    α3, β3 = deg2rad(dir3_deg), deg2rad(dip3_deg)

	    # Cosinus directeurs (normales aux plans)
	    l1 = sin(β1)*sin(α1); m1 = sin(β1)*cos(α1); n1 = cos(β1)
	    l2 = sin(β2)*sin(α2); m2 = sin(β2)*cos(α2); n2 = cos(β2)
	    l3 = sin(β3)*sin(α3); m3 = sin(β3)*cos(α3); n3 = cos(β3)

	    T = [
	        l1^2 m1^2 n1^2 2*m1*n1 2*n1*l1 2*l1*m1;
	        l2^2 m2^2 n2^2 2*m2*n2 2*n2*l2 2*l2*m2;
	        l3^2 m3^2 n3^2 2*m3*n3 2*n3*l3 2*l3*m3;
	        l2*l3 m2*m3 n2*n3 (m2*n3+m3*n2) (n2*l3+n3*l2) (l2*m3+l3*m2);
	        l3*l1 m3*m1 n3*n1 (m3*n1+m1*n3) (n3*l1+n1*l3) (l3*m1+l1*m3);
	        l1*l2 m1*m2 n1*n2 (m1*n2+m2*n1) (n1*l2+n2*l1) (l1*m2+l2*m1)
	    ]

	    S_global = T * S_loc * T'
	    return S_global * sigma_g
	end


	begin
		# Plans de symétrie (dip-direction / dip)
		dir1, dip1 = 120, 60
		dir2, dip2 = 30, 60
		dir3, dip3 = 210, 30

		# Calcul pour les 3 modèles
		def_ortho = calculer_deformations_modeles(dir1, dip1, dir2, dip2, dir3, dip3, S_ortho, sigma_g)

		def_ti    = calculer_deformations_modeles(dir1, dip1, dir2, dip2, dir3, dip3, S_ti_locale, sigma_g)

		# isotrope : rotation inutile
		def_iso   = S_iso_locale * sigma_g
		
		nothing
	end
end

# ╔═╡ 66abddf0-a5dc-4d14-a384-defa912b7dfe
md"""
> **Note :** La rotation des plans inclinés fait que les contraintes normales globales agissent obliquement, générant un cisaillement local et donc des γ ≠ 0, même avec seulement des contraintes normales dans le système global.
"""

# ╔═╡ 0a3b67c5-5053-4751-99eb-ad46695f9d63
md"""

## Conclusions

Les résultats démontrent l’impact de la prise en compte de la variabilité directionnelle dans les modèles mécaniques constitutifs. Dans un contexte géologique, cette variabilité souligne la nécessité de représenter les éléments structuraux comme partie intégrante des modèles mécaniques.

En comparant, par exemple, le modèle isotrope et le modèle TI, on observe que même en considérant uniquement le plan d’isotropie, il existe des divergences dans les déformations. Cela s’explique par le fait que le modèle TI distingue $E_v$ et $E_h$ et prend en compte l’effet combiné de ces deux paramètres, qui est “homogénéisé” par le modèle isotrope.

De manière générale, on peut considérer que, d’une part, dans les modèles isotropes, une contrainte normale entraîne exclusivement une déformation normale, et la réponse aux contraintes de cisaillement est directe. Dans le cas des modèles TI, une contrainte normale peut induire des déformations de cisaillement résultant précisément de ce comportement relatif entre les paramètres élastiques dans le plan d’isotropie et perpendiculairement à celui-ci. Cet effet est appelé couplage constitutif.

Dans le modèle orthotrope, le couplage devient encore plus prononcé en raison de l’existence de plusieurs plans de symétrie. Dans ce cas, le facteur directionnel du rake et l’ajout de paramètres élastiques supplémentaires rendent le modèle plus complexe.
"""

# ╔═╡ c147a4da-344f-4dc0-886f-ca4d4198b9f8
md"""
# Littérature recommandée

[1] Jaeger, J. C., Cook, N. G. W., & Zimmerman, R. W. (2007). Linear elasticity. In *Fundamentals of rock mechanics* (4th ed., pp. 106–144). Wiley-Blackwell.

[2] Jing, L., & Stephansson, O. (2007). Constitutive models of rock fractures and rock masses – The basic. In *Fundamentals of discrete element methods for rock engineering: Theory and applications* (pp. 47–109). Elsevier Science.

[3] Malvern, L. E. (1977). *Introduction to the mechanics of a continuous medium*. Pearson.
"""

# ╔═╡ 77d669f4-64f0-4105-9963-ba8526b7f2ad
md"""##### Fonctions de support"""

# ╔═╡ 9f2d8af0-9d04-4ffd-b78f-d1716e1ccdcb
begin
	hint(text) = (Markdown.MD(Markdown.Admonition("hint", 
		"Astuce", [text])))
	nothing
end

# ╔═╡ 8d9d75ea-e451-41cc-adbe-f9b078647ec6
hint(md"""
Pour les matériaux transversalement isotropes, $G_{hh} = \frac{E_h}{2(1 +\nu_{hh})}$
""")

# ╔═╡ 531f17d4-6fc4-4029-b34f-0164955453b0
hint(md"""
Pour les matériaux transversalement isotropes, l'indice d'anisotropie $A_E$ est donné par : $A_E = \frac{E_v}{E_h}$
""")

# ╔═╡ a7433aef-0954-4f60-afea-87a95062f1fd
hint(md"""
As equações apresentadas anteriormente esperam ângulos em _strike/dip_, sendo assim, é necessário realizar a conversão: 

$strike=dip-direction−90\degree$
""")

# ╔═╡ 6a5549da-1610-46df-ac8d-9360604cec0f
hint(md"""
Pour la solution, nous allons considérer l'Éq.10 comme référence, ce qui donne :

$$\begin{bmatrix} 
\epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ \gamma_{23} \\ \gamma_{31} \\ \gamma_{12} 
\end{bmatrix} = 
\frac{1}{35000} \begin{bmatrix} 
1 & -0.25 & -0.25 & 0 & 0 & 0 \\ 
-0.25 & 1 & -0.25 & 0 & 0 & 0 \\ 
-0.25 & -0.25 & 1 & 0 & 0 & 0 \\ 
0 & 0 & 0 & 2.5 & 0 & 0 \\ 
0 & 0 & 0 & 0 & 2.5 & 0 \\ 
0 & 0 & 0 & 0 & 0 & 2.5 
\end{bmatrix} 
\begin{bmatrix} 
15 \\ 20 \\ 30 \\ 0 \\ 0 \\ 0 
\end{bmatrix}$$

N'oubliez pas de convertir les unités !
""")

# ╔═╡ ac30cba0-b35e-4adb-aab6-a0b8f857e1ec
hint(md"""
1) Le problème établit un état de déformation plane, c’est-à-dire que le problème est invariant le long de "y" et $\epsilon_{yy} = 0$. 

2) En développant l'Éq.10, on obtient : 

$\varepsilon_{yy} = \frac{1}{E} \left( \sigma_{yy} - \nu (\sigma_{xx} + \sigma_{zz}) \right)$

En substituant $\epsilon_{yy} = 0$ et en isolant la variable d’intérêt, on obtient :

$\sigma_{yy}= \nu(\sigma_{xx} + \sigma_{zz})$
""")

# ╔═╡ 0cf577f7-f22e-4ada-9000-11147813ac8e
begin
	ans(text) = (Markdown.MD(Markdown.Admonition("ans", 
		"Réponse", [text])))
	nothing
end

# ╔═╡ 51851d8a-63ca-41cb-9548-6a5559f7f1b2
ans(md"""
Dans les matériaux élastiques, l’énergie est stockée lors du chargement et libérée lors du déchargement. Ainsi, elle peut être décrite comme un **potentiel scalaire** $\bar{U}$ dont la dérivée donne la contrainte $\tau_{ij}$.

$$\tau_{ij} = \frac{\partial \bar{U}}{\partial \epsilon_{ij}}$$
	
La matrice constitutive relie la contrainte et la déformation ; elle peut donc être décrite comme la dérivée seconde de l’énergie accumulée :
	
$$C_{ijkl} = \frac{\partial \tau_{ij}}{\partial \epsilon_{kl}} = \frac{\partial^2 \bar{U}}{\partial \epsilon_{kl} \partial \epsilon_{ij}}$$

Selon le **théorème de Schwarz**, l’ordre des dérivées partielles n’affecte pas le résultat, c’est-à-dire :

$$\frac{\partial^2 \bar{U}}{\partial \epsilon_{kl} \partial \epsilon_{ij}} = \frac{\partial^2 \bar{U}}{\partial \epsilon_{ij} \partial \epsilon_{kl}} \therefore C_{ijkl} = C_{klij}$$
""")

# ╔═╡ 154a5c9c-a793-47e1-a8b0-03c4fe72c75d
ans(md"""
En termes pratiques, trois limitations principales peuvent être observées :

* **Limitations procédurales :** les équipements et capteurs disponibles, ainsi que la préparation d’échantillons représentatifs de toutes les directions d’anisotropie, sont difficiles voire irréalisables en pratique ;

* **Couplage des déformations :** dans les matériaux anisotropes, l’application d’une traction uniaxiale produit à la fois des déformations normales et des distorsions (associées au cisaillement). Isoler une composante spécifique de $C_{ijkl}$ devient donc difficile et expérimentalement instable ;

* **Conservation de l’énergie :** la prémisse fondamentale de la méthodologie repose sur la conservation de l’énergie élastique. Cependant, en pratique, la dissipation d’énergie sous d’autres formes, la formation de microfissures ou encore des effets viscoélastiques peuvent violer cette hypothèse.
""")

# ╔═╡ 0d598744-5bc2-4189-bb89-0552004c9f47
ans(md"""
On s’attend à ce que $E_3 < E_2$ ($E_i$ est inversement proportionnel à l’intensité de fracturation). L’énoncé indique que $F2$ est considérablement plus fréquent que $F1$, ce qui implique une plus grande déformabilité associée à un glissement relatif plus important entre les blocs.
""")

# ╔═╡ ff2daffc-8939-4e66-a0a0-8e0e60a2e874
ans(md"""
Non. L’Eq.11 présente la matrice de compliance d’un matériau orthotrope, dans laquelle la relation suivante est imposée :

$\frac{\nu_{23}}{E_2} = \frac{\nu_{32}}{E_3}$

Comme nous l’avons vu dans la question précédente, la configuration structurale du massif impose que $E_3 < E_2$, donc $\nu_{23} \neq \nu_{32}$.
""")

# ╔═╡ df890513-be63-405e-80ce-ca95c9e81322
ans(md"""
En considérant la relation constitutive orthotrope uniquement pour les directions principales, nous avons :

$$\begin{bmatrix} \epsilon_1 \\ \epsilon_2 \\ \epsilon_3 \end{bmatrix} = 
\begin{bmatrix} 
1/E_1 & -\nu_{21}/E_2 & -\nu_{31}/E_3 \\ 
-\nu_{12}/E_1 & 1/E_2 & -\nu_{32}/E_3 \\ 
-\nu_{13}/E_1 & -\nu_{23}/E_2 & 1/E_3 
\end{bmatrix} 
\begin{bmatrix} \sigma_1 \\ \sigma_2 \\ \sigma_3 \end{bmatrix}$$

Si l’on considère une compression uniaxiale, nous avons :

$\sigma_1 = 0, \quad \sigma_2 = 0 \quad et \quad \sigma_3 = \sigma$

Ainsi, la multiplication matricielle se résume à :
	
$$\epsilon_3 = \left( -\frac{\nu_{13}}{E_1} \cdot 0 \right) + \left( -\frac{\nu_{23}}{E_2} \cdot 0 \right) + \left( \frac{1}{E_3} \cdot \sigma \right) \implies \epsilon_3 = \frac{\sigma}{E_3}$$

Et de manière analogue :

$\epsilon_2 = \frac{\sigma}{E_2}$

D’après l’item A, nous savons que $E_3 < E_2$, donc :

$\epsilon_3 = \frac{\sigma}{E_3} \quad > \quad \epsilon_2 = \frac{\sigma}{E_2}$

Ainsi, nous pouvons conclure qu’une déformation axiale plus grande est attendue dans la direction $x_3$ en raison de la réduction de la rigidité effective.
""")

# ╔═╡ a3a38f66-3b4b-4ba6-bcb7-fcce3ae5431c
ans(md"""
Les effets dynamiques, tels que la propagation sismique et le régime transitoire expérimenté par la roche, rendent ce problème plus complexe que la forme purement considérée ici. Cependant, étant donné la limitation de portée, nous pouvons associer que le plan qui sera préférentiellement activé par le champ de contraintes de la détonation (pas nécessairement le plan de rupture de la roche) sera le plan de rigidité la plus faible, donc $x_3$.
""")

# ╔═╡ 151cb52b-7c3b-465d-adf0-80470394657f
ans(md"""
La condition d’isotropie ou d’isotropie transversale peut être adoptée dans le cas où :
* Un des plans de discontinuité prédomine sur les autres – cas d’isotropie transversale.
* Distribution statistiquement égale de fractures dans toutes les directions avec un comportement mécanique similaire – approximation du cas isotrope.
* Homogénéisation due à l’échelle – *Representative Elementary Volume*.
""")

# ╔═╡ 03b941fc-6368-4115-94c3-b1df3796c5bc
ans(md"""

| Composante | Valeur |
|-----------|-------|
| ε₁₁ | $(round(epsilon_ti[1], digits=7)) |
| ε₂₂ | $(round(epsilon_ti[2], digits=7)) |
| ε₃₃ | $(round(epsilon_ti[3], digits=7)) |
| γ₂₃ | $(round(epsilon_ti[4], digits=7)) |
| γ₃₁ | $(round(epsilon_ti[5], digits=7)) |
| γ₁₂ | $(round(epsilon_ti[6], digits=7)) |
""")

# ╔═╡ cd716300-0821-4b01-b79d-3dd2c835324c
ans(md"""L'indice d'anisotropie est de $A_e.""")

# ╔═╡ 8835568c-dba0-4714-9422-f58d0dad9cc2
ans(md"""
Résultats :

|   | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| **1** | $(round(S_global[1,1], digits=9)) | $(round(S_global[1,2], digits=9)) | $(round(S_global[1,3], digits=9)) | $(round(S_global[1,4], digits=9)) | $(round(S_global[1,5], digits=9)) | $(round(S_global[1,6], digits=9)) |
| **2** | $(round(S_global[2,1], digits=9)) | $(round(S_global[2,2], digits=9)) | $(round(S_global[2,3], digits=9)) | $(round(S_global[2,4], digits=9)) | $(round(S_global[2,5], digits=9)) | $(round(S_global[2,6], digits=9)) |
| **3** | $(round(S_global[3,1], digits=9)) | $(round(S_global[3,2], digits=9)) | $(round(S_global[3,3], digits=9)) | $(round(S_global[3,4], digits=9)) | $(round(S_global[3,5], digits=9)) | $(round(S_global[3,6], digits=9)) |
| **4** | $(round(S_global[4,1], digits=9)) | $(round(S_global[4,2], digits=9)) | $(round(S_global[4,3], digits=9)) | $(round(S_global[4,4], digits=9)) | $(round(S_global[4,5], digits=9)) | $(round(S_global[4,6], digits=9)) |
| **5** | $(round(S_global[5,1], digits=9)) | $(round(S_global[5,2], digits=9)) | $(round(S_global[5,3], digits=9)) | $(round(S_global[5,4], digits=9)) | $(round(S_global[5,5], digits=9)) | $(round(S_global[5,6], digits=9)) |
| **6** | $(round(S_global[6,1], digits=9)) | $(round(S_global[6,2], digits=9)) | $(round(S_global[6,3], digits=9)) | $(round(S_global[6,4], digits=9)) | $(round(S_global[6,5], digits=9)) | $(round(S_global[6,6], digits=9)) |

*Valeurs en $MPa^{-1}$.*
""")

# ╔═╡ 266603d1-7e2d-4d24-b195-3dfc81820440
ans(md"""
Résultats :

| Composante | Valeur |
|:-----------|:-------|
| $\epsilon_{11}$ | $(round(epsilon_global[1], digits=7)) |
| $\epsilon_{22}$ | $(round(epsilon_global[2], digits=7)) |
| $\epsilon_{33}$ | $(round(epsilon_global[3], digits=7)) |
| $\gamma_{23}$   | $(round(epsilon_global[4], digits=7)) |
| $\gamma_{31}$   | $(round(epsilon_global[5], digits=7)) |
| $\gamma_{12}$   | $(round(epsilon_global[6], digits=7)) |
""")

# ╔═╡ 039a6d98-e40f-4ea8-860f-5df6aa48f86b
ans(md"""

| Composante | Valeur |
|-----------|-------|
| ε₁₁ | $(round(epsilon[1], digits=7)) |
| ε₂₂ | $(round(epsilon[2], digits=7)) |
| ε₃₃ | $(round(epsilon[3], digits=7)) |
| γ₂₃ | $(round(epsilon[4], digits=7)) |
| γ₃₁ | $(round(epsilon[5], digits=7)) |
| γ₁₂ | $(round(epsilon[6], digits=7)) |
""")

# ╔═╡ 52d47c1c-9ccb-4e3d-b717-cce533ababcd
ans(md"A contrainte imposée $\sigma_{yy}$ est de $σ_yy \text{ MPa}$.") 

# ╔═╡ 5a5e2246-b3a5-48b0-9ccc-59e5492a5abd
ans(md"""

| Composante | Isotrope | Transv. Isotrope (TI) | Orthotrope (Cible) | Écart (Iso x Ortho) | Écart (TI x Ortho) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **ϵ_xx** | $(round(def_iso[1], sigdigits=4)) | $(round(def_ti[1], sigdigits=4)) | $(round(def_ortho[1], sigdigits=4)) | $(round(((def_ortho[1]-def_iso[1])/def_iso[1])*100, digits=2))% | $(round(((def_ortho[1]-def_ti[1])/def_ti[1])*100, digits=2))% |
| **ϵ_yy** | $(round(def_iso[2], sigdigits=4)) | $(round(def_ti[2], sigdigits=4)) | $(round(def_ortho[2], sigdigits=4)) | $(round(((def_ortho[2]-def_iso[2])/def_iso[2])*100, digits=2))% | $(round(((def_ortho[2]-def_ti[2])/def_ti[2])*100, digits=2))% |
| **ϵ_zz** | $(round(def_iso[3], sigdigits=4)) | $(round(def_ti[3], sigdigits=4)) | $(round(def_ortho[3], sigdigits=4)) | $(round(((def_ortho[3]-def_iso[3])/def_iso[3])*100, digits=2))% | $(round(((def_ortho[3]-def_ti[3])/def_ti[3])*100, digits=2))% |
| **γ_yz** | $(round(def_iso[4], sigdigits=4)) | $(round(def_ti[4], sigdigits=4)) | $(round(def_ortho[4], sigdigits=4)) | -- | -- |
| **γ_xz** | $(round(def_iso[5], sigdigits=4)) | $(round(def_ti[5], sigdigits=4)) | $(round(def_ortho[5], sigdigits=4)) | -- | -- |
| **γ_xy** | $(round(def_iso[6], sigdigits=4)) | $(round(def_ti[6], sigdigits=4)) | $(round(def_ortho[6], sigdigits=4)) | -- | -- |
""")

# ╔═╡ e010405a-88d5-4811-8f5a-193d9434e9ab
begin
	# Summary
	PlutoUI.TableOfContents(aside=true, title="Sommaire",
							indent=true, depth=2)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.41.5"
PlutoUI = "~0.7.79"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.8"
manifest_format = "2.0"
project_hash = "3dc50b702c5bead7f5f043742971bceee102ab5b"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a21c5464519504e41e0cbc91f0188e8ca23d7440"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c9cbeda6aceffc52d8a0017e71db27c7a7c0beaf"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "28145feabf717c5d65c1d5e09747ee7b1ff3ed13"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.3"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e015f211ebb898c8180887012b938f3851e719ac"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.55+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─b1f454c0-1be0-11f1-1707-8dc83ac1495c
# ╟─22a2d086-00f9-43fb-ab60-c6c7df65f0c9
# ╟─27c9ebff-4e97-4a7d-aba0-1c26583d2de8
# ╟─4942766d-0e7c-44fd-9dcd-379e8be02cdb
# ╟─9c03c6b6-63d8-46a3-8c84-8a0e9f424055
# ╟─d7887b31-b1ed-447a-9887-2237574c093f
# ╟─51851d8a-63ca-41cb-9548-6a5559f7f1b2
# ╟─95d91406-842e-442c-9a07-0109252942cc
# ╟─154a5c9c-a793-47e1-a8b0-03c4fe72c75d
# ╟─ad7b4db9-1116-4e26-a861-4e8b42342c28
# ╟─5718c387-4455-4a49-933d-63ffef23fffa
# ╟─0bcd566e-b7dd-40f8-8bb6-ae67db90bc41
# ╟─4bce3618-e8a3-4940-a21e-fcdbd61fc1d2
# ╟─0d598744-5bc2-4189-bb89-0552004c9f47
# ╟─916a9876-b586-46a3-848f-af54e2e92a83
# ╟─ff2daffc-8939-4e66-a0a0-8e0e60a2e874
# ╟─535aaadf-0402-421f-9035-fdb502e20f7b
# ╟─df890513-be63-405e-80ce-ca95c9e81322
# ╟─85f9d2b4-e450-4f8e-bec9-f1ed6f15a2e9
# ╟─a3a38f66-3b4b-4ba6-bcb7-fcce3ae5431c
# ╟─70cc842d-6b86-47d2-a356-5746cd0f7524
# ╟─151cb52b-7c3b-465d-adf0-80470394657f
# ╟─6e83a91d-bba7-4c84-9e18-429719b1ee0f
# ╟─72ba316d-6105-443b-a2b9-6bb8d49c9ece
# ╟─8d9d75ea-e451-41cc-adbe-f9b078647ec6
# ╠═ea9cf99c-63df-47fe-bc76-1a8d74ce7d24
# ╟─03b941fc-6368-4115-94c3-b1df3796c5bc
# ╟─5bb957f1-ba8d-4ce3-98ca-09ea9a423870
# ╟─531f17d4-6fc4-4029-b34f-0164955453b0
# ╠═34f2592c-55bf-4eb0-ad55-9eaec705c41b
# ╟─cd716300-0821-4b01-b79d-3dd2c835324c
# ╟─c079af61-8a8c-46fb-ac7d-5575204bf40c
# ╟─ccd56e88-96b6-4e74-bc9e-0e41b3273088
# ╟─a7433aef-0954-4f60-afea-87a95062f1fd
# ╟─f9e4626d-f4bd-4916-98ef-956836082321
# ╠═f081361b-d177-48ed-91e9-85d4e40e450f
# ╠═dbee0666-71cb-4f8d-981f-679e5529296b
# ╠═eeb79a78-c9dc-4ab6-b652-a91416ea75be
# ╟─8835568c-dba0-4714-9422-f58d0dad9cc2
# ╟─4cbc7429-75c9-4d4d-aa1a-585d6a0f5163
# ╠═f3078e9d-ba6f-4ad4-a683-34e0cae10d54
# ╟─266603d1-7e2d-4d24-b195-3dfc81820440
# ╟─f5c0608b-1150-4e6f-bbb1-18dff2f0e4d5
# ╟─c1cbf5e0-3dd4-449f-af7f-4c48502f840c
# ╠═d3a30c94-e288-42dc-a14c-9ef0ce26df98
# ╟─c2d00603-e869-4109-a730-87183db9a9c6
# ╠═4e7dd0cb-a2df-460f-a21f-dc1d097e4281
# ╠═805e087a-a4de-4de4-949b-2afc3db7533e
# ╟─d891b81d-82ea-4139-b228-deb0966c11b8
# ╠═a980c0ec-10f6-497f-b9f1-3b77580d0eec
# ╟─3e84285a-8629-4f31-960d-9918abc9064a
# ╟─4f236e1e-4857-4c1b-a3a3-e0d25bc02c75
# ╟─6a5549da-1610-46df-ac8d-9360604cec0f
# ╟─db6ea286-b90a-473b-8689-578e82f59aa1
# ╠═d57de5a6-82bc-4713-a06b-e91fc29bb49a
# ╟─039a6d98-e40f-4ea8-860f-5df6aa48f86b
# ╟─aef76a33-3cdb-4a6b-965e-ada9f468e898
# ╟─ac30cba0-b35e-4adb-aab6-a0b8f857e1ec
# ╠═9b4bf833-99dd-418c-9db1-455213680e4c
# ╟─52d47c1c-9ccb-4e3d-b717-cce533ababcd
# ╟─efd223aa-0faf-4d85-9c02-cef7482b5cf9
# ╠═d0893a8f-1d26-4e54-97c1-6d9a81b1c515
# ╠═67bbd68f-159c-451a-84df-2fb2fbe77cc4
# ╠═f8a43654-5385-47bc-b285-d442e4858635
# ╟─5a5e2246-b3a5-48b0-9ccc-59e5492a5abd
# ╟─66abddf0-a5dc-4d14-a384-defa912b7dfe
# ╟─0a3b67c5-5053-4751-99eb-ad46695f9d63
# ╟─c147a4da-344f-4dc0-886f-ca4d4198b9f8
# ╟─77d669f4-64f0-4105-9963-ba8526b7f2ad
# ╠═9f2d8af0-9d04-4ffd-b78f-d1716e1ccdcb
# ╠═0cf577f7-f22e-4ada-9000-11147813ac8e
# ╠═e010405a-88d5-4811-8f5a-193d9434e9ab
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
