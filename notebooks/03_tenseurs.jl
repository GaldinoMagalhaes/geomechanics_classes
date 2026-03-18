### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ b5719bd8-4f04-4d42-a081-7504d946129f
begin
	using PlutoUI
	using CSV
	using DataFrames
	using WGLMakie
	using LinearAlgebra
	using Colors
end

# ╔═╡ a27f367e-5290-4a78-b244-14ba9b9b7c2d
html"""
<div style="display: flex; justify-content: center; background-color: white; padding: 10px; border-radius: 8px;">
	<img src="https://raw.githubusercontent.com/GaldinoMagalhaes/geomechanics_classes/98574ff7d01e1b5c722538b70605a89c3fceddfb/images/logo_GRMM.png" width="800">
</div>
"""

# ╔═╡ c44c680a-1c0a-4fbe-af3a-af512ce272f4
html"""
<p style="background-color:lightgray" xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">&nbsp&nbsp💻&nbsp<b>Anisotropy (Part III)</b></span> par <span property="cc:attributionName">Magalhães et Cacciari (2026)</span> est licencié sous <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:10px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>
"""

# ╔═╡ 915727cc-9923-40ab-b0ed-bd178528ddbe
begin
	md"""
	# Considérations initiales
	
	Le matériel suivant a été développé comme support pour les cours de premier cycle du cours de Mécanique des Roches à Polytechnique Montréal. Ce notebook fait partie d'une série en trois parties sur le comportement anisotrope des roches, composée de :
	
	* Partie I - Déformation des roches anisotropes
	* Partie II - Critères de rupture anisotropes
	* Partie III - Tenseurs de structure

	Le contenu de ce notebook — _Partie III_ — descreve duas abordagens para análise de anisotropia a partir de abordagens tensoriais. O primeiro é baseado no fabric tensor proposto por Oda (1986) e o segundo consiste no microstructural tensor, de acordo com Pietruszczak (1999).
	
	Les matériaux sont élaborés par le groupe de recherche Geomechanics & Rock Mass Modeling et peuvent être utilisés à des fins éducatives à condition que la référence appropriée soit donnée. Pour des suggestions, corrections ou pour tout contact général, veuillez utiliser l'adresse courriel [gabrielgaldinodm@gmail.com](mailto:gabrielgaldinodm@gmail.com). 
	
	💡 **Le contenu est principalement destiné à la révision et à la pratique des concepts, et les étudiants sont encouragés à continuer à consulter la littérature.** 

	"""
end

# ╔═╡ e27bb47b-6851-4f24-b8cb-290c81d6c0e3
md"""
# 📚 Petite révision...

Des concepts introdutifs de mécanique et d’algèbre linéaire sont indispensables pour mieux comprendre les approches tensorielles. Répondez aux questions suivantes avant d'entrer dans les critères de rupture. 

> 1. Quelle est la limitation des critères de rupture qu'utilisent des directions spécifiques du tenseur des tensions (*e.g.* $\sigma_x, \tau_{xy}, \sigma_x \ge \sigma_{crit}$) ?
"""

# ╔═╡ b353f542-1b10-4310-be00-5240842c6f89
md"""> 2. Comment se représente un tenseur de second ordre dans l'espace ? Quelles sont les informations nécessaires pour le définir complètement ?"""

# ╔═╡ e819782b-d982-4d88-8180-02d0e3b63c16
md"""> 3. Comment la rotation d’un tenseur est-elle calculée selon un axe quelconque ?"""

# ╔═╡ 2f9b2337-a9e8-4651-b85f-05fe88535fb1
md"""
# 📖 Fabric tensor

Le tenseur de structure (ou fabric tensor) décrit la distribution géométrique des fractures dans un massif rocheux. Son avantage est de fournir une description continue de l’orientation et de l’intensité des fractures. La démarche suivante est basée sur Oda (1986, 1992, 1993) et Jing (2007).

##### 1) Distribution de taille et d’orientation des fractures :

Les fractures sont modélisées comme des plans de géométrie circulaire (rayon = $r$). Les dimensions des fractures peuvent être calculées à partir de la fonction de probabilité suivante :

$\int_{0}^{\infty} f(r) dr = 1 \quad\quad (Eq. 1)$

Afin de prendre en compte simultanément la taille et l’orientation des fractures, on définit ensuite une distribution conjointe taille-orientation, $E(\boldsymbol{n}, r)$, qui dépend du rayon et de l'orientation du vecteur normal sortant. Cela est normalisé sur l’espace des orientations et des tailles.

$\int_{0}^{\infty} \iint_{\Omega} E(\boldsymbol{n}, r) \, d\Omega \, dr = 1 \quad\quad (Eq. 2)$

où $d\Omega$ est l’incrément d’angle solide dans le système de coordonnées sphériques. 

##### 2) Densité des fractures:
Ici, nous allons considérer la probabilité qu’une scanline intercepte une fracture dans l’espace. De manière intuitive, on peut penser que cette probabilité dépend de l’aire de la fracture ($2\pi r^2$) ainsi que de la relation angulaire entre la scanline et le disque. Dans cette optique, nous considérons le produit scalaire entre le vecteur normal unitaire de la scanline ($\mathbf{s}$) et le vecteur normal unitaire du disque ($\mathbf{n}$). Toutefois, les fractures ne sont pas toutes identiques ; il est donc nécessaire de prendre en compte également les distributions de probabilité de taille et d’orientation.

$\eta = (\text{Aire projetée} \times \text{Densité de probabilité})$
$\eta = 2\pi r^2 |\vec{n} \cdot \vec{s}| E(\mathbf{n}, r) d\Omega dr \quad\quad (Eq. 3)$

La systématisation jusqu’ici repose sur la nature probabiliste de l’interception. Pour passer à une grandeur physique dans l’espace, nous introduisons la densité spatiale de fractures par $\rho = M_V/V$, où $M_V$ est le nombre total de fractures contenues dans le volume $V$.

##### 3) Effet directionnel:
Afin de prendre en compte l’effet directionnel dans le plan de fracture, l’auteur définit le vecteur $\mathbf{m} = 2r\mathbf{n_j}$, en fonction de la direction orthogonale ($n_j$) et du diamètre, ce qui conduit à :

$\eta \cdot \mathbf{m} = 4\pi r^3 n_i n_j E(\mathbf{n}, r) d\Omega dr \quad\quad (Eq. 4)$

##### 4) Tenseur de structure:
Ainsi, l'intégrale sur tout l'espace donne le tenseur de structure selon l'équation suivante :

$$F_{ij} = 4\pi\rho \int_{0}^{\infty} \int_{\Omega/2} r^3 n_i n_j E(\mathbf{n}, r) \, d\Omega \, dr\quad\quad (Eq. 5)$$

"""

# ╔═╡ 90cc5cd8-54f6-425c-98d2-e5ebb4f6a0f6
md"""
## 🖋️ Activité 1
La proposition de cette exercise est le calcul du tenseur de structure a parti d'un modèle de reseau de fractures (DFN). Quelles sont les étapes pour tracer un réseau de fractures discret et calculer le tenseur de structure ?

"""

# ╔═╡ 35f1d4ee-090b-4fe6-b049-5552e38fa2b6
md"""#### Chargement de la base de données

- **(x, y, z)** → Coordonnées du centre du plan de discontinuité  
- **fdip** → Pendage de la fracture  
- **fddir** → Direction du pendage  
- **fdia** → Diamètre du disque de la fracture  
- **fname** → Nom de l’ensemble de discontinuités (joint set)

"""

# ╔═╡ fe1def5e-f6ee-446d-8fb9-95e544566e82
begin
	chemin_notebook = @__DIR__
	chemin_base_donnees = joinpath(chemin_notebook, "..", "database", "fractures.txt")
	
	# Définir les noms ici
	noms_colonnes = ["id", "x", "y", "z", "fdip", "fdia", "fddir", "fname"]
	
	# La correction est ici : nous utilisons header pour passer les noms
	# et indiquons que les données commencent à la ligne 1
	df = CSV.read(chemin_base_donnees, DataFrame; 
                  header = noms_colonnes, 
                  skipto = 1,
                  delim = ' ')
end

# ╔═╡ 065f3d89-59cc-4f1d-bc4c-d43d772335f3
md"""#### Visualisation du modèle DFN"""

# ╔═╡ b196b879-1743-4bc8-bfcf-1f00d784ae77
begin
    # --- Fonctions auxiliaires ---
    function points_disque(rayon, npoints=30)
        θ = range(0, 2π, length=npoints)
        return [Point3f(rayon * cos(t), rayon * sin(t), 0) for t in θ]
    end

    function faire_rotations_points(points, dip, dipdir)
        θ, φ = deg2rad(dip), deg2rad(dipdir)
        n = [0.0, 0.0, 1.0]
        nouveau_n = [sin(θ)*sin(φ), sin(θ)*cos(φ), cos(θ)]
        axe = cross(n, nouveau_n)
        
        if norm(axe) < 1e-6
            R = I(3)
        else
            axe /= norm(axe)
            c = dot(n, nouveau_n)
            s = norm(cross(n, nouveau_n))
            K = [0.0 -axe[3] axe[2]; axe[3] 0.0 -axe[1]; -axe[2] axe[1] 0.0]
            R = I(3) + K + K^2 * ((1 - c)/(s^2))
        end
        return [Point3f(R * [p[1], p[2], p[3]]) for p in points]
    end

    # --- Traitement des données ---
    liste_disques = Vector{Point3f}[]
    liste_couleurs = RGBAf[]

    for i in 1:nrow(df)
        pts = points_disque(df.fdia[i] / 2)
        pts_rot = faire_rotations_points(pts, df.fdip[i], df.fddir[i])
        centre = Point3f(df.x[i], df.y[i], df.z[i])
        push!(liste_disques, [p + centre for p in pts_rot])
        
        couleur = df.fname[i] == "J1" ? RGBAf(1.0, 0.5, 0.0, 0.6) :
                  df.fname[i] == "J2" ? RGBAf(0.0, 0.8, 0.2, 0.6) :
                                         RGBAf(0.2, 0.4, 1.0, 0.6)
        push!(liste_couleurs, couleur)
    end

    # --- Visualisation dans un espace réduit ---
    
    # 1. Définir figure_padding à 0 ou très faible
    fig = Figure(size=(700, 450), figure_padding=0)
    
    # 2. Utiliser alignmode = Outside() pour que l'axe occupe tout l'espace
    ax = Axis3(fig[1, 1], 
               aspect=:data, 
               perspectiveness=0.5,
               title="", 
               titlevisible=false,
               xlabel="X", ylabel="Y", zlabel="Z",
               alignmode = Outside(0)) # Supprime le padding interne du layout

    poly!(ax, liste_disques, 
          color = liste_couleurs, 
          strokecolor = :black, 
          strokewidth = 0.5, 
          transparency = true)

    # 3. Forcer le cadrage sans marges supplémentaires
    autolimits!(ax)
    
    fig
end

# ╔═╡ d7ae3d7f-56c4-4f8b-9180-e41762297aa6
md"""#### Calculer les vecteurs normaux"""

# ╔═╡ 352dc524-b134-4b45-97a0-7793539b8a36
md"""
Pour effectuer la conversion des angles en composantes du vecteur normal, on utilise la formule suivante. Les mesures doivent être converties des degrés en radians avant son application.

$$\mathbf{n} = \begin{bmatrix} n_x \\ n_y \\ n_z \end{bmatrix} = \begin{bmatrix} \sin(\delta) \cdot \sin(\alpha) \\ \sin(\delta) \cdot \cos(\alpha) \\ -\cos(\delta) \end{bmatrix} \quad\quad (Eq.04)$$

où $\delta$ est l'angle de pendage ($fdip$), et $\alpha$ est la direction du pendage ($fddir$).
"""

# ╔═╡ 15b1edbf-fb22-491d-9860-f602006059ab
begin
	dip_rad = deg2rad.(df.fdip)
	ddir_rad = deg2rad.(df.fddir)
	
	# 3. Calcul des composantes du vecteur normal
	# Remarque : nous utilisons sin(dip) pour les composantes horizontales et cos(dip) pour la verticale.
	# Pour s'aligner avec le Nord (Y), nous utilisons sin(ddir) pour X et cos(ddir) pour Y.
	
	df.nx = sin.(dip_rad) .* sin.(ddir_rad)
	df.ny = sin.(dip_rad) .* cos.(ddir_rad)
	df.nz = -cos.(dip_rad) # Négatif car le vecteur normal pointe généralement vers le bas/dedans
	
	# Si vous préférez une seule colonne avec les vecteurs regroupés :
	df.vecteur_normal = collect.(zip(df.nx, df.ny, df.nz))
	
	df
end

# ╔═╡ b8d40fb0-8aa8-45fa-b16a-fcd3ddcfbcd1
md"""#### Calculer le tenseur de strcuture"""

# ╔═╡ 741ba561-8d47-4346-9828-2e8dd50e1463
md"""
La implementation numérique est basée dans l'equation suivante:

$$F_{ij} = \frac{2 \pi}{V} \sum_{k=1}^{N} (r^3 n_i n_j)_k$$
"""

# ╔═╡ b18ba78a-39ff-4e7b-a0b3-ba3abf097b09
md""" Resultat:"""

# ╔═╡ 59d6354c-7baa-4540-8a32-d9e0e32f8b45
begin
	# 1. volume (fixe)
	V = 10 * 10 * 10
	
	# 2. initialiser une matrice de zéros
	Fij = zeros(3, 3)
	
	# 3. contribution de chaque fracture
	for i in 1:nrow(df)
	    r_k = df.fdia[i] / 2
	    
	    # vecteur normal
	    n = [df.nx[i], df.ny[i], df.nz[i]]
	    
	    # produit externe
	    Fij += (r_k^3) * (n * n')
	end
	
	# 4. scalaire
	Fij *= 2π / V

	Fij
end

# ╔═╡ 3587624a-3d13-42e6-acdb-b01491ab3f86
md"""
## 🖋️ Activité 2
Calculer les vecteurs propres et les valeurs propres, et aprés, représente le tenseur de structure dans l'espace.
"""

# ╔═╡ eeb00d40-c0a5-42d6-b7e2-06f48c646bcc
md"""#### Valeurs propes
Par un tenseur de second ordre, on peut calculer les valeurs propres avec l'équation caractéristique. 

$det(A - \lambda I)=0$

Chaque racine du polynôme caractéristique va retourner une valeur propre. 
"""

# ╔═╡ 227a5f34-8580-4621-87dd-2992ff61c3a3
begin
    E = eigen(Symmetric(Fij))
    eigen_values = E.values
	nothing
end

# ╔═╡ 19e73cc3-68e8-471b-a7da-61dddfbd5d4c
md"""#### Vecteurs propres
Pour calculer les vecteurs propres on doit utiliser $\lambda$ et résoudre le système d'équations linéaires :

$(A- \lambda_i I) \mathbf{v} = 0$
"""

# ╔═╡ 47be2b68-8c7c-45a6-8dcd-683c492095b8
begin
	eigen_vectors = E.vectors
	nothing
end

# ╔═╡ c7c86f07-6de4-42ab-b20d-f0d936601655
md"""#### Visualisation du tenseur de structure"""

# ╔═╡ 2e0022ae-ba32-47d4-9fb7-e0e2df9abcdd
begin
	
	# Active le moteur WebGL pour une interactivité complète dans Pluto
	WGLMakie.activate!() 
	
	# 1. Données d'exemple (Valeurs propres et vecteurs propres)
	valeurs_propres = [4.0, 2.0, 1.0]
	vecteurs_propres = Float64.(I(3)) 
	
	# --- Construction de la géométrie ---
	u = range(0, 2π, length=40)
	v = range(0, π, length=20)
	
	# Matrice de transformation : Rotation * Échelle
	A = vecteurs_propres * diagm(sqrt.(valeurs_propres))
	
	# Génération des coordonnées X, Y, Z séparément pour éviter les erreurs de type
	xs_rot = [ (A * [cos(ui)*sin(vi), sin(ui)*sin(vi), cos(vi)])[1] for ui in u, vi in v ]
	ys_rot = [ (A * [cos(ui)*sin(vi), sin(ui)*sin(vi), cos(vi)])[2] for ui in u, vi in v ]
	zs_rot = [ (A * [cos(ui)*sin(vi), sin(ui)*sin(vi), cos(vi)])[3] for ui in u, vi in v ]

	# --- Visualisation ---
	fig2 = Figure(size=(700, 450), figure_padding=2)
	
	ax2 = Axis3(fig2[1, 1], 
	            aspect=:data, 
	            perspectiveness=0.5,
	            xlabel="X", ylabel="Y", zlabel="Z",
	            alignmode = Outside(0))

	# Surface fournit un volume cliquable. 
	# shading = true remplace l'ancien FastShading
	surface!(ax2, xs_rot, ys_rot, zs_rot, 
	         color=:lightblue, 
	         alpha=0.4, 
	         shading=true, 
	         transparency=true)
	
	# Wireframe au-dessus pour visualiser la grille
	wireframe!(ax2, xs_rot, ys_rot, zs_rot, 
	           color=:black, 
	           linewidth=0.5, 
	           transparency=true)

	autolimits!(ax2)
	
	fig2
end

# ╔═╡ c7f0eef4-cbb6-46fe-90a4-83381b76f01c
md"""
# 📖 Microstructural tensor
"""

# ╔═╡ ed241fbf-80c7-4959-a085-602fb8167637
md"""
Les matériaux cohésivo-frictionnels, dont le comportement est dicté par la cohésion et l'angle de friction, peuvent présenter une anisotropie due aux relations microstructurales. Cette anisotropie débute généralement par l'existence d'une orientation minérale et se développe ensuite par la propagation de microfissures. 

## Développement du microstructural tensor :

Considérons une sphère englobant un volume de matériau rocheux. Cette sphère de roche est traversée par une ligne théorique qui intercepte les pores et les minéraux de la matrice. Nous pouvons étiqueter chaque portion de cette ligne selon le contenu de l'intersection, ce qui donne une longueur totale correspondant aux pores et une longueur totale correspondant aux minéraux de la matrice.   

L'image suivante illustre l'expérience proposée :
"""

# ╔═╡ 95af722b-de22-4fa6-9631-6a9855f09dd1
html"""
<div style="display: flex; justify-content: center; background-color: white; padding: 10px; border-radius: 8px;">
	<img src="https://raw.githubusercontent.com/GaldinoMagalhaes/geomechanics_classes/e2b234f6af5263bd60d0971fe2b8d19d4f46a302/images/micro_tensor_en.png" width="800">
</div>
"""

# ╔═╡ fb0452c2-268d-4b6a-ac70-0689c4b1cb50
md"""
Maintenant, considérons que cette sphère est interceptée non pas par une seule, mais par plusieurs lignes test de longueur $\bar{L}$ et d'orientation $v_i$ qui remplissent complètement le volume de la roche. Nous pouvons décrire la fraction de $\bar{L}$ occupée par les vides par :

$L(v_i) = l(v_i)\bar{L}, \quad l(v_i) = \sum_{k} l_k(v_i) \quad (\text{Eq.05})$

où $l(v_i)$ représente le total des intersections de ces lignes avec les vides. La quantité homogénéisée de $L(v_i)$ dans le domaine $S$ (sphère) est donnée par :

$L_{av} = \frac{1}{4\pi} \oint_S L(v_i) g(v_i)dS \quad (\text{Eq.06})$

où $g(v_i)$ est une fonction scalaire décrivant la distribution spatiale uniforme des lignes test. En intégrant sur toute la surface, on peut encore affirmer que :

$\frac{1}{4\pi} \oint_S g(v_i)dS = 1 \quad (\text{Eq.07})$

Comme les lignes test sont utilisées pour compter les vides, on peut comprendre que les équations précédentes se rapportent à la porosité elle-même de la roche. 

$n = L_{av},\quad \bar{n}(v_i)\equiv L(v_i) \quad (\text{Eq.08})$
"""

# ╔═╡ 3963a903-6400-4c1d-9260-ca689b32bf4d
md"""
## 🖋️ Activité 3
Quelle est la différence entre la porosité globale et directionnelle ?
"""

# ╔═╡ 07f7f88e-608e-45f8-93dc-669c12c1009b
md"""
La relation qui décrit la porosité dans l'équation précédente peut être restructurée en fonction d’un tenseur ($\Omega_{ij}$) qui décrit l’anisotropie de la porosité via un ellipsoïde :

$n(v_i) = n_0(1 + \Omega_{ij}v_i v_j)\quad (\text{Eq.09})$

De manière intuitive, on peut considérer que $n_0$ est la composante isotrope, qui décrit la porosité moyenne dans n’importe quelle direction. Les composantes $v_i$ et $v_j$ indiquent la direction dans laquelle la mesure de porosité sera orientée.  

Pour des raisons de normalisation et d’intégration dans les modèles constitutifs, il est courant de convertir cette expression en Fabric Tensor ($A_{ij}$). Voici la déduction étape par étape :

**1) Définition du Fabric Tensor et du Delta de Kronecker :**

Le Fabric Tensor $A_{ij}$ est défini de sorte que la somme de ses composantes principales (sa trace) soit égale à l’unité ($A_{ii} = 1$). On utilise le **Delta de Kronecker** ($\delta_{ij}$) pour le relier au tenseur de déviations $\Omega_{ij}$ :

$\delta_{ij} = \begin{cases} 1, & \text{si } i = j \\ 0, & \text{si } i \neq j \end{cases}$

La relation de transformation est :

$A_{ij} = \frac{1}{3}(\delta_{ij} + \Omega_{ij})$

> **Remarque :** En trois dimensions, la trace du delta de Kronecker est $\delta_{ii} = 3$. Comme la trace de $\Omega_{ij}$ est nulle, cette division par 3 garantit que la trace de $A_{ij}$ est exactement 1 ($3/3 = 1$).

**2) Isolation du Tenseur d’Anisotropie :**

Pour substituer dans l’équation originale de la porosité, on isole $\Omega_{ij}$ :

$\Omega_{ij} = 3 A_{ij} - \delta_{ij}$

**3) Substitution et simplification tensorielle :**

En substituant cette expression dans l’Éq. 09 :

$n(v_i) = n_0(1 + (3 A_{ij} - \delta_{ij}) v_i v_j)$

$n(v_i) = n_0(1 + 3 A_{ij} v_i v_j - \delta_{ij} v_i v_j)$

Ici, on applique une propriété fondamentale du calcul tensoriel : le produit $\delta_{ij} v_i v_j$ représente le produit scalaire du vecteur unitaire par lui-même ($v \cdot v$). 

Comme le vecteur d’orientation $v$ est, par définition, unitaire :

$\delta_{ij} v_i v_j = v_i v_i = |v|^2 = 1$

**4) Résultat final :**

En substituant $\delta_{ij} v_i v_j = 1$ dans l’équation :

$n(v_i) = n_0(1 + 3 A_{ij} v_i v_j - 1)$

$n(v_i) = 3 n_0 A_{ij} v_i v_j\quad (\text{Eq.10})$

L’équation finale (Éq.10) permet de calculer la porosité dans n’importe quelle direction $v$ en connaissant uniquement la porosité moyenne $n_0$ et le Fabric Tensor $A_{ij}$ de la roche.

**5) Invariance à la rotation :**

Enfin, on sait que la porosité est invariante par rotation. Ainsi, en faisant tourner le système, notre ancien $v_p$ est relié au vecteur tourné $v^*_i$ par :

$v^*_i = \omega_{ip} v_p$

En remplaçant cela dans l’équation de porosité, on obtient :

$n(v^*) = 3n_0 A^*_{ij} v^*_i v^*_j$

Finalement, l’invariance est formalisée par :

$A_{pq} = \omega_{ip} \omega_{jq} A^*_{ij}$
"""

# ╔═╡ 326517af-10b0-45e8-b5d7-5ae823d9bb12
md"""
## 🖋️ Activité 3
Considere que uma dada rocha tem porosidade global $n_0 = 0.15$ e fabric tensor descrito por:

$A_{ij} = 
\begin{bmatrix}
0.5 & 0.1 & 0.0 \\
0.1 & 0.3 & 0.0 \\
0.0 & 0.0 & 0.2
\end{bmatrix}$

"""

# ╔═╡ 53184e30-ecbc-4df4-8f40-c43e617ccdaa
md"""
A) Calcule o tensor $\Omega_{ij}$ correspondente.
"""

# ╔═╡ 44683e23-0477-4617-8d1f-e1a957266dd6
# Input data
begin
	n₀ = 0.15
	A_ij = [0.5 0.1 0.0;
		 0.1 0.3 0.0;
		 0.0 0.0 0.2]
end

# ╔═╡ da823f5e-7947-4f7f-adea-a356bf2013c7
begin
	δ = Matrix{Float64}(I, 3, 3)  # I é a matriz identidade 3x3
	Ω = 3 .* A_ij .- δ
end

# ╔═╡ 72972205-f69f-4e50-a2bf-3fb3ff8b3d5c
md"""
B) Calcule a porosidade $n(\mathbf{v})$ para os seguintes vetores unitários:

$\mathbf{v}_1 = [1, 0, 0]$
$\mathbf{v}_2 = [0, 1, 0]$
$\mathbf{v}_3 = \frac{1}{\sqrt{2}} [1, 1, 0]$
"""

# ╔═╡ eecd032c-e1ad-4d31-a889-fa90f0235b34
begin
	# Vetores direcionais
	v₁ = [1.0, 0.0, 0.0]
	v₂ = [0.0, 1.0, 0.0]
	v₃ = [1/sqrt(2), 1/sqrt(2), 0.0]
end

# ╔═╡ a0836c3f-d262-46e6-822c-2ffcbd3ebf0c
begin
	# Equação de porosidade direcional
	function porosidade_dir(A, n₀, v)
    	return 3 * n₀ * dot(v, A * v)  # n(v) = 3 n₀ v^T A v
	end
end

# ╔═╡ 81103522-594c-4b8e-b025-95a1b12bb2ff
begin
	# Calcular porosidade para cada vetor
	n_v1 = porosidade_dir(A_ij, n₀, v₁)
	n_v2 = porosidade_dir(A_ij, n₀, v₂)
	n_v3 = porosidade_dir(A_ij, n₀, v₃)
	nothing
end

# ╔═╡ 43b0f1b5-bc42-43f8-ae43-625bb39d4b98
md"""
| Autovetor | Porosidade Direcional ($n_{v}$) |
|:---:|:---:|
| $v_1$ | $(round(n_v1, digits=4)) |
| $v_2$ | $(round(n_v2, digits=4)) |
| $v_3$ | $(round(n_v3, digits=4)) |

"""

# ╔═╡ 88232c7e-70c4-4daf-8fed-ea2aa12576dc
md"""
Explique qual direção apresenta maior porosidade e o que isso indica sobre a anisotropia do material.
"""

# ╔═╡ ddddaf7b-2a0b-4ed1-a572-e791646e8af1
md"""
# Littérature recommandée

[1] L. Jing and O. Stephansson, Fundamentals of Discrete Element Methods for Rock Engineering: Theory and Applications, Elsevier Science, pp. 47–109, 2007.

[2] M. Oda, “Fabric tensor for discontinuous geological materials,” Soils and Foundations, vol. 22, no. 4, pp. 96–108, 1982.

[3] M. Oda, “Permeability tensor for discontinuous rock masses,” Geotechnique, vol. 35, no. 4, pp. 483–495, 1985.

[4] S. Pietruszczak, “On inelastic behaviour of anisotropic frictional materials,” Mechanics of Cohesive-Frictional Materials, vol. 4, pp. 281–293, 1999.

[5] S. Pietruszczak, Fundamentals of Plasticity in Geomechanics, 1st ed., John Wiley & Sons, 2010.

[6] S. Pietruszczak and Z. Mroz, “Formulation of anisotropic failure criteria incorporating a microstructure tensor,” Computers and Geotechnics, vol. 26, no. 2, pp. 105–112, 2000.

"""

# ╔═╡ 94026925-bb5f-4ec4-8eef-0b72e3b0b133
md"""##### Fonctions de support"""

# ╔═╡ ea0a1a2c-19c2-4fdb-8f1e-8cf099996a7f
begin
	# Auxiliare functions
	# Hint function
	hint(text) = (Markdown.MD(Markdown.Admonition("hint", 
		"Astuce", [text])))
end

# ╔═╡ 63d60cc8-165a-4522-9139-087fb99d27f8
hint(md"""
Les critères de rupture doivent être invariants par rapport à une transformation orthogonale. Les critères qui considèrent une direction spécifique (comme les composantes spécifiques du tenseur de contrainte) peuvent varier selon la direction. Les modèles à formulation tensorielle présentés ici sont invariants.
""")

# ╔═╡ c11da0d3-ca06-4ecf-8044-5ea57f36ba07
hint(md"""
Comme un ellipsoïde. Il est nécessaire de calculer les vecteurs propres et les valeurs propres.
""")

# ╔═╡ f271ca5b-e211-4ef4-bf73-97011e669abf
hint(md"""
Étape 1 : Donner un vecteur unitaire $u = (u_x, u_y, u_z)$ et un angle $\theta$ qui décrivent l’axe de rotation.

Étape 2 : Calculer la matrice de rotation à l’aide de la formule de Rodrigues.

$$\mathbf{R} = \cos\theta \mathbf{I} + \sin\theta [\hat{u}]_\times + (1 - \cos\theta) (\hat{u} \otimes \hat{u})$$

Étape 3 : Utiliser la matrice de rotation pour calculer
$$\mathbf{T}' = \mathbf{R} \mathbf{T} \mathbf{R}^\mathsf{T}$$.
""")

# ╔═╡ 8ab84c6c-d67d-4758-b6d4-3009bed804ea
hint(md"""
1. **Chargement** de la base de données  
2. **Visualisation** du modèle DFN  
3. **Calculer** les vecteurs normaux  
4. **Calculer** le tenseur de structure
""")

# ╔═╡ f5e9f6fc-071f-4c0f-82a9-479882176113
hint(md"""
Utilize a relação: $\Omega_{ij} = 3A_{ij} - \delta_{ij}$
""")

# ╔═╡ 3d06fa33-5130-4d7b-8384-34bc5e2a6836
hint(md""" Utilize a relação: 
$n(v_i) = 3 n_0 A_{ij} v_i v_j$
""")

# ╔═╡ c3d0fbc7-64e2-41e5-9c9c-05b198911f43
begin
	ans(text) = (Markdown.MD(Markdown.Admonition("ans", 
		"Réponse", [text])))
end

# ╔═╡ bad96012-f41a-4afb-8eda-4f1569175f69
ans(md"""
| Index|Value |
|--------|-------|
| λ₁ | $(round(eigen_values[1], sigdigits=4)) |
| λ₂ | $(round(eigen_values[2], sigdigits=4)) |
| λ₃ | $(round(eigen_values[3], sigdigits=4)) |
""")

# ╔═╡ b2ef6d19-538a-4347-b59c-8125aa9c5105
ans(md"""
|      | v₁ | v₂ | v₃ |
|------|----|----|----|
| x    | $(round(eigen_vectors[1,1], sigdigits=4)) | $(round(eigen_vectors[1,2], sigdigits=4)) | $(round(eigen_vectors[1,3], sigdigits=4)) |
| y    | $(round(eigen_vectors[2,1], sigdigits=4)) | $(round(eigen_vectors[2,2], sigdigits=4)) | $(round(eigen_vectors[2,3], sigdigits=4)) |
| z    | $(round(eigen_vectors[3,1], sigdigits=4)) | $(round(eigen_vectors[3,2], sigdigits=4)) | $(round(eigen_vectors[3,3], sigdigits=4)) |
""")

# ╔═╡ 5b57d8f2-fceb-407b-b499-67f93754f308
ans(md"""
La porosité globale se réfère à une valeur qui exprime le rapport entre la matrice et les vides dans une roche de manière générale, une valeur scalaire moyenne, tandis que la porosité directionnelle, comme son nom l’indique, capture cette relation selon une orientation préférentielle. Il est intéressant de noter que pour une roche isotrope, toute ligne de scan suffit à obtenir une valeur représentative de porosité. Pour les roches anisotropes, en revanche, il existe une dépendance de la porosité vis-à-vis de l’orientation structurelle.
""")

# ╔═╡ eec763be-c7aa-41f2-91e6-c49ffabb7dad
# Célula de Exibição Markdown
ans(md"""
| $\Omega$ | c₁ | c₂ | c₃ |
|:---:|:---:|:---:|:---:|
| **x** | $(round(Ω[1,1], digits=2)) | $(round(Ω[1,2], digits=2)) | $(round(Ω[1,3], digits=2)) |
| **y** | $(round(Ω[2,1], digits=2)) | $(round(Ω[2,2], digits=2)) | $(round(Ω[2,3], digits=2)) |
| **z** | $(round(Ω[3,1], digits=2)) | $(round(Ω[3,2], digits=2)) | $(round(Ω[3,3], digits=2)) |
""")

# ╔═╡ 1d675e29-386f-42ed-b3e9-b92ad27ed6bc
ans(md"""

O material apresenta maiores valores de porosidade direcional ao longo dos eixos $v_1$ e $v_3$. Isso indica que a rocha apresenta anisotropia e que ao longo das direções selecionadas existem mais poros, implicando naturalmente na forma com fluidos em geral atravessam essa rocha. 

""")

# ╔═╡ a5c3fc7c-98d8-487f-a37c-51bfb8cf76cd
	# Summary
	PlutoUI.TableOfContents(aside=true, title="Sommaire",
							indent=true, depth=2)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
WGLMakie = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"

[compat]
CSV = "~0.10.15"
Colors = "~0.13.1"
DataFrames = "~1.8.1"
PlutoUI = "~0.7.78"
WGLMakie = "~0.13.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.8"
manifest_format = "2.0"
project_hash = "8f28157f8996d65860b6ef8c37d975fff7034534"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "SIMD", "TranscodingStreams"]
git-tree-sha1 = "a8f503e8e1a5f583fbef15a8440c8c7e32185df2"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BaseDirs]]
git-tree-sha1 = "bca794632b8a9bbe159d56bf9e31c422671b35e0"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.3.2"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bonito]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "JSON", "LinearAlgebra", "Markdown", "MbedTLS", "MsgPack", "Observables", "OrderedCollections", "Random", "RelocatableFolders", "SHA", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "bb43f72801f703ad3c66833bd02b8f54c7328238"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
version = "4.2.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

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
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "76dab592fa553e378f9dd8adea16fe2591aa3daa"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.1.6"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "c55f5a9fd67bdbc8e089b5a3111fe4292986a8e8"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.6.6"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "7bebc8aad6ee6217c78c5ddcf7ed289d65d0263e"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.6"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

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

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6522cfb3b8fe97bec632252263057996cbd3de20"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.18.0"
weakdeps = ["HTTP"]

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

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

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "IterTools", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "1f5a80f4ed9f5a4aada88fc2db456e637676414b"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.10"

    [deps.GeometryBasics.extensions]
    GeometryBasicsGeoInterfaceExt = "GeoInterface"

    [deps.GeometryBasics.weakdeps]
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6570366d757b50fabae9f4315ad74d2e40c0560a"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "5.2.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "6b4d2dc81736fe3980ff0e8879a9fc7c33c44ddf"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "93d5c27c8de51687a2c70ec0716e6e76f298416f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

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

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "696144904b76e1ca433b886b4e7edd067d76cbf7"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.9"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "02b61501dbe6da3b927cc25dacd7ce32390ee970"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.2"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

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

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

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

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "d1b974f376c24dad02c873e951c5cd4e351cd7c2"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.24.8"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "7eb8cdaa6f0e8081616367c10b31b9d9b34bb02a"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "567515ca155d0020a45b05175449b499c63e7015"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.29+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "df9b7c88c2e7a2e77146223c526bf9e236d5f450"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.4.4+0"

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

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "cf181f0b1e6a18dfeb0ee8acc4a9d1672499626c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.4"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

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

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "c5a07210bd060d6a8491b0ccdee2fa0235fc00bf"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.1.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "e24dc23107d426a096d3eae6c165b921e74c18e4"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "818554664a2e01fc3784becb2eb3a82326a604b6"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "0494aed9501e7fb65daba895fb7fd57cc38bc743"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.5"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eee1b9ad8b29ef0d936e3ec9838c7ec089620308"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.16"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a3c1536470bf8c5e02096ad4853606d7c8f62721"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.2"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

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

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "98b9352a24cb6a2066f9ababcc6802de9aed8ad8"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.6"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

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

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
    Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.WGLMakie]]
deps = ["Bonito", "Colors", "FileIO", "FreeTypeAbstraction", "GeometryBasics", "Hyperscript", "LinearAlgebra", "Makie", "Observables", "PNGFiles", "PrecompileTools", "RelocatableFolders", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "32801246eb6c7afb0e1e49509b3ffebecb538657"
uuid = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
version = "0.13.8"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

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

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

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

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "c1733e347283df07689d71d61e14be986e49e47a"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.5+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "4e4282c4d846e11dce56d74fa8040130b7a95cb3"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.0+0"

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
"""

# ╔═╡ Cell order:
# ╠═b5719bd8-4f04-4d42-a081-7504d946129f
# ╟─a27f367e-5290-4a78-b244-14ba9b9b7c2d
# ╟─c44c680a-1c0a-4fbe-af3a-af512ce272f4
# ╟─915727cc-9923-40ab-b0ed-bd178528ddbe
# ╟─e27bb47b-6851-4f24-b8cb-290c81d6c0e3
# ╟─63d60cc8-165a-4522-9139-087fb99d27f8
# ╟─b353f542-1b10-4310-be00-5240842c6f89
# ╟─c11da0d3-ca06-4ecf-8044-5ea57f36ba07
# ╟─e819782b-d982-4d88-8180-02d0e3b63c16
# ╟─f271ca5b-e211-4ef4-bf73-97011e669abf
# ╟─2f9b2337-a9e8-4651-b85f-05fe88535fb1
# ╟─90cc5cd8-54f6-425c-98d2-e5ebb4f6a0f6
# ╟─8ab84c6c-d67d-4758-b6d4-3009bed804ea
# ╟─35f1d4ee-090b-4fe6-b049-5552e38fa2b6
# ╠═fe1def5e-f6ee-446d-8fb9-95e544566e82
# ╟─065f3d89-59cc-4f1d-bc4c-d43d772335f3
# ╠═b196b879-1743-4bc8-bfcf-1f00d784ae77
# ╟─d7ae3d7f-56c4-4f8b-9180-e41762297aa6
# ╟─352dc524-b134-4b45-97a0-7793539b8a36
# ╠═15b1edbf-fb22-491d-9860-f602006059ab
# ╟─b8d40fb0-8aa8-45fa-b16a-fcd3ddcfbcd1
# ╟─741ba561-8d47-4346-9828-2e8dd50e1463
# ╟─b18ba78a-39ff-4e7b-a0b3-ba3abf097b09
# ╠═59d6354c-7baa-4540-8a32-d9e0e32f8b45
# ╟─3587624a-3d13-42e6-acdb-b01491ab3f86
# ╟─eeb00d40-c0a5-42d6-b7e2-06f48c646bcc
# ╠═227a5f34-8580-4621-87dd-2992ff61c3a3
# ╟─bad96012-f41a-4afb-8eda-4f1569175f69
# ╟─19e73cc3-68e8-471b-a7da-61dddfbd5d4c
# ╠═47be2b68-8c7c-45a6-8dcd-683c492095b8
# ╟─b2ef6d19-538a-4347-b59c-8125aa9c5105
# ╟─c7c86f07-6de4-42ab-b20d-f0d936601655
# ╟─2e0022ae-ba32-47d4-9fb7-e0e2df9abcdd
# ╟─c7f0eef4-cbb6-46fe-90a4-83381b76f01c
# ╟─ed241fbf-80c7-4959-a085-602fb8167637
# ╟─95af722b-de22-4fa6-9631-6a9855f09dd1
# ╟─fb0452c2-268d-4b6a-ac70-0689c4b1cb50
# ╟─3963a903-6400-4c1d-9260-ca689b32bf4d
# ╟─5b57d8f2-fceb-407b-b499-67f93754f308
# ╟─07f7f88e-608e-45f8-93dc-669c12c1009b
# ╟─326517af-10b0-45e8-b5d7-5ae823d9bb12
# ╟─53184e30-ecbc-4df4-8f40-c43e617ccdaa
# ╟─f5e9f6fc-071f-4c0f-82a9-479882176113
# ╠═44683e23-0477-4617-8d1f-e1a957266dd6
# ╠═da823f5e-7947-4f7f-adea-a356bf2013c7
# ╟─eec763be-c7aa-41f2-91e6-c49ffabb7dad
# ╟─72972205-f69f-4e50-a2bf-3fb3ff8b3d5c
# ╟─3d06fa33-5130-4d7b-8384-34bc5e2a6836
# ╠═eecd032c-e1ad-4d31-a889-fa90f0235b34
# ╠═a0836c3f-d262-46e6-822c-2ffcbd3ebf0c
# ╠═81103522-594c-4b8e-b025-95a1b12bb2ff
# ╟─43b0f1b5-bc42-43f8-ae43-625bb39d4b98
# ╟─88232c7e-70c4-4daf-8fed-ea2aa12576dc
# ╟─1d675e29-386f-42ed-b3e9-b92ad27ed6bc
# ╟─ddddaf7b-2a0b-4ed1-a572-e791646e8af1
# ╟─94026925-bb5f-4ec4-8eef-0b72e3b0b133
# ╟─ea0a1a2c-19c2-4fdb-8f1e-8cf099996a7f
# ╟─c3d0fbc7-64e2-41e5-9c9c-05b198911f43
# ╠═a5c3fc7c-98d8-487f-a37c-51bfb8cf76cd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
