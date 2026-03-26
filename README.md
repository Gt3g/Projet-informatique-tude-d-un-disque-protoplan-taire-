# Projet-informatique-tude-d-un-disque-protoplan-taire-
Projet info Ludivine Baptiste Nils


DESCRIPTION PROGRAMMES:

-particules.h/.cpp : classe particule avec position (x,y,z), vitesse(vx,vy,vz), masse et charge(pas utilisé,utilisé comme moment_cin mais finallement non)

-nuage.h/.cpp : classe nuage (vecteur de particule), initialisation avec(nbr de particule,rayon_particule,masse_part,rayon nuage,masse_étoile,rayon étoile,vitesse de rotation de nuage)
                les part initialisées dans un volume de rayon ray_nuage exclue du cylindre centré sur l'axe z et de rayon_étoile, vitesse initial = distance à l'axe*vitesse de rotation.

-nuage_solver.h/.cpp : évolution du système en fonction du temps
                       -position et vitesse de chaque part dans une unique liste pour passer dans rk4
                       -système avec les différentes forces (gravit avec l'étoile, frottement selon z (sinon le système se stabilise pas)
                       -résolution rk4 et écriture dans le fichier evolution_nuage

-main_nuage_etoile.cpp : déclaration des constantes et des différents parramètres
                         calcul du temps de simulation en fonction du temps de chute libre
                         définition des paramètre temporelle ( t0, tEnd prop à tff, pas en fraction de tff)
                         initialisation du nuage
                         simulation rk4







# Simulation de formation d'un disque protoplanétaire

Simulation N-corps d'un nuage moléculaire sphérique homogène en rotation
s'effondrant sous sa propre gravité pour former un disque protoplanétaire.

## Physique simulée

1. **Nuage initial** : sphère uniforme en rotation solide autour de Z,
   avec un paramètre de rotation β = E_rot / |E_grav| (valeur typique : 0.04–0.07).
2. **Effondrement gravitationnel** : auto-gravité étoile ↔ grains (+ grain ↔ grain optionnel).
3. **Formation du disque** : conservation du moment cinétique → aplatissement.
4. **Dissipations différées** (activées à t > T_DISS_START × t_ff) :
   - amortissement vertical (modèle disque mince irradié)
   - circularisation des orbites (amortissement de la vitesse radiale)

L'intégration temporelle est réalisée par un schéma **Runge-Kutta d'ordre 4**.

## Structure des fichiers

```
.
├── main_nuage_etoile.cpp   # Point d'entrée, initialisation, lancement de la simulation
├── params.h                # Tous les paramètres physiques et numériques (modifier ici)
├── particule.h / .cpp      # Classe Particule (position, vitesse, masse)
├── nuage.h / .cpp          # Classe Nuage (ensemble de particules + initialisation)
├── nuage_solver.h / .cpp   # Solveur RK4, calcul des dérivées, diagnostics, I/O
├── forces.h / .cpp         # Gravité softened, amortissement vertical, circularisation
├── Makefile
└── README.md
```

## Compilation

### Prérequis

- Compilateur C++ compatible C++14 ou supérieur (`g++` ≥ 5 ou `clang++`)
- `make`

### Compiler et exécuter

```bash
make          # compile le projet
make run      # compile puis lance la simulation
make clean    # supprime les fichiers objets et l'exécutable
```

### Compiler en mode debug

```bash
make debug
```

## Paramètres

Tous les paramètres sont centralisés dans **`params.h`**. Recompiler après toute modification (`make`).

| Paramètre         | Valeur par défaut   | Description                                      |
|-------------------|---------------------|--------------------------------------------------|
| `N_PARTICULES`    | 500                 | Nombre de super-grains                           |
| `MASSE_ETOILE`    | 1 M☉                | Masse de l'étoile centrale                       |
| `FRACTION_DISK`   | 0.10                | Fraction de masse stellaire pour le disque       |
| `RAYON_NUAGE`     | 120 UA              | Rayon initial du nuage                           |
| `BETA_ROT`        | 0.5                 | Rapport E_rot / \|E_grav\|                       |
| `DUREE_EN_TFF`    | 20.0                | Durée totale en unités de t_ff                   |
| `PAS_EN_TFF`      | 300.0               | Nombre de pas RK4 par t_ff                       |
| `ENABLE_GRAIN_GRAIN` | false            | Active l'auto-gravité grain ↔ grain (O(N²))      |
| `T_DISS_START`    | 1.0                 | Activation des dissipations (en unités de t_ff)  |
| `EPS_FACTOR`      | 0.15                | Facteur de softening gravitationnel              |
| `F_CIRC`          | 1.1                 | Taux de circularisation (sans dimension)         |
| `F_DISS_DISK`     | 2.0                 | Coefficient d'amortissement vertical             |
| `PAS_ECRITURE`    | 5                   | Fréquence d'écriture des positions (en pas)      |
| `PAS_DIAG`        | 1                   | Fréquence d'écriture des diagnostics (en pas)    |

## Fichiers de sortie

| Fichier                | Contenu                                                   |
|------------------------|-----------------------------------------------------------|
| `evolution_nuage.fich` | Trajectoires : `t  i  x  y  z  vx  vy  vz` (SI)         |
| `diagnostics.fich`     | `t  E_cin  E_grav  E_tot  L_z  dLz_rel` (SI)             |

### Validation

Sans dissipation (`T_DISS_START` très grand), `E_tot` et `L_z` doivent être conservés.
Vérifier `dLz_rel` dans `diagnostics.fich` : doit rester petit (< ~1%).

## Visualisation (Python)

Exemple minimal pour visualiser les positions dans le plan :

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("evolution_nuage.fich", comments="#")
# Dernier pas de temps
t_last = data[:, 0].max()
snap = data[data[:, 0] == t_last]
plt.scatter(snap[:, 2] / 1.496e11, snap[:, 3] / 1.496e11, s=1)
plt.xlabel("x [UA]")
plt.ylabel("y [UA]")
plt.axis("equal")
plt.title(f"t = {t_last / 3.156e7:.1f} ans")
plt.show()
```

## Références

- Goodman et al. (1993), ApJ 406 — rotation des cœurs denses
- Papaloizou & Larwood (2000) — circularisation des orbites
- Tanaka & Ward (2004) — migration planétaire
- Fromang & Nelson (2006), A&A 457 — amortissement dans les disques
- Dehnen (2001) — softening gravitationnel
