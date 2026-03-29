# Simulation de formation d'un disque protoplanétaire

Simulation N-corps d'un nuage moléculaire sphérique homogène en rotation s'effondrant sous
sa propre gravité pour former un disque protoplanétaire, intégrée par Runge-Kutta d'ordre 4.

---

## Table des matières

1. [Contexte physique](#1-contexte-physique)
2. [Modèle mathématique](#2-modèle-mathématique)
   - 2.1 [État initial](#21-état-initial)
   - 2.2 [Équations du mouvement](#22-équations-du-mouvement)
   - 2.3 [Forces modélisées](#23-forces-modélisées)
3. [Méthode numérique](#3-méthode-numérique)
   - 3.1 [Intégrateur RK4](#31-intégrateur-rk4)
   - 3.2 [Représentation de l'état](#32-représentation-de-létat)
   - 3.3 [Méthode Particle-Mesh pour la pression](#33-méthode-particle-mesh-pour-la-pression)
4. [Structure du projet](#4-structure-du-projet)
   - 4.1 [Arborescence](#41-arborescence)
   - 4.2 [Description des fichiers](#42-description-des-fichiers)
   - 4.3 [Diagramme de dépendances](#43-diagramme-de-dépendances)
5. [Compilation et exécution](#5-compilation-et-exécution)
   - 5.1 [Prérequis](#51-prérequis)
   - 5.2 [Commandes make](#52-commandes-make)
6. [Paramètres de simulation](#6-paramètres-de-simulation)
   - 6.1 [Paramètres physiques](#61-paramètres-physiques)
   - 6.2 [Paramètres numériques](#62-paramètres-numériques)
   - 6.3 [Dissipations](#63-dissipations)
   - 6.4 [Pression (méthode PM)](#64-pression-méthode-pm)
7. [Fichiers de sortie](#7-fichiers-de-sortie)
8. [Validation](#8-validation)
9. [Limitations et pistes d'amélioration](#10-limitations-et-pistes-damélioration)

---

## 1. Contexte physique

Les étoiles se forment à partir de l'effondrement gravitationnel de nuages moléculaires
denses (typiquement 50–200 UA de rayon, quelques dizaines de températures kelvin). Si le
nuage possède un moment cinétique non nul — ce qui est quasi-universel dans les
observations — la conservation de ce moment cinétique lors de l'effondrement empêche la
matière de tomber directement sur l'étoile centrale. Un disque protoplanétaire se forme
alors autour de la protoétoile.

Ce code simule ce processus en représentant le gaz et les poussières du nuage par un
ensemble de **super-particules** (particules macro-physiques, chacune représentant une
fraction de la masse totale du disque), et en intégrant numériquement leurs équations du
mouvement sous l'effet :

- de la gravité de l'étoile centrale,
- de la gravité mutuelle grain↔grain (optionnel),
- d'une force de pression du gaz (méthode Particle-Mesh),
- de forces de dissipation différées (amortissement vertical, circularisation des orbites).

---

## 2. Modèle mathématique

### 2.1 État initial

Le nuage est initialisé comme une sphère homogène de rayon `RAYON_NUAGE` en rotation
solide autour de l'axe Z, avec une vitesse angulaire `ω₀` définie par le paramètre sans
dimension β (rapport énergie de rotation / énergie gravitationnelle) :

```
ω₀ = sqrt( β · 3 G M / R³ )
```

Les positions des grains sont tirées aléatoirement et uniformément dans la sphère
(méthode de rejet), en excluant la zone centrale occupée par l'étoile (`R_cyl < RAYON_ETOILE`).
Les vitesses initiales combinent la rotation solide et une agitation thermique gaussienne :

```
vx = -ω₀ · y  +  v_therm · N(0,1)
vy = +ω₀ · x  +  v_therm · N(0,1)
vz =             v_therm · N(0,1)
```

La graine du générateur aléatoire est fixée (42) pour garantir la reproductibilité.

### 2.2 Équations du mouvement

Pour chaque grain i (i ≠ étoile), le système d'EDO intégré est :

```
d/dt [xᵢ, yᵢ, zᵢ] = [vxᵢ, vyᵢ, vzᵢ]

d/dt [vxᵢ, vyᵢ, vzᵢ] = a_grav,i  +  a_pression,i  +  a_diss,i
```

L'étoile centrale (indice 0) est maintenue fixe à l'origine (masse fixe, pas d'équation
du mouvement calculée pour elle).

### 2.3 Forces modélisées

#### Gravité softened (Newton + adoucissement de Plummer)

Évite la divergence numérique en r → 0 :

```
a_grav(i ← j) = G · mⱼ · Δr / (|Δr|² + ε²)^(3/2)
```

Le paramètre de softening `ε` est calculé automatiquement :

```
ε = EPS_FACTOR · R_nuage / N^(1/3)
```

Ce choix assure que ε est de l'ordre de la distance inter-particulaire moyenne.

La gravité grain↔grain (boucle O(N²)) est activable via `ENABLE_GRAIN_GRAIN`. Elle est
physiquement négligeable lorsque la masse du disque représente seulement quelques pourcents de la
masse stellaire.

#### Pression du gaz (Particle-Mesh)

Voir section 3.3 pour le détail de l'algorithme.

Loi d'état polytropique : `P = K · ρ^γ`

- γ = 1 (isotherme) : `P = cs² · ρ`, modèle par défaut
- γ = 5/3 (adiabatique) : disque chaud

#### Amortissement vertical (dissipation dans le disque)

Modèle disque mince irradié passif : la température suit un profil en loi de puissance :

```
T(R) = T₀ · sqrt(R₀ / R)
```

La hauteur de pression locale du disque vaut :

```
H(R) = cs(R) / Ω_K(R)   ~  R^(5/4)
```

Si `|zᵢ| < H(Rᵢ)`, une force d'amortissement dissipe la vitesse verticale :

```
a_z,diss = -F_DISS_DISK · Ω_K · vz
```

#### Circularisation des orbites

Amortit la composante radiale de la vitesse en coordonnées cylindriques, sans affecter
la composante azimutale (moment cinétique L_z conservé) :

```
F_e = -(f_circ · Ω_K) · v_r · r̂_cyl
```

Le temps de circularisation typique est `t_e ~ (1/f_circ) / (2π)` périodes orbitales.

Les dissipations sont activées seulement à partir de `t = T_DISS_START · t_ff` pour ne
pas biaiser la phase d'effondrement sphérique.

---

## 3. Méthode numérique

### 3.1 Intégrateur RK4

L'intégration temporelle utilise un schéma de **Runge-Kutta d'ordre 4** classique, avec
un pas de temps fixe `h` :

```
k1 = f(s, t)
k2 = f(s + h/2 · k1, t + h/2)
k3 = f(s + h/2 · k2, t + h/2)
k4 = f(s + h · k3, t + h)

s(t+h) = s(t) + (h/6) · (k1 + 2·k2 + 2·k3 + k4)
```

L'erreur locale est en O(h⁵), l'erreur globale en O(h⁴). Le pas est exprimé en fraction
du temps de chute libre `t_ff` via `PAS_EN_TFF`.

### 3.2 Représentation de l'état

L'état complet du système de N particules est aplati dans un vecteur `std::vector<double>`
de taille 6N :

```
s = [ x₀, y₀, z₀, vx₀, vy₀, vz₀,
      x₁, y₁, z₁, vx₁, vy₁, vz₁,  ...  ]
```

Les masses sont stockées séparément dans le `Nuage` (invariantes au cours du temps).
Des fonctions `nuage_to_state` / `state_to_nuage` assurent la conversion entre les deux
représentations.

### 3.3 Méthode Particle-Mesh pour la pression

La pression du gaz est calculée en O(N + Ng³) au lieu de O(N²) grâce à une grille
cubique de résolution `Ng` (paramètre `N_GRID`), sur un domaine cylindrique
`|R_cyl| ≤ R_DISK_PM`, `|z| ≤ H_DISK_PM`.

**Algorithme en 4 étapes :**

1. **Dépôt Cloud-In-Cell (CIC)** : chaque particule distribue sa masse sur les 8 cellules
   voisines par interpolation trilinéaire (poids proportionnel au volume de recouvrement).
   → grille de densité `ρ[ix][iy][iz]`

2. **Loi d'état** : `P[ix][iy][iz] = K · ρ^γ` sur chaque cellule.

3. **Gradient de pression** : différences finies centrées à l'intérieur,
   condition de Neumann (gradient nul) aux bords.
   ```
   (∂P/∂x)[i] = (P[i+1] - P[i-1]) / (2 dx)
   ```

4. **Interpolation inverse CIC** : relecture du gradient aux positions des particules
   avec les mêmes poids CIC.
   ```
   a_press,i = -(1/ρᵢ) · ∇P(rᵢ)
   ```

Les particules hors du cylindre ne déposent pas de masse et ne reçoivent pas de force de
pression (la pression du gaz est négligeable en dehors du disque).

---

## 4. Structure du projet

### 4.1 Arborescence

```
.
├── main_nuage_etoile.cpp   # Point d'entrée : initialisation, affichage, lancement
├── params.h                # Tous les paramètres physiques et numériques (MODIFIER ICI)
├── particule.h / .cpp      # Classe Particule (position, vitesse, masse)
├── nuage.h / .cpp          # Classe Nuage (collection de particules + initialisation)
├── nuage_solver.h / .cpp   # Solveur RK4, dérivée du système, diagnostics, I/O
├── forces.h / .cpp         # Gravité softened, amortissement vertical, circ., PM
├── Makefile                # Compilation, debug, valgrind, plot
└── README.md               # Ce fichier
```

### 4.2 Description des fichiers

| Fichier | Rôle |
|---|---|
| `params.h` | Centralise **toutes** les constantes physiques et numériques sous forme de `constexpr`. C'est le seul fichier à modifier pour changer les paramètres. |
| `particule.h/cpp` | Classe `Particule` : stocke position, vitesse, masse et moment cinétique scalaire. Fournit `afficher()`, `dist_origine()`, `dist_axe_carre()`. |
| `nuage.h/cpp` | Classe `Nuage` : vecteur de `Particule`. Méthode `init_nuage_homogene_etoile()` qui place l'étoile en [0] et génère les grains par rejet dans la sphère. |
| `nuage_solver.h/cpp` | Fonction principale `rk4_nuage()`. Contient aussi : `derivee()` (calcul des accélérations), `rk4_step()` (un pas RK4), `compute_diagnostics()` (E_cin, E_grav, L_z), `write_step()` / `write_diag()` (écriture fichiers). |
| `forces.h/cpp` | `gravite_softened()`, `amortissement_vertical()`, `circularisation()`, `pression_grille()` (méthode PM). |
| `main_nuage_etoile.cpp` | Calcule ω₀ et t_ff, remplit `SimParams`, crée le `Nuage`, appelle `rk4_nuage()`. |

### 4.3 Diagramme de dépendances

```
main_nuage_etoile.cpp
    ├── params.h
    ├── nuage.h ──────── particule.h
    └── nuage_solver.h
            ├── nuage.h
            ├── params.h
            └── forces.h ─── nuage.h
                              nuage_solver.h
                              params.h
```

---

## 5. Compilation et exécution

### 5.1 Prérequis

- Compilateur C++ compatible C++14 : `g++` ≥ 5 ou `clang++` ≥ 3.4
- `make`
- Python 3 + `numpy` + `matplotlib` (uniquement pour `make plot`)
- Valgrind (uniquement pour `make valgrind`) : `sudo apt install valgrind`

### 5.2 Commandes make

| Commande | Description |
|---|---|
| `make` | Compile en mode release (`-O2`). Produit `./simulation_disque`. |
| `make run` | Compile (si nécessaire) puis lance la simulation. |
| `make debug` | Compile avec `-g -O0` + AddressSanitizer + UndefinedBehaviorSanitizer. Utile pour détecter les corruptions mémoire et comportements indéfinis. |
| `make valgrind` | Compile en debug puis lance sous Valgrind (`--leak-check=full`). Plus lent qu'ASan mais ne nécessite pas de recompilation spéciale. |
| `make plot` | Lance `visualisation.py` (voir §9). |
| `make clean` | Supprime les `.o`, l'exécutable, `evolution_nuage.fich` et `diagnostics.fich`. |
| `make help` | Affiche le résumé des cibles disponibles. |

**Exemple de session complète :**

```bash
# Modifier les paramètres si besoin
nano params.h

# Compiler et lancer
make run

# Visualiser
make plot
```

---

## 6. Paramètres de simulation

Tous les paramètres sont définis dans **`params.h`** sous forme de constantes `constexpr`.
Il suffit de modifier ce fichier et de relancer `make` pour mettre à jour l'exécutable.

### 6.1 Paramètres physiques

| Paramètre | Défaut | Unité | Description |
|---|---|---|---|
| `MASSE_ETOILE` | 1 M☉ = 1.989×10³⁰ | kg | Masse de la protoétoile |
| `RAYON_ETOILE` | 60 UA | m | Rayon d'exclusion autour de l'étoile |
| `FRACTION_DISK` | 0.10 | — | Fraction de masse stellaire du disque |
| `MASSE_DISK` | 0.10 M☉ | kg | Masse totale du disque (`FRACTION_DISK × MASSE_ETOILE`) |
| `RAYON_NUAGE` | 150 UA | m | Rayon initial du nuage moléculaire |
| `BETA_ROT` | 0.4 | — | Rapport E_rot / \|E_grav\| du nuage initial |
| `V_THERM` | 0.0 | m/s | Agitation thermique initiale (σ gaussien). 0 = désactivée |

**Note sur BETA_ROT :** Les observations donnent β ~ 0.02–0.07 pour les cœurs denses réels.
Une valeur plus grande (0.1–0.5) accélère la formation du disque et est plus pratique pour
les tests numériques. Dans notre cas l'étoile est déjà formée et le choix d'un β élevé est cohérent.

### 6.2 Paramètres numériques

| Paramètre | Défaut | Description |
|---|---|---|
| `N_PARTICULES` | 1000 | Nombre de super-grains (hors étoile) |
| `DUREE_EN_TFF` | 15.0 | Durée totale de la simulation en unités de t_ff |
| `PAS_EN_TFF` | 200.0 | Nombre de pas RK4 par t_ff (plus grand → plus précis, plus lent) |
| `PAS_ECRITURE` | 5 | Écriture des positions toutes les N étapes RK4 |
| `PAS_DIAG` | 1 | Écriture des diagnostics toutes les N étapes RK4 |
| `EPS_FACTOR` | 0.25 | Facteur de softening : ε = EPS_FACTOR × R_nuage / N^(1/3) |
| `ENABLE_GRAIN_GRAIN` | false | Active la gravité mutuelle grain↔grain (O(N²)) |

**Choix de N_PARTICULES :**
- N = 100–500 : tests rapides (quelques secondes à quelques minutes)
- N = 1000 : simulation standard (quelques minutes sans grain-grain)
- N > 2000 avec grain-grain : très lent (O(N²) ≈ 4×10⁶ paires)

**Choix de PAS_EN_TFF :**
Pour les orbites keplériennes les plus rapides (rayon ~ RAYON_ETOILE), la période est
T_K ~ 2π / Ω_K(R_etoile). Le pas doit résoudre cette période : `h < T_K / 10`.
Augmenter PAS_EN_TFF si les énergies divergent.

### 6.3 Dissipations

| Paramètre | Défaut | Description |
|---|---|---|
| `T_DISS_START` | 1.0 | Activation des dissipations (en t_ff). Mettre très grand pour les désactiver. |
| `F_DISS_DISK` | 1.5 | Coefficient d'amortissement vertical (sans dim., typique : 1–3) |
| `F_CIRC` | 0.05 | Taux de circularisation. t_e ≈ 1/(f_circ · Ω_K). |
| `T0_DISK` | 150.0 | K | Température du disque à R = 1 UA |
| `MU_MOL_DISK` | 2.0×10⁻³ | kg/mol | Masse molaire du gaz (H₂) |

**Note sur T_DISS_START :** Les dissipations n'ont de sens physique qu'une fois un disque
plan formé (t > quelques t_ff). Les activer trop tôt biaiserait l'effondrement sphérique.
Les désactiver complètement (`T_DISS_START` très grand) permet de valider la conservation
de E_tot et L_z.

### 6.4 Pression (méthode PM)

| Paramètre | Défaut | Description |
|---|---|---|
| `ENABLE_PRESSURE` | true | Active le calcul de pression PM |
| `N_GRID` | 16 | Résolution Ng de la grille (Ng³ cellules). Typique : 8–32 |
| `CS_PM` | 100.0 | m/s | Vitesse du son isotherme. À 1 UA dans H₂ à 150 K : ~700 m/s |
| `GAMMA_PM` | 1.0 | Indice polytropique (1 = isotherme, 5/3 = adiabatique) |
| `R_DISK_PM` | 0.75 × RAYON_NUAGE | m | Rayon du cylindre de la grille PM |
| `H_DISK_PM` | 0.15 × RAYON_NUAGE | m | Demi-hauteur du cylindre de la grille PM |

**Note sur N_GRID :** Augmenter Ng améliore la résolution spatiale de la pression mais
peut introduire des instabilités numériques (la pression réagit à des fluctuations de
densité de plus en plus petites). Ng = 16 est un bon compromis de départ.

---

## 7. Fichiers de sortie

### `evolution_nuage.fich`

Trajectoires de toutes les particules, écrites toutes les `PAS_ECRITURE` étapes.

Format (colonnes séparées par des espaces) :

```
t (s)   i   x (m)   y (m)   z (m)   vx (m/s)   vy (m/s)   vz (m/s)
```

- Ligne `i=0` : étoile centrale (toujours à l'origine)
- Lignes `i=1..N` : grains du disque

### `diagnostics.fich`

Diagnostics physiques, écrits toutes les `PAS_DIAG` étapes.

Format :

```
t (s)   E_cin (J)   E_grav (J)   E_tot (J)   L_z (kg·m²/s)   dLz_rel
```

- `dLz_rel = (L_z(t) - L_z(0)) / |L_z(0)|` : variation relative du moment cinétique
  (indicateur de précision numérique).

---

## 8. Validation

### Conservation de l'énergie et du moment cinétique

Sans dissipation (`T_DISS_START` très grand, `ENABLE_PRESSURE = false`) :

- `E_tot = E_cin + E_grav` doit être **conservée** : vérifier dans `diagnostics.fich`
  que la variation relative reste inférieure à ~1%.
- `L_z` doit être **conservée** : `dLz_rel` doit rester petit (< ~0.1–1% selon N et h).

Si ces quantités dérivent significativement, réduire le pas de temps (`PAS_EN_TFF` plus grand).

### Avec dissipation

- `L_z` doit rester conservée (la circularisation conserve le moment cinétique).
- `E_tot` doit décroître monotonement (dissipation irréversible).
- Les particules doivent progressivement s'établir sur des orbites circulaires dans le
  plan z = 0.

---


## 9. Limitations et pistes d'amélioration

### Limitations actuelles

- **Complexité O(N²)** pour la gravité grain↔grain : impraticable au-delà de N ~ 2000.
  Solution : implémenter un arbre de Barnes-Hut (O(N log N)) ou une méthode PM pour
  la gravité (comme pour la pression).

- **Pas de temps fixe** : inefficace si les orbites ont des périodes très différentes
  (orbites proches de l'étoile très rapides, orbites externes lentes).
  Solution : intégrateur adaptatif (RK45, DOP853) ou pas de temps individuel par
  particule.

- **Grille PM cubique** : la grille PM couvre un cube mais le disque est cylindrique.
  Des particules à grand R_cyl mais z ~ 0 peuvent sortir de la grille sans recevoir de
  force de pression.

- **Pas de collisions** : les particules se traversent (pas de section efficace de
  collision). Acceptable pour de la dynamique orbitale, non pour la coagulation.

### Pistes d'amélioration

- Amélioration du stochage de donnés (pipe entre python et c++)
- Arbre de Barnes-Hut pour la gravité N-corps
- Intégrateur symplectique (Leapfrog, Yoshida) mieux adapté aux systèmes hamiltoniens
- Pas de temps adaptatif
- Loi de pression plus réaliste (EOS avec refroidissement radiatif)
- Parallélisation OpenMP des boucles de forces


---
