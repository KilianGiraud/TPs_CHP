# TP3 – Générateurs pseudo-aléatoires de qualité et intervalles de confiance

> Version anglaise : [TP3.md](TP3.md)

Deux générateurs de nombres pseudo-aléatoires de qualité, Mersenne Twister
MT19937 et MRG32k3a de L'Ecuyer, sont vérifiés par rapport à des valeurs de
référence puis chronométrés. Ils servent ensuite à estimer le volume d'une
sphère, et π, à l'aide de réplications Monte Carlo indépendantes calculées en
parallèle (OpenMP), avec un intervalle de confiance autour de la moyenne (loi
de Student).

## Principe

### Q1 – Générateurs

| Générateur | Auteur | Sortie | Période | État |
|---|---|---|---|---|
| MT19937 | Matsumoto & Nishimura (1998) | `genrand_real1` : [0, 1] | 2¹⁹⁹³⁷ − 1 | 624 × 32 bits |
| MRG32k3a | L'Ecuyer (1999) | ]0, 1[ | ≈ 2¹⁹¹ | 6 × 32 bits |

Les deux générateurs sont écrits sous forme réentrante : l'état est une
structure passée en paramètre au lieu de variables globales, pour que chaque
réplication (thread) possède son propre générateur. Les algorithmes eux-mêmes
sont inchangés et sont vérifiés par rapport à des valeurs de référence au
lancement de `benchRNG` (voir Résultats).

### Q2 – Volume de la sphère

Des points `(x, y, z)` sont tirés uniformément dans le cube `[0,R]³` (un
octant, par symétrie). La proportion `f` de points vérifiant
`x² + y² + z² ≤ R²` converge vers `(4/3 π R³ / 8) / R³ = π/6`, donc

```
V = 4/3 π R³ ≈ 8 R³ f          π ≈ 6 f
```

Avec `dim = 2`, le programme calcule le quart de disque du TP2 (`π ≈ 4 f`).
Écart-type d'une estimation de π avec `N` points (loi binomiale) :

```
3D : σ = 6 √(p(1−p)/N),  p = π/6   ≈ 3.00 / √N
2D : σ = 4 √(p(1−p)/N),  p = π/4   ≈ 1.64 / √N
```

L'estimateur 3D de π est 1.8 fois plus bruité que l'estimateur 2D et demande
3 nombres aléatoires par point au lieu de 2.

### Q2 – Intervalle de confiance (annexe du sujet)

Pour `n` réplications indépendantes `X1 … Xn` (chacune avec son propre flux
aléatoire) :

```
X(n)   = Σ Xi / n                              moyenne
S²(n)  = Σ (Xi − X(n))² / (n − 1)              variance sans biais
S(n)/√n                                        erreur standard de la moyenne
R      = t(n−1, 1−α/2) × √(S²(n) / n)          rayon de confiance
IC     = [X(n) − R ; X(n) + R]                 au niveau 1 − α = 95 %
```

`t(n−1, 0.975)` est lu dans la Table 1 du sujet avec **ν = n − 1 degrés de
liberté** (voir Remarques). Un intervalle à 99 % (α = 0.01, `t(n−1, 0.995)`)
est également affiché.

## Fichiers

| Fichier | Question | Description |
|---|---|---|
| `mt19937ar.h` | Q1-A | Mersenne Twister MT19937, version réentrante du `mt19937ar.c` original |
| `mrg32k3a.h` | Q1-B | MRG32k3a issu de `combmrg2.c` (réentrant) + saut de 2¹²⁷ pas pour obtenir des flux indépendants |
| `benchRNG.c` | Q1-C | Vérification des deux générateurs par rapport aux valeurs de référence, puis chronométrage de 10⁹ appels (+ `rand()` pour comparaison) |
| `RunBench` | Q1-C | Lance le benchmark compilé en `-O0`, `-O2` et `-O3 -march=native` |
| `sphereCI.c` | Q2 | Volume de la sphère et π, réplications en parallèle (OpenMP), intervalles de confiance à 95 % et 99 % |
| `RunCI` | Q2 | 10 et 40 réplications × 10⁹ points avec les deux générateurs + deux petites études de sensibilité |
| `results/` | Q1, Q2 | Sorties de `RunBench` et `RunCI` |
| `openMPcomputePI.c` | – | Fourni : exemple OpenMP en 2D (générateur xorshift) |
| `combmrg2.c` | – | Fourni : MRG32k3a, MRG32k5a et MRG63k3a (référence uniquement, non compilable tel quel) |

## Prérequis

- Linux (ou macOS / WSL sous Windows)
- `gcc` avec le support d'OpenMP (`libgomp`, installé avec gcc sur la plupart
  des distributions)
- Optionnel : `htop` pour observer les threads

## Compilation

```bash
gcc -O0 -Wall -o benchRNG_O0 benchRNG.c
gcc -O2 -Wall -o benchRNG_O2 benchRNG.c
gcc -O3 -march=native -Wall -o benchRNG_O3 benchRNG.c
gcc -O2 -Wall -fopenmp -o sphereCI sphereCI.c -lm
chmod +x RunBench RunCI
```

Les générateurs sont des fichiers d'en-tête (fonctions `static inline`) : il
n'y a rien d'autre à compiler ni à lier. Les binaires sont listés dans
`.gitignore`.

## Utilisation

### Q1 – Générateurs : vérifications et performances

```bash
./benchRNG_O2          # vérifications de référence + 10^9 appels par générateur
./benchRNG_O2 1e8      # <nbCalls>
./RunBench             # -O0, -O2, -O3 l'un après l'autre -> results/bench_O*.txt
```

Les trois exécutions de `RunBench` sont volontairement séquentielles : deux
benchmarks lancés en même temps se partageraient le CPU et fausseraient les
mesures.

### Q2 – Intervalles de confiance

```bash
./sphereCI mt 10                   # MT19937, 10 réplications, 10^9 points, 3D
./sphereCI mrg 40                  # MRG32k3a, 40 réplications
./sphereCI mt 10 1e6               # <mt|mrg> <nbReplicates> [nbPoints] [dim]
./sphereCI mt 10 1e6 2             # 2D (quart de disque, comme au TP2)
OMP_NUM_THREADS=4 ./sphereCI mt 10 # choix du nombre de threads
./RunCI                            # toutes les exécutions des Résultats -> results/ (~8 min)
```

- `nbReplicates` : de 2 à 50 (la table de Student est utilisée jusqu'à ν = 40).
- `nbPoints` accepte la notation scientifique (`1e9`).
- Chaque réplication possède son générateur, dont le flux ne dépend que du
  numéro de réplication :
  - `mt` : `init_by_array` avec la clé originale de Makoto Matsumoto
    `{0x123, 0x234, 0x345, 0x456}`, complétée par le numéro de réplication en
    5ᵉ mot ;
  - `mrg` : la réplication *r* utilise le flux *r*, c'est-à-dire la graine par
    défaut `{12345 × 6}` avancée de *r* × 2¹²⁷ pas.
- Les résultats sont donc **identiques quel que soit le nombre de threads**
  (vérifié avec 1, 5 et 12 threads). La boucle utilise
  `schedule(dynamic, 1)` : un thread prend la réplication suivante dès qu'il
  est libre.

## Résultats

Machine : AMD Ryzen 5 5600H (6 cœurs / 12 threads), gcc 16.2.1.

### Q1 – Vérifications de référence

| Vérification | Attendu (référence) | Résultat |
|---|---|---|
| MT `init_by_array({0x123, 0x234, 0x345, 0x456})`, 5 premières sorties | `1067595299 955945823 477289528 4107218783 4228976476` (`mt19937ar.out`) | OK |
| MT `init_genrand(5489)`, 10000ᵉ sortie | `4123659995` (norme C++11, `std::mt19937`) | OK |
| MRG32k3a, graine `{12345 × 6}`, première sortie | `0.1270111220` (calcul exact en entiers) | OK |
| MRG32k3a, saut de 2¹²⁷ pas | 2ᵉ flux de RngStreams : `{3692455944, 1366884236, 2968912127, 335948734, 4161675175, 475798818}` | OK |

Les vérifications passent aux trois niveaux d'optimisation, et les moyennes
des tirages sont identiques (à la précision affichée) entre `-O0`, `-O2` et
`-O3` : l'optimisation change la vitesse, pas les séquences aléatoires.

### Q1 – Performances, 10⁹ appels (ns par appel = temps total en s)

| Générateur | `-O0` | `-O2` | `-O3 -march=native` |
|---|---|---|---|
| `rand()` (TP2, référence) | 7.43 | 7.40 | 7.40 |
| MT19937 `genrand_real1` | 9.15 | 4.30 | 2.05 |
| MRG32k3a | 48.74 | 28.71 | 27.14 |

- **MT19937** gagne ×2.1 en `-O2` et ×4.5 en `-O3 -march=native` : une fois
  optimisé, c'est de loin le plus rapide, comme l'annonce le sujet.
- **`rand()`** ne change pas : il se trouve dans la bibliothèque C, déjà
  compilée, donc nos options de compilation ne s'appliquent pas à lui.
- **MRG32k3a** ne gagne que ×1.7 et reste 5 à 13 fois plus lent que MT :
  chaque appel enchaîne 2 divisions et 3 branchements qui dépendent des
  données (voir Remarques).

### Q2 – Intervalles de confiance, 10⁹ points par réplication (sphère 3D)

σ d'une réplication (théorie) : 9.48·10⁻⁵. π exact = 3.1415926536.

| Générateur | n | π moyen | S(n) | t(n−1 ; 0.975) | rayon R à 95 % | IC à 95 % | π ∈ IC | temps |
|---|---|---|---|---|---|---|---|---|
| MT19937 | 10 | 3.1416330228 | 5.99·10⁻⁵ | 2.262 | 4.29·10⁻⁵ | [3.1415901724 ; 3.1416758732] | oui | 15.1 s |
| MT19937 | 40 | 3.1415992318 | 7.92·10⁻⁵ | 2.023 | 2.53·10⁻⁵ | [3.1415738960 ; 3.1416245677] | oui | 57.5 s |
| MRG32k3a | 10 | 3.1415737494 | 8.99·10⁻⁵ | 2.262 | 6.43·10⁻⁵ | [3.1415094622 ; 3.1416380366] | oui | 82.9 s |
| MRG32k3a | 40 | 3.1415963326 | 1.07·10⁻⁴ | 2.023 | 3.41·10⁻⁵ | [3.1415622124 ; 3.1416304529] | oui | 329.9 s |

Volume de la sphère avec 40 réplications (exact : 4/3 π = 4.1887902048) :
MT 4.1887989758 ± 3.38·10⁻⁵, MRG32k3a 4.1887951102 ± 4.55·10⁻⁵ (95 %).

- **Plus de réplications ⇒ rayon de confiance plus petit.** De 10 à 40
  réplications, le rayon est divisé par 1.7 (MT) et par 1.9 (MRG). Avec σ
  connu, le facteur attendu est (2.262 / 2.023) × √(40/10) = 2.24 (rayons
  théoriques 6.78·10⁻⁵ et 3.03·10⁻⁵). L'écart vient de S(n) : avec 10
  réplications, c'est une estimation bruitée de σ (5.99·10⁻⁵ pour MT, basse
  par hasard) ; avec 40, elle se rapproche de 9.48·10⁻⁵. C'est précisément
  cette incertitude sur S que la loi de Student prend en compte (t = 2.262 au
  lieu de 1.960).
- Les 4 intervalles contiennent π. Avec 40 réplications (4·10¹⁰ points),
  l'erreur sur π vaut 6.6·10⁻⁶ (MT) et 3.7·10⁻⁶ (MRG32k3a).
- Les deux générateurs donnent des résultats statistiquement équivalents ;
  MRG32k3a prend 5.5 à 5.7 fois plus de temps.

### Q2 – Effet du nombre de réplications (MT, 10⁶ points chacune)

| n | t(n−1 ; 0.975) | rayon R à 95 % | R avec σ théorique | IC à 95 % | π ∈ IC |
|---|---|---|---|---|---|
| 5 | 2.776 | 2.47·10⁻³ | 3.72·10⁻³ | [3.1383042930 ; 3.1432469070] | oui |
| 10 | 2.262 | 2.50·10⁻³ | 2.14·10⁻³ | [3.1397010252 ; 3.1446965748] | oui |
| 20 | 2.093 | 1.30·10⁻³ | 1.40·10⁻³ | [3.1405810524 ; 3.1431841476] | oui |
| 30 | 2.045 | 1.10·10⁻³ | 1.12·10⁻³ | [3.1404858257 ; 3.1426841743] | oui |
| 40 | 2.023 | 8.94·10⁻⁴ | 9.59·10⁻⁴ | [3.1404780850 ; 3.1422656150] | oui |
| 50 | 2.021 | 7.47·10⁻⁴ | 8.56·10⁻⁴ | [3.1405881746 ; 3.1420816654] | oui |

Le rayon décroît comme `t(n−1) / √n`. Avec peu de réplications, il fluctue
(de 5 à 10 : le S des 5 premières réplications est bas par hasard) : en
dessous d'une dizaine de réplications, `t` et l'incertitude sur S sont tous
deux grands. Les réplications 1 … n sont les mêmes flux pour chaque n : ces
intervalles sont donc des échantillons emboîtés, et non des expériences
indépendantes.

### Q2 – Effet du nombre de points par réplication (MT, 10 réplications)

| N | rayon R à 95 % | R avec σ théorique | IC à 95 % | π ∈ IC |
|---|---|---|---|---|
| 10³ | 7.59·10⁻² | 6.78·10⁻² | [3.0614978493 ; 3.2133021507] | oui |
| 10⁴ | 2.24·10⁻² | 2.14·10⁻² | [3.1130336956 ; 3.1578063044] | oui |
| 10⁵ | 6.95·10⁻³ | 6.78·10⁻³ | [3.1329965676 ; 3.1468914324] | oui |
| 10⁶ | 2.50·10⁻³ | 2.14·10⁻³ | [3.1397010252 ; 3.1446965748] | oui |
| 10⁷ | 6.61·10⁻⁴ | 6.78·10⁻⁴ | [3.1407424139 ; 3.1420633861] | oui |
| 10⁸ | 1.57·10⁻⁴ | 2.14·10⁻⁴ | [3.1412230140 ; 3.1415368140] | **non** |
| 10⁹ | 4.29·10⁻⁵ | 6.78·10⁻⁵ | [3.1415901724 ; 3.1416758732] | oui |

- Le rayon est divisé par environ √10 ≈ 3.16 à chaque fois que N est
  multiplié par 10 : c'est la convergence lente en `1/√N` du Monte Carlo, un
  chiffre correct de plus coûte 100 fois plus de points.
- **L'intervalle à 10⁸ points ne contient pas π, et c'est attendu** : un
  intervalle à 95 % manque la vraie valeur dans 5 % des cas en moyenne (1 fois
  sur 20). Sur les 15 intervalles distincts de ce TP, 1 la manque. Utiliser
  plus de points ne donne pas « plus de certitude », seulement un intervalle
  plus étroit au même niveau de confiance de 95 %.

### Q2 – Accélération parallèle (12 réplications × 5·10⁷ points, 3D)

```bash
for gen in mt mrg; do
    for t in 1 6 12; do OMP_NUM_THREADS=$t ./sphereCI $gen 12 5e7 | grep time; done
done
```

| Threads | MT19937 | accélération | MRG32k3a | accélération |
|---|---|---|---|---|
| 1 | 6.26 s | 1.0 | 50.24 s | 1.0 |
| 6 | 1.05 s | 6.0 | 8.36 s | 6.0 |
| 12 | 0.82 s | 7.6 | 4.18 s | 12.0 |

Avec un thread par cœur physique (6), l'accélération est idéale pour les deux
générateurs : les réplications sont totalement indépendantes. Les 6 threads
matériels supplémentaires (SMT) apportent peu à MT, qui occupe déjà les
unités d'exécution, mais doublent la vitesse de MRG32k3a : son code passe
l'essentiel de son temps à attendre (divisions, branchements mal prédits), et
le second thread matériel de chaque cœur utilise les unités inoccupées.

## Remarques

- **Lecture de la Table 1 du sujet** : la première colonne est le nombre de
  degrés de liberté ν = n − 1, et non le nombre de réplications n (12.706 est
  la valeur de t pour ν = 1, c'est-à-dire 2 réplications ; il n'y a pas
  d'intervalle avec une seule réplication). Avec n = 10 réplications, la bonne
  valeur est donc t(9) = 2.262, et non 2.228.
  Entre ν = 30 et ν = 40, la table est interpolée linéairement en 1/ν
  (t(39) = 2.023, valeur exacte 2.0227) ; de ν = 40 à 49, on utilise la valeur
  de ν = 40 (2.021), légèrement supérieure à la valeur exacte, ce qui rend
  l'intervalle légèrement conservateur. Le nombre de réplications est limité à
  50. La valeur de la table pour ν = 8 (2.308) est conservée telle quelle, la
  valeur exacte étant 2.306.
- **α = 0.01** : la table pour α = 0.01 demandée dans l'annexe se trouve dans
  `sphereCI.c` (`T995`, de 63.657 pour ν = 1 à 2.704 pour ν = 40) ;
  l'intervalle à 99 % est affiché après celui à 95 %.
- **Flux indépendants** : pour MRG32k3a, les flux sont espacés de 2¹²⁷ tirages
  (matrices de saut du package RngStreams de L'Ecuyer) et une réplication
  utilise 3·10⁹ ≈ 2³¹·⁵ tirages : les flux ne peuvent pas se chevaucher,
  l'indépendance est garantie par construction. Pour MT19937, les réplications
  utilisent des clés d'initialisation différentes : elles sont pratiquement
  indépendantes (période 2¹⁹⁹³⁷ − 1, `init_by_array` mélange l'intégralité de
  l'état de 624 mots), mais l'absence de chevauchement n'est pas formellement
  garantie. La technique rigoureuse pour MT est le saut (jump-ahead) ou
  l'utilisation d'états précalculés et sauvegardés (sequence splitting).
- **double, et non float** : tous les calculs se font en `double` et les
  compteurs sont des `long long`. Un compteur `float` cesse d'augmenter à
  2²⁴ = 16 777 216 (ajouter 1 ne le modifie plus), bien en dessous des 10⁹
  points d'une réplication ; et avec ~7 chiffres significatifs, un `float` est
  à la limite de la précision de 10⁻⁶ atteinte sur π.
- **Variance en deux passes** : S²(n) est calculée à partir des écarts à la
  moyenne, comme dans l'annexe, et non avec `Σ Xi²/n − X(n)²`, qui soustrait
  deux nombres ≈ 9.87 pour obtenir une variance ≈ 10⁻⁸ : environ la moitié des
  chiffres significatifs d'un `double` serait perdue (et tous avec un
  `float`).
- **Pourquoi MRG32k3a est lent ici** : chaque appel enchaîne 2 divisions en
  `double` et des conversions `double ↔ long`, plus 3 branchements qui
  dépendent des valeurs aléatoires (`if (p < 0)` deux fois et la combinaison
  finale), très souvent mal prédits. Une variante de test où ces branchements
  sont remplacés par des masques (même récurrence, sorties identiques bit à
  bit) tourne à 13.6 ns par appel au lieu de 28.1 ns en `-O2`. Le code fourni
  est conservé tel quel dans ce TP.
- **`-march=native`** : le binaire `-O3` utilise les instructions du CPU sur
  lequel il a été compilé (AVX2, FMA…) et peut ne pas fonctionner sur une
  autre machine ; `-O2` est le choix portable.
- **`M_PI`** : `sphereCI.c` définit `_DEFAULT_SOURCE` / `_USE_MATH_DEFINES`
  avant d'inclure `<math.h>`, comme au TP2.

## Dépannage

| Problème | Solution |
|---|---|
| `undefined reference to 'omp_get_wtime'` / `omp.h: No such file` | Compiler avec `-fopenmp` (installer `libgomp` si besoin) |
| `undefined reference to 'sqrt'` / `'pow'` | Ajouter `-lm` à la fin de la commande gcc |
| `Permission non accordée` / `Permission denied` | `chmod +x RunBench RunCI` |
| `./benchRNG_O0: No such file or directory` | Compiler d'abord les binaires (section Compilation), dans le répertoire `TP3` |
| `Illegal instruction` avec `benchRNG_O3` | Compilé avec `-march=native` sur un autre CPU : le recompiler sur cette machine |
| `Usage: ./sphereCI ...` | Générateur `mt` ou `mrg`, 2 à 50 réplications, dim 2 ou 3 |
| `Reference check FAILED` | Un générateur a été modifié : comparer avec `mt19937ar.h` / `mrg32k3a.h` |
| `/bin/bash^M: bad interpreter` | Fins de ligne Windows : `sed -i 's/\r$//' RunBench RunCI` |
