# TP2 – Simulation Monte Carlo et réplications parallèles

> Version anglaise : [TP2.md](TP2.md)

Estimation de π par la méthode de Monte Carlo : d'abord en une seule
exécution, puis avec des réplications indépendantes exécutées séquentiellement,
et enfin avec des réplications lancées en parallèle sur une machine SMP grâce
au multitâche Unix (modèle SPMD).

## Principe

Des points `(xr, yr)` sont tirés uniformément dans le carré unité `[0,1]²`. La
proportion de points tombant dans le quart de disque `x² + y² ≤ 1` converge
vers `π/4`, donc

```
π ≈ 4 × inDisk / nbPoints
```

Chaque point nécessite deux nombres aléatoires. L'écart-type d'une estimation
vaut `σ = √(π(4−π)/n) ≈ 1.64/√n` : chaque chiffre correct supplémentaire coûte
100 fois plus de points.

## Fichiers

| Fichier | Question | Description |
|---|---|---|
| `simPi.c` | Q1 | Estimation unique de π avec `rand()`, tests avec 10³, 10⁶, 10⁹ points et étude de précision |
| `replicates.c` | Q2 | 10 réplications indépendantes exécutées séquentiellement avec Mersenne Twister, moyenne, erreurs absolue et relative |
| `simPI.c` | Q3 | Une réplication, graine donnée en ligne de commande, résultat affiché sur stdout |
| `Sim10times` | Q3 | Script shell lançant 10 processus `simPI` en parallèle avec des graines différentes |

## Prérequis

- Linux (ou macOS / WSL sous Windows)
- `gcc`
- Optionnel : `htop` ou `top` pour observer les processus

## Compilation

```bash
gcc -O2 -Wall -o simPi      simPi.c      -lm
gcc -O2 -Wall -o replicates replicates.c -lm
gcc -O2 -Wall -o simPI      simPI.c
chmod +x Sim10times
```

`-O2` est important : sans optimisation, les exécutions à 10⁹ points sont
plusieurs fois plus lentes.

## Utilisation

### Q1 – Estimation unique

```bash
./simPi            # tests avec 1e3, 1e6, 1e9 points + étude de précision
./simPi 1000000    # une seule exécution avec le nombre de points donné
```

### Q2 – Réplications séquentielles

```bash
./replicates              # 10 réplications pour 1e3, 1e6 et 1e9 points
./replicates 10 1000000   # <nbReplicates> <nbPoints>
```

Le générateur est initialisé **une seule fois** avec la clé originale du code
de référence de Makoto Matsumoto (`{0x123, 0x234, 0x345, 0x456}`) et n'est
jamais réinitialisé entre deux réplications : chaque réplication utilise la
suite du même flux aléatoire, ce qui garantit l'indépendance.

### Q3 – Réplications parallèles (SPMD)

```bash
./simPI 19937              # une réplication, 1e9 points par défaut
./simPI 19937 1000000      # <seed> [nbPoints]
./simPI 19937 > PI1        # résultat redirigé dans un fichier

./Sim10times               # 10 réplications en parallèle -> PI1 … PI10
time ./Sim10times          # mesure du temps d'exécution en parallèle
```

Chaque ligne de `Sim10times` lance `./simPI <seed> > PIi &` en arrière-plan ;
Linux répartit les 10 processus sur les cœurs disponibles. Un `wait` final
fait que le script ne se termine que lorsque toutes les réplications sont
finies.

Afficher les résultats et une moyenne rapide :

```bash
for i in $(seq 1 10); do echo "PI$i: $(cat PI$i)"; done
cat PI* | awk '{s+=$1; n++} END {printf "mean = %.10f\n", s/n}'
```

## Résultats

### Q1 – Précision en fonction du nombre de points (théorie)

| Précision visée | n (1σ, ~68 %) | n (2σ, ~95 %) |
|---|---|---|
| 10⁻² (3.14) | ≈ 2.7·10⁴ | ≈ 1.1·10⁵ |
| 10⁻³ | ≈ 2.7·10⁶ | ≈ 1.1·10⁷ |
| 10⁻⁴ (3.1415) | ≈ 2.7·10⁸ | ≈ 1.1·10⁹ |

### Q2 – 10 réplications séquentielles (Mersenne Twister)

| nbPoints | meanPI | erreur absolue | erreur relative |
|---|---|---|---|
| 10³ | 3.1480000 | 6.4·10⁻³ | 2.0·10⁻³ |
| 10⁶ | 3.1412544 | 3.4·10⁻⁴ | 1.1·10⁻⁴ |
| 10⁹ | 3.1415644 | 2.8·10⁻⁵ | 9.0·10⁻⁶ |

### Q3 – 10 réplications parallèles, 10⁹ points chacune

| Graine | Fichier | π estimé |
|---|---|---|
| 19937 | PI1 | 3.1415328240 |
| 52149 | PI2 | 3.1416400000 |
| 31415 | PI3 | 3.1416069640 |
| 27182 | PI4 | 3.1416070920 |
| 14142 | PI5 | 3.1416196400 |
| 17320 | PI6 | 3.1416013360 |
| 44257 | PI7 | 3.1415138560 |
| 70913 | PI8 | 3.1415850600 |
| 63113 | PI9 | 3.1416275640 |
| 86226 | PI10 | 3.1416601640 |

| Statistique | Valeur |
|---|---|
| π moyen | 3.1415994500 |
| erreur absolue | 6.8·10⁻⁶ |
| erreur relative | 2.2·10⁻⁶ |
| écart-type | ≈ 4.5·10⁻⁵ (théorie : 5.2·10⁻⁵) |
| erreur standard de la moyenne | ≈ 1.4·10⁻⁵ |

Les résultats sont reproductibles bit à bit : les mêmes graines donnent les
mêmes valeurs sur n'importe quelle machine.

## Remarques

- **`rand()` n'est utilisé qu'en Q1**, comme première approche. Il est
  statistiquement faible et lent (la glibc prend un verrou à chaque appel) ;
  Q2 et Q3 utilisent Mersenne Twister MT19937 (période 2¹⁹⁹³⁷ − 1).
- **L'accélération parallèle dépend du nombre de cœurs logiques** (`nproc`).
  Avec 10 cœurs ou plus, les 10 réplications prennent à peu près le temps
  d'une seule ; avec moins, le temps total vaut environ `10 / nproc` fois
  celui d'une exécution.
- **Indépendance des flux parallèles** : des graines différentes avec
  `init_genrand` donnent des flux pratiquement indépendants, mais l'absence de
  chevauchement n'est pas formellement garantie. L'approche rigoureuse consiste
  à utiliser des états du générateur précalculés et largement espacés
  (jump-ahead), un par processus.
- **`M_PI` non défini** : `M_PI` est défini par POSIX, pas par l'ISO C.
  `simPi.c` définit une valeur de repli, `replicates.c` définit
  `_DEFAULT_SOURCE` / `_USE_MATH_DEFINES` avant d'inclure `<math.h>`.

## Dépannage

| Problème | Solution |
|---|---|
| `Permission non accordée` / `Permission denied` | `chmod +x Sim10times` |
| Toujours refusé après `chmod` | Le répertoire est peut-être monté en `noexec` (vérifier avec `findmnt -T .`) : copier les fichiers dans un répertoire local |
| `./simPI: No such file or directory` | Compiler `simPI.c` dans le même répertoire que le script |
| `/bin/bash^M: bad interpreter` | Fins de ligne Windows : `sed -i 's/\r$//' Sim10times` |
