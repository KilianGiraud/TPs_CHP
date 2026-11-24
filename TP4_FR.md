# TP4 – Attention à l'arithmétique flottante

> Version anglaise : [TP4.md](TP4.md)

Chaque exemple de code du sujet est exécuté et expliqué (un programme par
section, en C ou en C++ comme dans le sujet), avec des expériences
supplémentaires pour comprendre les résultats. La partie « To go a little
further » est ensuite traitée : un softmax stable pour les étudiants en IA et,
pour le HPC, la sommation de Kahan et d'autres sommes compensées, mesurées sur
de très grands vecteurs aléatoires (résultats et temps de calcul), ainsi que la
reproductibilité d'une réduction parallèle avec OpenMP. Ce fichier constitue la
trace du TP demandée à la question 4.

## Sommaire

| Section du sujet | Programme | Observation principale |
|---|---|---|
| How does it impact my computing? | `compare.c` | `0.7f` promu en double est < 0.7 et 0.5 est exact : « only 1 is right » ; le sens de l'arrondi dépend de la valeur (0.6f > 0.6) |
| Don't test for equality, Relative errors | `equality.cpp` | `sqrt(10)²` = 10 + 1 ulp ; une tolérance absolue est fausse pour les grands comme pour les tout petits nombres |
| Worry about addition and subtraction (1) | `sumMul.c` | 10⁷ additions de 0.1f donnent 1 087 937 (+8.8 %), 10⁸ donnent 2 097 152 (−79 %) ; un produit n'est arrondi qu'une fois |
| Worry about addition and subtraction (2) | `derivative.cpp` | sortie du sujet reproduite chiffre pour chiffre ; valeurs négatives de gcc expliquées (x87) ; pas complexe : exact |
| Floating point numbers have finite ranges | `ranges.cpp` | 2²⁴ + 1 et 2⁵³ + 1 sont des cas d'égalité arrondis au pair ; DBL_EPSILON est l'écart après 1 |
| Use logarithms to avoid overflow and underflow | `logFactorial.cpp` | 200!/(190! 10!) = NaN en direct, erreur relative de 10⁻¹³ avec les logarithmes, exact avec un produit de rapports |
| Numeric operations don't always return numbers | `notNumbers.cpp` | inf, NaN, drapeaux IEEE ; `x - x == 0` échoue aussi pour inf ; un NaN ne se propage pas toujours |
| Out of order execution | `associativity.c` | (a + b) + c ≠ a + (b + c) ; gcc fusionne a·b + c en FMA avec `-mfma` ; ordre d'une longue somme |
| Pour aller plus loin – IA | `softmax.c` | le softmax naïf donne NaN pour des logits de ±1000 ; soustraire le max ; log-softmax |
| Pour aller plus loin – HPC, Q1 à Q3 | `summation.c`, `exactSum.c`, `sumBench.c` | Kahan exact au dernier bit sur 10⁷–10⁸ valeurs, 3.9× plus lent sur un cœur, gratuit avec 8 accumulateurs ou avec des threads ; variantes jusqu'à la somme exacte |
| Phase de réduction (HPC) | `reproSum.c` | une réduction OpenMP change avec le nombre de threads ET d'une exécution à l'autre ; les réductions compensées ou exactes non |

## Fichiers

| Fichier | Section / question | Description |
|---|---|---|
| `compare.c` | How does it impact… | Les 2 codes, les 2 corrections, champs IEEE et valeurs exactes de 0.1, 0.5, 0.6, 0.7 |
| `equality.cpp` | Égalité, erreurs relatives | sqrt(10)², ulps, tolérance absolue ou relative, `almostEqual` |
| `sumMul.c` | Addition (1) | 10 × 0.1f contre 0.1f × 10, puis n = 10 … 10⁸ additions |
| `derivative.cpp` | Soustraction (2) | Dérivée de sin : différence avant, centrée, pas complexe ; reproduction x87 |
| `ranges.cpp` | Plages finies | 2²⁴, 2⁵³, DBL_EPSILON, écart entre nombres consécutifs, limites de `float.h` |
| `logFactorial.cpp` | Logarithmes | 200!/(190! 10!) calculé de 4 façons, probabilité binomiale B(2000, ½) |
| `notNumbers.cpp` | Résultats non numériques | inf, NaN, drapeaux d'exception IEEE, comparaisons et propagation des NaN |
| `associativity.c` | Out of order execution | Non-associativité, FMA, série harmonique dans les deux sens |
| `softmax.c` | Pour aller plus loin (IA) | Softmax naïf et stable, log-softmax |
| `summation.h`, `summation.c` | Pour aller plus loin (HPC) Q2, Q3 | naïf, par paires, Kahan, Neumaier, Sum2, Klein, long double, versions à 8 accumulateurs et versions float |
| `exactSum.c` | Q2, Q3 | Somme exacte (Shewchuk / `math.fsum` de Python) : la référence |
| `mt19937ar.h`, `mt19937ar.c` | Q2 | MT19937 du TP3 + `genrand_res53`, découpé en `.h` / `.c` |
| `sumBench.c` | Q2, Q3 | Précision et temps de chaque algorithme sur 4 jeux de données (double, float), étude en fonction du conditionnement |
| `reproSum.c` | Phase de réduction | Ordre des additions, réductions OpenMP de 1 à 12 threads, exécutions répétées |
| `Makefile` | – | Compile tout, y compris les 3 binaires de démonstration |
| `RunAll` | – | Lance tout, sorties enregistrées dans `results/` |
| `results/` | – | Sorties de `./RunAll big` |
| `*.pdf` | – | Fournis : sujet, Goldberg (1991), guides de style C / C++ |

## Prérequis

- Linux (ou WSL / macOS avec gcc), `gcc` et `g++`, `make`
- OpenMP (`libgomp`, installé avec gcc sur la plupart des distributions) pour
  `reproSum`
- Un processeur x86-64 : les remarques sur le `long double` (x87) sont propres
  au x86 ; `associativity_fma` demande les instructions FMA (Intel depuis 2013,
  AMD depuis 2012)
- Environ 200 Mo de mémoire libre (1.3 Go pour `./sumBench 1e8`)

## Compilation

```bash
make            # tous les programmes
make clean      # supprime les binaires et les fichiers objets
```

Commandes équivalentes :

```bash
gcc -O2 -Wall -Wextra -o compare       compare.c       -lm
gcc -O2 -Wall -Wextra -o sumMul        sumMul.c        -lm
gcc -O2 -Wall -Wextra -o associativity associativity.c -lm
gcc -O2 -Wall -Wextra -o softmax       softmax.c       -lm
g++ -O2 -Wall -Wextra -o equality      equality.cpp
g++ -O2 -Wall -Wextra -o derivative    derivative.cpp
g++ -O2 -Wall -Wextra -o ranges        ranges.cpp
g++ -O2 -Wall -Wextra -o logFactorial  logFactorial.cpp
g++ -O2 -Wall -Wextra -o notNumbers    notNumbers.cpp
gcc -O2 -Wall -Wextra -c summation.c exactSum.c mt19937ar.c sumBench.c
gcc -o sumBench sumBench.o summation.o exactSum.o mt19937ar.o -lm
gcc -O2 -Wall -Wextra -fopenmp -o reproSum reproSum.c summation.o exactSum.o mt19937ar.o -lm

# binaires de démonstration, compilés exprès avec des options risquées
g++ -O2 -Wall -Wextra -ffast-math -o notNumbers_fast notNumbers.cpp
gcc -O2 -Wall -Wextra -mfma -o associativity_fma associativity.c -lm
gcc -O3 -ffast-math -Wall -Wextra -c -o summation_fast.o summation.c
gcc -o sumBench_fast sumBench.o summation_fast.o exactSum.o mt19937ar.o -lm
```

Tout est compilé en `-O2`, sans `-ffast-math` et sans `-march=native` :
arithmétique IEEE 754 stricte, mêmes résultats sur n'importe quelle machine
x86-64. Le `#warning` affiché à la compilation de `summation_fast.o` est
normal.

## Utilisation

```bash
./compare  ./equality  ./sumMul  ./derivative  ./ranges        # partie 1 : sans argument
./logFactorial  ./notNumbers  ./associativity  ./softmax
./notNumbers_fast  ./associativity_fma          # mêmes sources, options risquées

./sumBench           # N = 1e7 valeurs par jeu de données (~8 s)
./sumBench 1e8       # N = 1e8 valeurs (~70 s, 1.3 Go)
./sumBench_fast      # idem, summation.c compilé en -O3 -ffast-math
./reproSum           # N = 1e7 valeurs, OpenMP de 1 à 12 threads (~6 s)
./reproSum 1e8
./RunAll             # tout -> results/*.txt (~25 s)
./RunAll big         # + sumBench 1e8 (~1 min 35 s)
```

## Résultats et observations

Machine : AMD Ryzen 5 5600H (6 cœurs / 12 threads, Zen 3, AVX2 + FMA), 7 Go de
RAM, Fedora Linux, gcc / g++ 16.2.1. Toutes les sorties sont dans `results/`.

### 1. How does it impact my computing? – `compare.c`

Les deux codes affichent **`only 1 is right`**.

| Constante | stockée dans un float (32 bits) | stockée dans un double (64 bits) |
|---|---|---|
| 0.7 | 0.699999988079071044921875 | 0.6999999999999999555910790149937383830547332763671875 |
| 0.5 | 0.5 | 0.5 |

- `a < 0.7` : `0.7` est une constante **double**, donc `a` est promu en
  double. Le float est un arrondi plus grossier de 0.7 (24 bits au lieu de 53)
  et il est arrondi **vers le bas** : 0.69999998808 < 0.69999999999999995559
  → vrai.
- `b < 0.5` : 0.5 = 2⁻¹ est exact dans les deux formats → faux.
- Avec `.7f` / `.5f`, ou avec des variables double : `0 conditions are right`,
  le comportement attendu naïvement.
- La Figure 1 du sujet est vérifiée : 0.1f = `0x3DCCCCCD`, 0.1 =
  `0x3FB999999999999A`, avec les mêmes bits de signe / exposant / mantisse.
- **Le résultat dépend de la valeur.** 0.7f est arrondi vers le bas, mais 0.1f
  (0.100000001490…) et 0.6f (0.600000023841…) sont arrondis **vers le haut** :
  avec 0.6 au lieu de 0.7, `a < 0.6` serait faux. La règle n'est pas « un float
  est plus petit que le double », mais « le float et le double sont deux
  arrondis différents de la constante décimale ».
- **Combien de chiffres afficher ?** `%0.18f` affiche 0.699999999999999956, qui
  n'est toujours pas la valeur stockée (elle a 52 décimales). 17 chiffres
  significatifs identifient un double de façon unique (`%.17g`,
  `DBL_DECIMAL_DIG`), 9 un float (`%.9g`), et `%a` affiche la valeur exacte en
  hexadécimal (`0x1.6666666666666p-1`).

### 2. Don't test for equality, Relative errors – `equality.cpp`

- Le code du sujet affiche `-1.77636e-15` (précision par défaut de `cout` :
  6 chiffres). L'écart exact vaut −2⁻⁴⁹ = −1.7763568394002505·10⁻¹⁵ :
  `sqrt(10)²` = 10.000000000000002 est **le double qui suit 10** (1 ulp). La
  valeur du sujet, `-1.778636e-015`, est probablement une coquille pour
  `-1.776357e-015`.
- `FLT_EPSILON` = 2⁻²³ ≈ 1.19·10⁻⁷, `DBL_EPSILON` = 2⁻⁵² ≈ 2.22·10⁻¹⁶.

| Paire (attendu, calculé) | x − y | test absolu (< 10⁻¹²) | erreur relative | test relatif (≤ 10⁻¹²) |
|---|---|---|---|---|
| 10, sqrt(10)² | −1.78·10⁻¹⁵ | égaux | 1.8·10⁻¹⁶ | égaux |
| 2·10¹⁰, sqrt(2·10¹⁰)² | −3.81·10⁻⁶ | **différents** (faux) | 1.9·10⁻¹⁶ | égaux |
| 0.3, 0.1 + 0.2 | −5.55·10⁻¹⁷ | égaux | 1.9·10⁻¹⁶ | égaux |
| 10⁻¹⁵, 2·10⁻¹⁵ | −1.00·10⁻¹⁵ | **égaux** (faux : facteur 2) | 1.0 | différents |
| 0, sin(π) | −1.22·10⁻¹⁶ | égaux | inf | **différents** (faux) |

- Une tolérance absolue ne convient qu'aux nombres proches de 1 (trop petite
  pour 10¹⁰, trop grande pour 10⁻¹⁵). Une tolérance relative convient à toutes
  les grandeurs, sauf près de 0 où l'erreur relative n'a plus de sens. Une
  comparaison robuste combine les deux :
  `|a − b| ≤ max(rel × max(|a|, |b|), abs)` (`almostEqual`, la règle de
  `math.isclose` en Python). Les valeurs de `rel` et `abs` dépendent toujours
  de la précision attendue du calcul, comme le dit le sujet.
- Erreur relative du sujet : résultat 99.5, attendu 100 → 0.005.

### 3. Worry about addition (more than multiplication) – `sumMul.c`

Code du sujet :
`sum = 1.000000119209290, mul = 1.000000000000000, mul2 = 1.000000000000000`.

- Les 10 additions de 0.1f font 10 arrondis qui s'accumulent : 1 + 2⁻²³
  (1 ulp au-dessus de 1). Le produit `f * 10` n'est arrondi **qu'une fois** :
  1.0000000149 arrondi en float donne 1.
- `mul2` (`f * 10` passé à `printf`) est lui aussi calculé en float sur x86-64
  (10 est converti en float, `FLT_EVAL_METHOD = 0`). Avec `f * 10.0`, le
  produit est un double : 1.000000014901161, la valeur exacte de 0.1f × 10.

| n | somme float | erreur rel. | produit float n × 0.1f | somme de Kahan en float |
|---|---|---|---|---|
| 10 | 1.0000001 | 1.0·10⁻⁷ | 1.0 | 1.0 |
| 10³ | 99.999 | 9.6·10⁻⁶ | 100.0 | 100.0 |
| 10⁵ | 9 998.557 | 1.4·10⁻⁴ | 10 000.0 | 10 000.0 |
| 10⁶ | 100 958.34 | 9.6·10⁻³ | 100 000.0 | 100 000.0 |
| 10⁷ | 1 087 937 | **8.8·10⁻²** | 1 000 000.0 | 1 000 000.0 |
| 10⁸ | 2 097 152 | **7.9·10⁻¹** | 10 000 000.0 | 10 000 000.0 |

- À 10⁸ additions, la somme float est **bloquée à 2²¹ = 2 097 152** : à cet
  endroit, deux floats consécutifs sont espacés de 0.25, 0.1 est inférieur à la
  moitié de cet écart et `s + 0.1` est arrondi à `s` (absorption).
- Le produit et la somme de Kahan (partie 2) donnent la valeur correctement
  arrondie pour tout n (l'erreur relative de 1.5·10⁻⁸ est seulement l'arrondi
  en float de la valeur exacte n × 0.1f).
- En double, la même somme de 0.1 a une erreur relative de 1.9·10⁻⁹ pour
  n = 10⁸, qui croît **linéairement** avec n : le terme est toujours le même,
  donc les erreurs d'arrondi ont le même signe et ne se compensent pas.

### 4. Worry about subtraction : dérivée de sin – `derivative.cpp`

Le programme reproduit la sortie donnée dans le sujet chiffre pour chiffre
(colonne Visual C++), y compris les quatre `0` de la fin :

| h | (sin(1+h) − sin(1))/h | erreur | différence centrée | erreur | pas complexe | erreur |
|---|---|---|---|---|---|---|
| 10⁻¹ | `0.4..............` | 4.3·10⁻² | `0.53.............` | 9.0·10⁻⁴ | `0.541............` | 9.0·10⁻⁴ |
| 10⁻⁵ | `0.5402...........` | 4.2·10⁻⁶ | `0.54030230585....` | **1.1·10⁻¹¹** | `0.54030230587....` | 9.0·10⁻¹² |
| 10⁻⁸ | `0.540302302......` | **3.0·10⁻⁹** | `0.540302308......` | 2.6·10⁻⁹ | `0.540302305868140` | **0** |
| 10⁻¹² | `0.54034..........` | 4.3·10⁻⁵ | `0.5402...........` | 1.2·10⁻⁵ | `0.540302305868140` | 0 |
| 10⁻¹⁵ | `0.55.............` | 1.5·10⁻² | `0.55.............` | 1.5·10⁻² | `0.540302305868140` | 0 |
| 10⁻¹⁶ … 10⁻¹⁹ | `0` | 0.54 | `0` (dès 10⁻¹⁷) | 0.54 | `0.540302305868140` | 0 |

- Erreur ≈ h/2·|sin''| (troncature) + 2ε·|sin|/h (annulation dans la
  soustraction) : minimum pour h ≈ √ε ≈ 10⁻⁸, avec une erreur d'environ 10⁻⁸,
  comme observé. Pour h ≤ 10⁻¹⁶, 1 + h == 1 en double et le résultat vaut
  exactement 0.
- **Astuces qui changent le problème**, comme le conseille le sujet : la
  différence centrée (f(x+h) − f(x−h))/2h a une erreur de troncature en h², le
  meilleur h est donc ≈ ε^(1/3) ≈ 10⁻⁵ (erreur 1.1·10⁻¹¹). Le **pas complexe**
  f'(x) ≈ Im f(x + ih)/h ne fait **aucune soustraction** : pas d'annulation, et
  cos(1) est obtenu au dernier chiffre près pour tout h ≤ 10⁻⁸
  (`std::complex`).
- **Pourquoi « gcc sous Linux » donnait −0.017…, −0.17…, −1.7… et 17….** Avec
  gcc 16 en x86-64, on obtient 0, comme Visual C++. Les anciennes valeurs
  viennent du **FPU x87** utilisé par les compilateurs x86 32 bits
  (`FLT_EVAL_METHOD = 2`) : les résultats intermédiaires restent dans des
  registres de 80 bits (mantisse de 64 bits). Si sin(1+h) est arrondi en
  double (stocké en mémoire) alors que sin(1) reste dans un registre de 80
  bits, le numérateur vaut l'erreur d'arrondi de sin(1) en double,
  −1.79·10⁻¹⁸, divisée par h. Reproduit avec `long double` : −0.0179, −0.179,
  −1.79, −17.9 (la 4ᵉ valeur est négative elle aussi). Le même source donne des
  résultats différents selon le compilateur, ses options et le matériel.

### 5. Floating point numbers have finite ranges – `ranges.cpp`

- `cout << f << " " << f+1` affiche `1.67772e+07 1.67772e+07` : « 16777216
  deux fois » demande au moins `setprecision(8)`. Ensuite `f + 1` = 16777216 :
  `1` est converti en float et 2²⁴ + 1 n'est pas un float (mantisse de 24
  bits). C'est un **cas d'égalité** entre 2²⁴ et 2²⁴ + 2, arrondi à la mantisse
  paire. `f + 1.0` est calculé en double : 16777217.
- Pour x = 2⁵³ : `(x + 1) − x = 0` (même cas d'égalité en double),
  `(x + 2) − x = 2`, `(x + 3) − x = 4` (égalité arrondie à 2⁵³ + 4).
- `Sorry!` est affiché. Plus précisément, `DBL_EPSILON` = 2⁻⁵² est **l'écart
  entre 1 et le double suivant**. 1 + ε/2 est exactement à mi-chemin entre 1 et
  1 + ε et il est arrondi au pair, c'est-à-dire à 1. Mais 1 + 0.51·ε ≠ 1 : le
  plus petit e tel que 1 + e ≠ 1 est juste au-dessus de ε/2, pas ε. Et
  1 − ε/2 ≠ 1, car l'écart sous 1 vaut ε/2.

| Valeur | écart au float suivant | écart au double suivant |
|---|---|---|
| 1 | 1.19·10⁻⁷ | 2.22·10⁻¹⁶ |
| 10⁶ | 0.0625 | 1.16·10⁻¹⁰ |
| 2²⁴ = 16 777 216 | 2 | 3.73·10⁻⁹ |
| 10⁹ | 64 | 1.19·10⁻⁷ |
| 2⁵³ ≈ 9.0·10¹⁵ | 1.07·10⁹ | 2 |
| 10³⁰ | 7.6·10²² | 1.4·10¹⁴ |

La précision **absolue** diminue quand la grandeur augmente (la précision
relative est constante) : un petit nombre ajouté à un grand est perdu
(`sumMul.c`, série harmonique de `associativity.c`, sommes float de la partie
2).

| | float | double | long double (x87) |
|---|---|---|---|
| Bits de la mantisse | 24 | 53 | 64 |
| Chiffres toujours conservés / nécessaires pour un aller-retour | 6 / 9 | 15 / 17 | 18 / 21 |
| EPSILON | 1.19·10⁻⁷ | 2.22·10⁻¹⁶ | 1.08·10⁻¹⁹ |
| Plus petit normalisé / sous-normal | 1.18·10⁻³⁸ / 1.40·10⁻⁴⁵ | 2.23·10⁻³⁰⁸ / 4.94·10⁻³²⁴ | 3.36·10⁻⁴⁹³² / 3.65·10⁻⁴⁹⁵¹ |
| Plus grand | 3.40·10³⁸ | 1.80·10³⁰⁸ | 1.19·10⁴⁹³² |
| Entiers consécutifs jusqu'à | 2²⁴ ≈ 1.7·10⁷ | 2⁵³ ≈ 9.0·10¹⁵ | 2⁶⁴ ≈ 1.8·10¹⁹ |

### 6. Use logarithms to avoid overflow and underflow – `logFactorial.cpp`

170! = 7.26·10³⁰⁶ mais 171! = inf (`DBL_MAX` = 1.80·10³⁰⁸). Calculé
directement, 200!/(190! 10!) = inf/inf = **NaN** (affiché `-nan`).

| Méthode | Résultat | Erreur relative |
|---|---|---|
| Exacte (entiers 64 bits) | 22 451 004 309 013 280 | – |
| `exp(logFactorial(200) − logFactorial(190) − logFactorial(10))` | 22 451 004 309 011 012 | 1.0·10⁻¹³ |
| `exp(lgamma(201) − lgamma(191) − lgamma(11))` | 22 451 004 309 011 172 | 9.4·10⁻¹⁴ |
| Produit des rapports (190 + i)/i, i = 1..10 | 22 451 004 309 013 280 | 0 |

- Les logarithmes évitent le dépassement mais coûtent 3 chiffres :
  log 200! = 863.23 et log 190! = 810.48 ont une précision absolue d'environ
  10⁻¹³ (l'ulp de 863), la soustraction donne 37.65 avec la même erreur
  absolue, et cette erreur devient l'erreur **relative** de exp(37.65). Même
  `lgamma` ne peut pas l'éviter. Quand c'est possible, une formule sans valeurs
  intermédiaires énormes (produit des rapports) est meilleure.
- Dépassement **et** sous-dépassement dans la même formule : P(X = 1000) pour
  X ~ B(2000, ½) vaut C(2000, 1000) × 0.5²⁰⁰⁰, avec C(2000, 1000) ≈ 2·10⁶⁰⁰ =
  inf et 0.5²⁰⁰⁰ ≈ 8.7·10⁻⁶⁰³ = 0, donc inf × 0 = NaN. Avec `lgamma` :
  0.017839011145854369, erreur relative 2.7·10⁻¹⁵ (valeur exacte
  0.01783901114585432, calculée avec les fractions de Python).

### 7. Numeric operations don't always return numbers – `notNumbers.cpp`

`2 * DBL_MAX` affiche `inf`, `sqrt(-2)` affiche `-nan` (la glibc affiche le
signe : le NaN par défaut du x86 a son bit de signe à 1).

| x | `x - x == 0` | `x == x` | `std::isnan(x)` |
|---|---|---|---|
| 1 | vrai | vrai | faux |
| inf | **faux** (inf − inf = NaN) | vrai | faux |
| NaN | faux | **faux** | vrai |

`x - x == 0` teste donc « x est fini », et `x == x` teste « x n'est pas un
NaN » (infinis compris), comme l'indique le sujet.

Drapeaux d'exception IEEE 754 (`<cfenv>`), levés à l'exécution :

| Opération | Résultat | Drapeaux |
|---|---|---|
| 2 × DBL_MAX, exp(1000) | inf | OVERFLOW INEXACT |
| 1 / 0, log(0) | inf, −inf | DIVBYZERO |
| 0 / 0, sqrt(−2), log(−1), inf − inf, 0 × inf | NaN | INVALID |
| inf + 1 | inf | – |
| exp(−1000), DBL_TRUE_MIN / 2 | 0 | UNDERFLOW INEXACT |
| DBL_MIN / 2 | 1.11·10⁻³⁰⁸ (sous-normal) | – (résultat exact : pas de drapeau underflow) |
| DBL_MIN / 3 | 7.42·10⁻³⁰⁹ (sous-normal) | UNDERFLOW INEXACT |
| 0.1 + 0.2 | 0.30000000000000004 | INEXACT |
| 0.5 + 0.25 | 0.75 | – |

- Toute comparaison avec un NaN est fausse (`NaN == NaN`, `NaN < 1`,
  `NaN >= 1`) sauf `!=` : `!(x < 1)` n'est donc pas équivalent à `x >= 1`.
- **Un NaN ne se propage pas toujours** : `pow(NaN, 0) = 1`,
  `fmax(NaN, 1) = 1` (Goldberg §3.2.2 : un maximum doit ignorer les NaN),
  `hypot(inf, NaN) = inf`, et `std::max(NaN, 1.0) = NaN` mais
  `std::max(1.0, NaN) = 1` : le résultat dépend de l'ordre des arguments.
- **Avec `-ffast-math` (`notNumbers_fast`)** : `x == x` est vrai pour NaN,
  `std::isnan(NaN)` est faux, `x - x == 0` est vrai pour inf et NaN,
  `NaN * 0 = 0`, et les sous-normaux sont remplacés par 0 (DBL_MIN / 2 = 0,
  mode flush-to-zero activé au lancement du programme). Le compilateur suppose
  qu'il n'y a ni NaN ni infini : les tests de cette section sont supprimés.

### 8. Out of order execution : les opérations ne sont pas associatives – `associativity.c`

| Expression | Résultat |
|---|---|
| (10⁻³ + 1) − 1 | 0.000999999999999890 |
| 10⁻³ + (1 − 1) | 0.001000000000000000 |
| (2⁶⁰ − 2⁶⁰) + 1 | 1 |
| (2⁶⁰ + 1) − 2⁶⁰ | 0 |
| (0.1 + 0.2) + 0.3 | 0.6000000000000001 |
| 0.1 + (0.2 + 0.3) | 0.6 |
| (0.1 × 0.2) × 0.3 | 0.006000000000000001 |
| 0.1 × (0.2 × 0.3) | 0.006 |
| (10³⁰⁸ × 10) × 0.1 | inf |
| 10³⁰⁸ × (10 × 0.1) | 10³⁰⁸ |

- Les produits ne sont pas associatifs non plus, mais l'écart reste de 1 ulp,
  sauf si un résultat intermédiaire déborde.
- **FMA** : `fma(0.1, 10, −1)` = 5.55·10⁻¹⁷ (un seul arrondi : c'est l'erreur
  exacte du produit 0.1 × 10), alors que `p = 0.1 * 10; p − 1` = 0 (deux
  arrondis). Compilé avec `-mfma` (`associativity_fma`) ou `-march=native`,
  **gcc fusionne de lui-même `p = a * b; p - c` en un FMA**
  (`-ffp-contract=fast` est le défaut en mode GNU C) et le même source affiche
  5.55·10⁻¹⁷. `-ffp-contract=off` ou `-std=c17` interdisent cette contraction.
  Le FMA est plus précis, mais les résultats changent avec les options de
  compilation et le processeur cible.
- **Ordre d'une longue somme** : série harmonique H(10⁷) en float. Dans
  l'ordre croissant des indices (1 → n), la somme **cesse d'augmenter à
  i = 2 097 152 = 2²¹**, quand 1/i devient plus petit que la moitié de l'écart
  entre floats autour de 15.4 : 15.4036827 (erreur 7.7 %). Dans l'ordre inverse
  (petits termes d'abord) : 16.6860313 (5.6·10⁻⁴). Exact :
  16.695311365859855.
- **À propos de l'exécution dans le désordre dans le processeur.** Le moteur
  d'exécution dans le désordre d'un cœur exécute les instructions
  indépendantes dans un autre ordre, mais chaque instruction reçoit les mêmes
  opérandes et les résultats sont validés dans l'ordre du programme : sur un
  seul thread, un binaire donné calcule des résultats identiques bit à bit à
  chaque exécution (vérifié : `sumBench` affiche les mêmes 17 chiffres à
  chaque exécution). Il accélère même certaines sommes sans les modifier
  (`pairwise` est plus rapide que `naive` parce que ses blocs de 128 valeurs
  sont des chaînes d'additions indépendantes). Ce qui change réellement l'ordre
  des opérations flottantes, c'est (1) le **compilateur** (vectorisation
  autorisée par `-ffast-math`, contraction en FMA), (2) le **support
  d'exécution parallèle** (découpage de la boucle selon le nombre de threads,
  ordre de combinaison des résultats partiels selon la vitesse des threads) et
  (3) les bibliothèques qui choisissent leur code à l'exécution (noyaux BLAS
  selon le processeur, par exemple). Cette non-reproductibilité dynamique est
  mesurée avec `reproSum` (partie 2).

### 9. Pour aller plus loin – étudiants en IA : softmax – `softmax.c`

En float (32 bits), comme dans les frameworks d'apprentissage profond :

| Logits | Softmax naïf | Softmax stable |
|---|---|---|
| (1, 2, 3, 4) | 0.032 0.087 0.237 0.644 | identique |
| (1000, 1001, 1002, 1003) | NaN NaN NaN NaN (inf / inf) | 0.032 0.087 0.237 0.644 |
| (−1000, −1001, −1002, −1003) | NaN NaN NaN NaN (0 / 0) | 0.644 0.237 0.087 0.032 |
| (100, 50, 0, −50) | NaN 0 0 0 (exp(100) > FLT_MAX ; correct en double) | 1 1.9·10⁻²² 3.8·10⁻⁴⁴ 0 |

- exp déborde pour z > 88.7 en float (709.8 en double, seulement 11.1 en demi-
  précision, utilisée pour l'entraînement en précision mixte) et vaut 0 pour
  des logits très négatifs.
- Le softmax ne change pas quand on soustrait la même constante à tous les
  logits. Avec m = max zⱼ, chaque exp(zᵢ − m) est ≤ 1 (pas de dépassement) et
  le dénominateur est ≥ 1 (pas de 0/0) : les logits (1000 … 1003) donnent
  exactement le résultat de (1 … 4).
- Entropie croisée −log pᵢ pour les logits (0, 0, 0, −200) et la dernière
  classe : en naïf p₃ = 0 (sous-dépassement), donc la perte vaut **inf** ; avec
  le log-softmax zᵢ − m − log Σⱼ exp(zⱼ − m) (astuce log-sum-exp), la perte vaut
  **201.099** (exact : 200 + ln 3 = 201.0986).

### 10. Pour aller plus loin – HPC : sommes compensées

#### Q1 – Le professeur William Kahan

William Morton Kahan (né en 1933 à Toronto), professeur émérite de
mathématiques et d'informatique à l'Université de Californie à Berkeley, est
**le principal architecte de la norme IEEE 754 (1985)**, d'où son surnom de
« père de la virgule flottante ». Il a reçu le **prix Turing en 1989** pour
ses contributions fondamentales à l'analyse numérique. À partir de 1976, il a
travaillé avec Intel sur le coprocesseur flottant 8087, dont l'arithmétique a
servi de base à la norme. Il est aussi l'auteur du programme « paranoia »
(1983), qui teste l'arithmétique flottante d'une machine, et de formules
numériquement robustes (aire des triangles très aplatis, racines d'une
équation du second degré, …).

Ce qu'il a proposé pour que les erreurs d'arrondi de la norme puissent être
compensées, ou au moins détectées :

- **La sommation compensée (1965)** : W. Kahan, « Further remarks on reducing
  truncation errors », Communications of the ACM 8(1). L'erreur de chaque
  addition est calculée exactement et réinjectée dans le terme suivant
  (Théorème 8 de Goldberg).
- **La précision étendue** dans la norme (le format 80 bits du x87, avec une
  mantisse de 64 bits), pour calculer les résultats intermédiaires avec plus de
  précision que les données stockées.
- **L'arrondi correct** de +, −, ×, / et √, qui rend exactes les
  transformations sans erreur comme TwoSum ; le **sous-dépassement progressif**
  (nombres sous-normaux : x − y = 0 ⇔ x = y) ; **inf, NaN et les drapeaux
  d'exception**, pour qu'un calcul puisse continuer et que le programme puisse
  tester ce qui s'est passé ; les **modes d'arrondi dirigés**, utiles aussi pour
  estimer la sensibilité d'un résultat aux erreurs d'arrondi en relançant le
  calcul avec un autre mode d'arrondi.

#### Q2 – Sommation de Kahan contre sommation classique : résultats et temps de calcul – `sumBench.c`

Algorithme de Kahan (Goldberg, Théorème 8), dans `summation.c` :

```c
s = 0.0; c = 0.0;
for (i = 0; i < n; i++)
{
   y = x[i] - c;        /* terme corrigé des bits perdus avant          */
   t = s + y;           /* les bits de poids faible de y sont perdus ici */
   c = (t - s) - y;     /* = -(bits perdus), calculé exactement         */
   s = t;
}
```

Bornes d'erreur (ε = 2⁻⁵³) : naïve ≤ (n − 1)·ε·Σ|xᵢ|, Kahan
≤ (2ε + O(nε²))·Σ|xᵢ|. L'erreur de la somme de Kahan ne croît plus avec n.

**Protocole.** 4 jeux de N = 10⁷ doubles (80 Mo, plus que les 16 Mo du cache
L3) tirés avec MT19937 (`genrand_res53`, 53 bits aléatoires par valeur), un
flux par jeu de données : U[0,1] ; 10^U[−10,10] (20 ordres de grandeur) ;
U[−1,1] (annulations) ; ±10^U[−10,10]. La référence est la somme exacte arrondie
une seule fois (`exactSum.c`, vérifiée bit à bit sur les 12 cas de test de
`math.fsum` de CPython). Temps = meilleur de 5 exécutions. Les valeurs
calculées sont identiques à chaque exécution.

Erreur relative, double, N = 10⁷ (0 = la somme correctement arrondie, bit à
bit) :

| Algorithme | U[0,1] | 10^U[−10,10] | U[−1,1] (cond 5.8·10³) | ±10^U[−10,10] (cond 898) | ns / valeur |
|---|---|---|---|---|---|
| naïf | 1.1·10⁻¹⁴ | **5.8·10⁻¹²** | 6.5·10⁻¹⁴ | 9.8·10⁻¹⁴ | 0.73 |
| naïf, 8 accumulateurs | 5.4·10⁻¹⁵ | 7.1·10⁻¹³ | 3.5·10⁻¹⁴ | 9.3·10⁻¹⁵ | 0.40 |
| **Kahan** | **0** | **0** | **0** | **0** | **2.87** |
| Kahan, 8 accumulateurs | 0 | 0 | 0 | 0 | 0.67 |
| exact (fsum) | 0 | 0 | 0 | 0 | 5.3 – 16.4 |

Avec N = 10⁸ (`results/sumBench_1e8.txt`) : naïf 4.6·10⁻¹³ (U[0,1]) et
6.2·10⁻¹¹ (10^U[−10,10]), Kahan toujours exact, mêmes temps par valeur.

Simple précision (les mêmes valeurs arrondies en float), N = 10⁷ :

| Algorithme | U[0,1] | 10^U[−10,10] | U[−1,1] | ±10^U[−10,10] | ns / valeur |
|---|---|---|---|---|---|
| naïf (accumulateur float) | 2.7·10⁻⁵ | **3.0·10⁻³** | 2.5·10⁻⁵ | 7.8·10⁻⁵ | 0.72 |
| Kahan (accumulateur float) | 3.6·10⁻⁸ | 2.8·10⁻⁸ | 2.2·10⁻⁹ | 3.2·10⁻⁹ | 2.87 |
| naïf dans un accumulateur double | 5.4·10⁻¹⁵ | 5.8·10⁻¹² | 2.6·10⁻¹⁶ | 5.2·10⁻¹⁴ | 0.71 |
| naïf float, **N = 10⁸** | **16 777 216 au lieu de 5.0·10⁷ (−66 %)** | 3.3·10⁻² | 1.3·10⁻⁴ | 2.8·10⁻⁶ | 0.72 |

- **Résultats.** Sur tous les jeux de données, Kahan donne **la somme
  correctement arrondie** : erreur nulle en double ; en float, le float le
  plus proche de la somme exacte (3.6·10⁻⁸ n'est que l'arrondi de
  5 002 029.678 en 5 002 029.5). La somme naïve perd 2 à 4 chiffres en double,
  et jusqu'au résultat entier en float : avec N = 10⁸ valeurs U[0,1], la somme
  float s'arrête à 2²⁴ = 16 777 216, où deux floats sont espacés de 2 et où
  tout terme inférieur à 1 est absorbé (section « finite ranges » du sujet, au
  cœur d'une vraie réduction).
- **La plage des valeurs compte plus que N.** Avec 20 ordres de grandeur, la
  moitié des valeurs est inférieure à 1 et elles sont en partie absorbées par
  une somme courante d'environ 10¹⁵. Les bits perdus ont tous le même signe,
  donc l'erreur est systématique : 500 fois plus grande qu'avec U[0,1]
  (5.8·10⁻¹² contre 1.1·10⁻¹⁴). La somme calculée est inférieure à la somme
  exacte.
- **Temps de calcul.** Kahan est **3.9 fois plus lent** que la boucle naïve
  (2.87 contre 0.73 ns par valeur), alors qu'il ne fait que 4 additions au lieu
  d'une : la boucle est limitée par la **latence** d'une addition flottante
  (3 cycles sur Zen 3, soit ≈ 0.72 ns à 4.2 GHz). La boucle naïve est une
  chaîne d'additions dépendantes, une par valeur ; dans la boucle de Kahan, la
  compensation enchaîne 4 opérations dépendantes par valeur
  (y → t → t − s → c) : 4 × 0.72 = 2.88 ns, exactement ce qui est mesuré.
- **Casser la chaîne, garder la compensation.** Avec 8 accumulateurs
  indépendants (l'accumulateur j somme les valeurs i ≡ j mod 8), **Kahan prend
  0.67 ns par valeur : plus rapide que la boucle naïve séquentielle**, avec la
  même précision. La somme naïve à 8 accumulateurs prend 0.40 ns (limitée par
  la bande passante mémoire d'un cœur). En parallèle, les threads jouent le
  même rôle (voir `reproSum` : même temps pour les sommes naïve et de Kahan
  avec 12 threads).
- Pour des données float, la correction la moins chère est un **accumulateur
  double** : gratuit (0.71 ns) et précis à 10⁻¹⁴–10⁻¹⁶ près (la « méthode
  simple » de Goldberg §3.2.3). Kahan en float donne le float le plus proche
  de la somme exacte.

#### Le piège de `-ffast-math` – `sumBench_fast`

`summation.c` compilé en `-O3 -ffast-math` (ce fichier seulement ; la
référence `exactSum.c` refuse de compiler avec cette option) :

| Algorithme | U[0,1], `-O2` | U[0,1], `-O3 -ffast-math` | ns / valeur, `-O2` → fast-math |
|---|---|---|---|
| naïf | 1.1·10⁻¹⁴ | 3.4·10⁻¹⁵ | 0.73 → 0.38 |
| Kahan | 0 | **3.4·10⁻¹⁵ (= naïf)** | 2.87 → 0.38 |
| Neumaier, Sum2, Klein | 0 | **3.4·10⁻¹⁵ (= naïf)** | 0.95 – 1.36 → 0.38 |
| Kahan float, 10^U[−10,10] | 2.8·10⁻⁸ | **7.6·10⁻⁴ (= naïf float)** | 2.86 → 0.19 |

- Avec `-ffast-math`, le compilateur a le droit d'utiliser l'algèbre des
  réels : c = ((s + y) − s) − y = 0, la compensation disparaît et **toutes les
  sommes compensées deviennent la somme naïve, bit à bit** (Kahan à 8
  accumulateurs devient la somme naïve à 8 accumulateurs ; il en va de même sur
  chaque ligne de l'étude en fonction du conditionnement). Goldberg prévenait
  déjà en 1991 (§3.2.3) : « an optimizer that believed floating-point
  arithmetic obeyed the laws of algebra would conclude that C = 0, rendering
  the algorithm completely useless ».
- La somme naïve est vectorisée (plusieurs sommes partielles en parallèle) :
  elle est plus rapide et ici même un peu plus précise, mais c'est **un autre
  résultat** (un autre ordre) ; la petite somme 1 + 10¹⁰⁰ + 1 − 10¹⁰⁰ donne même
  2 au lieu de 0. `-Ofast` inclut `-ffast-math`, et les compilateurs Intel
  utilisent par défaut un modèle flottant « fast » : il faut y compiler les
  sommes compensées avec `-fp-model=precise`.

#### Q3 – Autres variantes de sommes compensées

| Algorithme | Auteur(s), année | Idée | Erreur (S = Σ\|xᵢ\|) | Temps ici |
|---|---|---|---|---|
| naïf | – | s += x | ≤ (n − 1)·ε·S | 0.73 ns |
| par paires (cascade) | utilisé par NumPy, Julia | somme des deux moitiés, récursivement (naïf en dessous de 128 valeurs) | ≈ ε·log₂n·S | 0.44 ns |
| **Kahan** | Kahan, 1965 | les bits perdus sont réinjectés dans le terme suivant | ≤ (2ε + O(nε²))·S | 2.87 ns |
| **Kahan–Babuška–Neumaier** | Babuška 1969, Neumaier 1974 | l'erreur de s + x est calculée à partir du plus grand opérande, les erreurs sont sommées à part | ≤ ε·\|Σx\| + O(n²ε²)·S | 0.95 ns |
| **Sum2 (TwoSum)** | Møller 1965, Knuth 1969 (TwoSum) ; Ogita, Rump, Oishi 2005 | erreur exacte de chaque addition sans aucun test (6 opérations) | comme la somme naïve calculée en précision double de la précision de travail : ≤ ε·\|Σx\| + γ²ₙ₋₁·S | 1.19 ns |
| **Klein** | Klein, 2006 | Kahan–Babuška d'ordre 2 : les erreurs de la compensation sont elles-mêmes compensées | observé ~10⁴ fois plus robuste que Neumaier | 1.35 ns |
| SumK, AccSum | Ogita, Rump, Oishi 2005–2008 | K − 1 passes sans erreur ; arrondi fidèle | comme K fois la précision | – |
| doublement compensée | Priest, 1992 | tri par \|x\| décroissant, puis double compensation | ≤ 2ε·\|Σx\| | tri O(n log n) |
| accumulateur plus large | – | long double (x87), double pour des données float, double-double, `__float128` | dépend du format | 1.57 ns (long double) |
| **exacte (fsum)** | Shewchuk 1997, `math.fsum` de Python | liste exacte de termes partiels sans recouvrement, arrondie une seule fois | correctement arrondie | 5 – 16 ns |
| accumulateur long, sommes reproductibles | Kulisch ; ReproBLAS (Demmel, Nguyen 2013) ; ExBLAS (Collange, Defour, Graillat, Iakymchuk 2015) | accumulateur en virgule fixe d'environ 4200 bits, ou pré-arrondi de toutes les valeurs sur une grille commune | exacte, ou identique bit à bit quel que soit l'ordre | – |

Les algorithmes en gras, le long double, la somme par paires et la somme exacte
sont implémentés dans `summation.c` / `exactSum.c` et mesurés ci-dessus.

**Petit exemple (le cas de Neumaier)** : 1 + 10¹⁰⁰ + 1 − 10¹⁰⁰ = 2. Le naïf, la
somme par paires, Kahan et le long double donnent **0** ; Neumaier, Sum2, Klein
et la somme exacte donnent **2**. Kahan échoue parce que le nouveau terme
(10¹⁰⁰) est plus grand que la somme courante : (t − s) − y n'est alors plus
l'erreur exacte.

**Sommes mal conditionnées.** N = 10⁶ valeurs : des paires (v, −v) avec |v|
jusqu'à 2^E, qui s'annulent exactement, plus de petites valeurs U[0,1) ; la
somme exacte est connue. cond = Σ|xᵢ| / |Σxᵢ| ; erreur relative (≥ 1 : plus
aucun chiffre correct) :

| cond | naïf | par paires | Kahan | Neumaier = Sum2 | Klein | long double | exact |
|---|---|---|---|---|---|---|---|
| 4 | 3.4·10⁻¹⁵ | 1.2·10⁻¹⁶ | 0 | 0 | 0 | 1.2·10⁻¹⁶ | 0 |
| 3.0·10⁵ | 5.5·10⁻¹² | 3.9·10⁻¹³ | 1.4·10⁻¹⁴ | 0 | 0 | 6.7·10⁻¹⁵ | 0 |
| 1.6·10¹¹ | 6.3·10⁻⁶ | 1.9·10⁻⁷ | 3.9·10⁻⁸ | 0 | 0 | 3.5·10⁻⁹ | 0 |
| 1.3·10¹⁴ | **0.90** | 1.7·10⁻⁴ | 1.9·10⁻⁵ | 2.9·10⁻¹⁴ | 0 | 5.1·10⁻⁶ | 0 |
| 1.1·10¹⁷ | 9.3 | 1.1 | **8.0·10⁻²** | 3.6·10⁻¹³ | 2.3·10⁻¹⁶ | 0.59 | 0 |
| 1.0·10²⁰ | 7.8·10³ | 2.7·10² | 73 | 1.1·10⁻¹⁰ | 3.2·10⁻¹³ | 7.4 | 0 |
| 8.3·10²⁵ | 1.3·10⁹ | 1.4·10⁸ | 4.3·10⁷ | 5.2·10⁻⁵ | 1.7·10⁻⁸ | 1.9·10⁶ | 0 |
| 7.3·10²⁸ | 1.0·10¹² | 3.6·10¹⁰ | 5.0·10¹⁰ | **0.98** | 4.3·10⁻⁶ | 6.8·10⁸ | 0 |
| 7.0·10³¹ | 7.6·10¹⁵ | 2.2·10¹⁴ | 2.1·10¹³ | 35 | **4.9·10⁻²** | 4.1·10¹² | 0 |
| 4.6·10⁵⁸ | 6.9·10⁴¹ | 2.1·10⁴¹ | 3.7·10⁴⁰ | 9.1·10²⁷ | 1.0 | 6.0·10³⁸ | **0** |

- L'erreur relative d'une somme est gouvernée par **ε × cond**, bien plus que
  par N. Chaque algorithme a son point de rupture : la somme naïve et le long
  double vers 10¹⁴–10¹⁷ ; Kahan vers 10¹⁶ (sa borne est relative à Σ|xᵢ| :
  c'est une méthode du premier ordre) ; Neumaier et Sum2 vers 10²⁸ (deux fois
  la précision de travail : l'erreur croît comme ε²·cond) ; Klein vers 10³² ;
  seule la somme exacte fonctionne toujours.
- **Neumaier et Sum2 donnent les mêmes résultats bit à bit** (sur tous les jeux
  de données) : les deux calculent l'erreur exacte de chaque addition
  (Fast2Sum sur le plus grand opérande pour Neumaier, TwoSum pour Sum2) et
  additionnent les mêmes erreurs dans le même ordre. Sum2 remplace le test par
  2 additions (pas de branchement).
- Kahan suffit pour les sommes bien conditionnées (tous les jeux de données du
  premier tableau). Neumaier et Sum2 sont moins chers (0.95 et 1.19 ns contre
  2.87 ns : leur compensation est en dehors de la chaîne d'additions
  dépendantes) et plus robustes : pour une somme séquentielle, ce sont les
  meilleurs compromis.

#### Phase de réduction et reproductibilité – `reproSum.c`

Les mêmes 10⁷ valeurs U[−1,1] (somme exacte −859.10864839042449, 1 ulp =
1.1·10⁻¹³).

**1. Sommes séquentielles, les mêmes valeurs dans 6 ordres** (croissant,
décroissant, |x| croissant, |x| décroissant, 2 permutations aléatoires) : la
somme naïve donne **6 résultats différents** (de 2.3·10⁻¹³ avec |x| croissant
à 2.9·10⁻¹⁰ avec |x| décroissant) ; Kahan, Neumaier et la somme exacte donnent
**le même résultat correctement arrondi** dans tous les ordres.

**2. OpenMP, de 1 à 12 threads** (somme − exact, puis temps) :

| Threads | naïf `reduction(+)` | Kahan par thread, combinaison exacte | exact (termes partiels par thread, fusionnés) | temps naïf | temps Kahan | temps exact |
|---|---|---|---|---|---|---|
| 1 | 5.6·10⁻¹¹ | 0 | 0 | 7.3 ms | 28.6 ms | 53.4 ms |
| 2 | 7.5·10⁻¹¹ | 0 | 0 | 3.8 ms | 14.5 ms | 27.2 ms |
| 3 | −1.2·10⁻¹⁰ | 0 | 0 | 2.6 ms | 9.7 ms | 18.7 ms |
| 4 | 2.9·10⁻¹¹ | 0 | 0 | 2.0 ms | 7.3 ms | 14.1 ms |
| 6 | 3.0·10⁻¹² | 0 | 0 | 2.0 ms | 4.9 ms | 9.8 ms |
| 8 | 7.1·10⁻¹¹ | 0 | 0 | 2.1 ms | 3.7 ms | 7.6 ms |
| 12 | −2.7·10⁻¹¹ | 0 | 0 | 2.6 ms | **2.5 ms** | 5.3 ms |
| résultats distincts | **7** | 1 | 1 | | | |

**3. Même binaire, mêmes données, 12 threads, 10 exécutions** : la
`reduction(+)` naïve avec `schedule(static)` donne **2 à 5 résultats
différents** sur 10 exécutions (cela dépend de l'exécution) ; avec
`schedule(dynamic)`, **10 résultats différents sur 10** ; Kahan par thread et
la somme exacte donnent toujours le même résultat.

- Le résultat d'une réduction parallèle naïve dépend du nombre de threads
  (découpage de la boucle) **et de l'exécution** : à la fin de la boucle,
  chaque thread ajoute sa somme partielle à la variable partagée au moment où
  il arrive, donc l'ordre de ces dernières additions dépend de la vitesse des
  threads (et avec `dynamic`, chaque thread reçoit même des blocs différents).
  C'est la non-reproductibilité dynamique décrite dans le sujet ; elle vient de
  l'ordonnancement des threads.
- Compenser à l'intérieur de chaque thread et combiner exactement les résultats
  partiels, dans un ordre fixe (numéro de thread), rend le résultat
  indépendant du nombre de threads, ici bit à bit. La somme parallèle exacte
  (termes partiels exacts par thread, fusionnés exactement) est reproductible
  **par construction**, quels que soient le découpage et l'ordre.
- **Coût.** La réduction naïve est limitée par la bande passante mémoire dès 4
  threads (80 Mo en 2 ms ≈ 40 Go/s). Kahan, 4 fois plus lent sur 1 thread,
  rattrape son retard : **avec 12 threads, Kahan prend le même temps que la
  somme naïve (2.5 ms)**. En HPC, où les réductions sont limitées par la
  mémoire, la compensation est presque gratuite ; la somme exacte coûte deux
  fois plus (5.3 ms).

## Remarques

- **Codes du sujet** : les guillemets typographiques `“ ”` du PDF sont
  remplacés par `"`, les `#include`, déclarations et `\n` manquants sont
  ajoutés, `using namespace std` est remplacé par des `std::` explicites ; à
  part cela, les codes sont inchangés.
- **Coquilles probables du sujet** : `-1.778636e-015` (valeur exacte
  −1.7763568·10⁻¹⁵) ; « affiche 16777216 deux fois » demande au moins
  `setprecision(8)` (`cout` affiche `1.67772e+07` par défaut) ; la 4ᵉ valeur de
  « gcc » est −17.9, négative comme les autres.
- **Options de compilation** : tout est compilé en `-O2`, sans `-ffast-math`
  (qui casse les sommes compensées et les tests de NaN) et sans `-march=native`
  (contraction en FMA, binaires non portables). Toutes les valeurs de ce
  fichier sont donc celles de l'arithmétique IEEE 754 stricte, identiques sur
  n'importe quelle machine x86-64, sauf les temps et les résultats dépendant
  de l'exécution de `reproSum` (section 3).
- **Valeurs aléatoires sur 53 bits** : `genrand_res53` (MT19937 du TP3,
  découpé en `mt19937ar.h` / `mt19937ar.c`). Avec `genrand_real1` (32 bits
  aléatoires), les valeurs seraient des multiples de 2⁻³² : des sommes allant
  jusqu'à environ 2²¹ d'entre elles seraient exactes en double, ce qui
  masquerait les erreurs d'arrondi étudiées ici.
- **La référence** : `exactSum.c` (algorithme de Shewchuk, comme `math.fsum`
  de Python) donne la somme exacte arrondie une seule fois. Elle est vérifiée
  bit à bit sur les 12 cas de test de CPython (`test_fsum`) au lancement de
  `sumBench`, et refuse de compiler avec `-ffast-math` (`#error`).
- **Guide de style ISIMA** (fourni avec le TP) : chaque fichier commence par un
  bloc de commentaires (nom et rôle), chaque fonction est précédée d'un bloc
  (rôle, entrée, sortie), les types utilisateur se terminent par `_t`, les
  constantes sont en majuscules (`const` en C++), les déclarations sont
  alignées, l'indentation est de 3 espaces, les lignes font au plus 80
  caractères, et les fichiers `.h` ne contiennent que des déclarations
  (contrairement au TP3, où les générateurs étaient des fonctions
  `static inline` dans des en-têtes).
- **`long double`** : en x86-64, c'est le format x87 de 80 bits (mantisse de
  64 bits). Il est ici deux fois plus lent que le double et n'est pas portable
  (c'est un simple double sous Windows/MSVC et sur les Mac ARM, un format 128
  bits logiciel sous Linux ARM64).

## Dépannage

| Problème | Solution |
|---|---|
| `make: command not found` | Installer `make`, ou utiliser les commandes de la section Compilation |
| `omp.h: No such file` / `undefined reference to 'omp_get_wtime'` | Compiler `reproSum` avec `-fopenmp` (installer `libgomp` si besoin) |
| `warning: #warning "summation.c built with -ffast-math…"` | Normal : c'est `summation_fast.o`, compilé exprès pour `sumBench_fast` |
| `#error "exactSum.c must not be built with -ffast-math"` | `-ffast-math` ou `-Ofast` a été ajouté à `CFLAGS` : le retirer |
| `Permission non accordée` / `Permission denied` sur `RunAll` | `chmod +x RunAll` |
| `malloc: Cannot allocate memory` | N trop grand pour la mémoire : `./sumBench 1e7` (80 Mo par tableau) |
| `Illegal instruction` avec `associativity_fma` | Processeur sans FMA : ignorer cette démonstration |
| Temps différents de ce fichier | Fréquence du processeur, autres programmes en cours : les temps varient ; les valeurs calculées doivent être identiques, sauf `reproSum` section 3 (dépend de l'exécution, volontairement) |
| `/bin/bash^M: bad interpreter` | Fins de ligne Windows : `sed -i 's/\r$//' RunAll` |
