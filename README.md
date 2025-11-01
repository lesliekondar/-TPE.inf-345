
import sympy as sp import numpy as np

Définition de la méthode de dichotomie def dich(f, a, b, e): m = e x = a y = b s = x

if f(y) == 0: return y

while abs(y - x) >= m and f(x) * f(y) < 0: s = (x + y) / 2 if f(x) * f(s) <= 0: y = s else: x = s return s

Méthode de balayage
def bal(f, a, b, e): x = a while f(x) * f(x + e) > 0: x += e if x > b: raise ValueError("Pas de changement de signe dans l'intervalle.") s = x + (e / 2) return s

Méthode de Lagrange
def lag(f, a, b, e): x = a y = b s = x if f(y) == 0: return y

while abs(y - x) >= e and f(x) * f(y) < 0: s = ((x * f(y)) - (y * f(x))) / (f(y) - f(x)) if f(x) * f(s) < 0: y = s

else: x = s return s

Méthode de Newton-Raphson
def new(f, x0, e): x = sp.Symbol('x') f_expr = f(x) df = sp.diff(f_expr, x)

f_lamb = sp.lambdify(x, f_expr, 'numpy') df_lamb = sp.lambdify(x, df, 'numpy')

s = x0 while abs(f_lamb(s)) > e: s = s - f_lamb(s) / df_lamb(s) return s

Méthode de dichotomie pour les #fonctions à plusieurs variables
def dichsup(f, a, b, e, mi=100): fa = f(*a) fb = f(*b) if fa * fb > 0: raise ValueError("f(a) et f(b) doivent avoir des signes contraires.")

for i in range(mi): m = (np.array(a) + np.array(b)) / 2 fm = f(*m) if abs(fm) < e: return m if fa * fm < 0: b = m fb = fm else: a = m fa = fm return m

import numpy as np

def newc(f, jf, a, e=1e-6, nm=100):

x = np.array(a, dtype=float)

for i in range(nm): fx = np.array(f(*x)) jfx = np.array(jf(*x))

if np.linalg.norm(fx) < e: return x

try: delta = np.linalg.solve(jfx, -fx) # Résolution J(x)Δx = -f(x) except np.linalg.LinAlgError: raise ValueError("Jacobienne non inversible à l'étape", i)

x = x + delta

if np.linalg.norm(delta) < e: return x

raise ValueError("Non convergence après {} itérations".format(nm))

import numpy as np

def gauss_pivot(A, b): """ Résout un système linéaire Ax = b par la méthode du pivot de Gauss. A : matrice carrée b : vecteur des constantes """ A = A.astype(float) b = b.astype(float) n = len(b)

Étape d'élimination
for k in range(n - 1): for i in range(k + 1, n): if A[k, k] == 0: raise ZeroDivisionError("Pivot nul, méthode échouée.") facteur = A[i, k] / A[k, k] A[i, k:] = A[i, k:] - facteur * A[k, k:] b[i] = b[i] - facteur * b[k]

Substitution arrière
x = np.zeros(n) for i in reversed(range(n)): x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i] retrn x

def gauss_jordan(A, b): """ Résout un système Ax = b par la méthode du pivot de Gauss-Jordan. Retourne la solution x. """ A = A.astype(float) b = b.astype(float) n = len(b)

Augmenter la matrice A avec b
aug =np.hstack((A,b.reshape(-1,1)))

for i in range(n): # Normaliser la ligne i if aug[i, i] == 0: raise ZeroDivisionError("Pivot nul détecté.") aug[i] = aug[i] / aug[i, i]

Annuler les autres lignes
for j in range(n): if i != j: aug[j] = aug[j] - aug[j, i] * aug[i]

return aug[:, -1] for j in range(n): if i != j: aug[j]-=aug[j,i]*aug[i]

return aug[:, -1]

def crout(A, b): """ Méthode de Crout pour résoudre Ax = b via décomposition LU """ n = len(b) L = np.zeros((n, n)) U = np.identity(n)

for j in range(n): for i in range(j, n): L[i, j] = A[i, j] - sum(L[i, k] * U[k, j] for k in range(j)) for i in range(j + 1, n): U[j, i] = (A[j, i] - sum(L[j, k] * U[k, i] for k in range(j))) / L[j, j]

Résolution de L*y = b
y = np.zeros(n) for i in range(n): y[i] = (b[i] - sum(L[i, j] * y[j] for j in range(i))) / L[i, i]

Résolution de U*x = y
x = np.zeros(n) for i in reversed(range(n)): x[i] = (y[i] - sum(U[i, j] * x[j] for j in range(i + 1, n))) / U[i, i]

return x

implémentation
import sympy as sp import numpy as np from main import dich, lag, new, bal, dichsup

a = float(input("Entrer la borne inférieure de votre intervalle : ")) b = float(input("Entrer la borne supérieure de votre intervalle : ")) e = float(input("Entrer la marge d'erreur : "))

n = int(input("Nombre de dimensions (1 pour 1D, >1 pour multi-dimensions) : "))

x = sp.symbols(' '.join([f'x{i+1}' for i in range(n)]))

f_str = input(f"Entrez la fonction f({', '.join([f'x{i+1}' for i in range(n)])}): ")

f_expr = sp.sympify(f_str)

#Conversion en fonction numérique utilisable avec numpy f = sp.lambdify(x, f_expr, 'numpy')

#Méthodes unidimensionnelles (seulement si n=1) if n == 1: s1 = dich(f, a, b, e) print("Solution par la méthode de dichotomie :", s1)

s2 = bal(f, a, b, e) print("Solution par la méthode de balayage :", s2)

s3 = lag(f, a, b, e) print("Solution par la méthode de Lagrange :", s3)

s4 = new(f, a, e) print("Solution par la méthode de Newton :", s4) else: # Méthodes multidimensionnelles a_vect = input(f"Entrez le vecteur a (séparé par des virgules, {n} valeurs) : ") b_vect = input(f"Entrez le vecteur b (séparé par des virgules, {n} valeurs) : ")

Conversion en numpy arrays de float
a_vect = np.array([float(val) for val in a_vect.split(",")]) b_vect = np.array([float(val) for val in b_vect.split(",")])

Vérification que les dimensions sont correctes
if len(a_vect) != n or len(b_vect) != n: print(f"Erreur : les vecteurs a et b doivent avoir exactement {n} dimensions.") else: # Appel de la méthode dichotomique multidimensionnelle s_multidim = dichsup(f, a_vect, b_vect, e) print("Solution par la méthode de dichotomie multidimensionnelle :", s_multidim)



