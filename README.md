
Entrer la borne inférieure de votre intervalle : 1
Entrer la borne supérieure de votre intervalle : 2
Entrer la marge d’erreur : 0.0001
Entrer la fonction f(x) : x**2 - 2

# -TPE.inf-345
test des fonctions
import sympy as sp
from sympy.abc import x
import numpy as np

# Importation des méthodes numériques depuis ton autre fichier
from main_071658 import *

# ============================================
# Application principale
# ============================================

# Saisie de l’intervalle et de la marge d’erreur
a = float(input("Entrer la borne inférieure de votre intervalle : "))
b = float(input("Entrer la borne supérieure de votre intervalle : "))
e = float(input("Entrer la marge d’erreur : "))

# Saisie de la fonction
fs = input("Entrer la fonction f(x) : ")

# Conversion en fonction sympy et fonction numérique (numpy)
f = sp.sympify(fs)
fn = sp.lambdify([x], f, 'numpy')

# ============================================
# Appel des différentes méthodes de calcul
# ============================================

# Méthode de la dichotomie
s1 = dich(fn, a, b, e)
print("La solution par la méthode de dichotomie est :", s1)

# Méthode du balayage
s2 = bal(fn, a, b, e)
print("La solution par la méthode de balayage est :", s2)

# Méthode de Lagrange
s3 = lag(fn, a, b, e)
print("La solution par la méthode de Lagrange est :", s3)

# Méthode de Newton
s4 = new(fn, a, e)
print("La solution par la méthode de Newton est :", s4)