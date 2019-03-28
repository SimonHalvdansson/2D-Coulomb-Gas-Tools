# -*- coding: utf-8 -*-

from scipy import *
from pylab import *
from scipy.integrate import *
from mpl_toolkits.mplot3d import Axes3D
from pylab import meshgrid
from decimal import *
import pickle
import time
import numpy as np
import matplotlib.pyplot as plt
from MLlib import *
from functools import lru_cache
import scipy.special as spec
import sys
import cProfile
from mpmath import *

mp.dps = 800

n = 600
k_s = 2

def Q(zeta):
   return abs(zeta**k_s - k_s**(-1/2))**2 - 1/k_s# -2*1/n*log(abs(zeta))
 
def inner(p, q):
   return nquad(lambda x, y: (lambda z: (p(z)*q(z).conjugate()*e**(-n*Q(z))).real)(x+y*1j), [[-4, 4], [-4, 4]])[0]/pi

def t():
   return time.time()

def innerzij_hypergeometric(i, j):
    a1 = gamma((i+1)/k_s)/factorial(round((i-j)/k_s))
    a2 = (mpf(n)/sqrt(mpf(k_s)))**((i-j)/k_s)
    return a1*a2/k_s/(mpf(n)**((i+1)/k_s))*hyp1f1((i+1)/k_s, 1 + round((i-j)/k_s), n/k_s)
   
def innerzij(i, j):
   #i = i+1
   #j = j+1
   if not ((i-j) % k_s == 0):
      return 0
    
   return innerzij_hypergeometric(max(i, j), min(i, j))
 
def genInnprods(start = 0, end = n, savefile = False):
    c = 0

    res = zeros(n)

    for i in range(0, n):
        for j in range(i, n):
            v = innerzij(i, j)
            res[i, j] = v
            res[j, i] = v
            c += 1
            if c % 500 == 0:
                print('{} % inner products calculated'.format(round(100*100*c/(n*(n+1)/2))/100))
        
    return res
   
def p(j, z):
   s = 0
   for i in range(j+1):
      s += float(ps[j][i])*z**i
   return s

def R_guess(z):
   return 1/2*math.erfc(-sqrt(2)*z.real)*abs(z)**2#/(1+1/(abs(z)**2))

def R(z):
   a = z.real
   b = z.imag
   return (1/(1.5+abs(z)**2))*1/2*math.erfc((-2*a**2 + 2*b**2)/sqrt(2))*E(abs(z)**2, 1/k_s, 1/k_s)*e**(-abs(z)**(2*k_s))

def E(z, a, b):
   return k_s*ml(z, 1/k_s, 1/k_s)

def frac(z):
   return R_guess(z)/Rn(z)#*(abs(z**2))

def Rn(z):
   zeta = z*n**(-1/(2*k_s))
   
   s = 0
   for j in range(n):
      s += abs(p(j, zeta)*e**(-0.5*n*Q(zeta)))**2
   return s*n**(-1/k_s)#/(abs(zeta)**2)

def operate_on_Narray(A, B, function):
    try:
        return [operate_on_Narray(a, b, function) for a, b in zip(A, B)]
    except TypeError as e:
        return function(A, B)

def plotComplex(func, a, b, c, d, step, contour=True):
    x = np.arange(a, b, step)
    y = np.arange(c, d, step)
    
    X, Y = meshgrid(x, y)
    Z = operate_on_Narray(X, Y, lambda a, b: func(a+b*1j))
    fig = plt.figure(figsize=(13,10))
    
    if contour:
        plt.contour(X, Y, Z, 100)
        plt.colorbar()
    else:    
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z, cstride=1, rstride=1)    
    plt.show()

def plotPj(j, a=-1.5, b=1.5, c=-1, d=1, step=0.01, contour=True):
    plotComplex(lambda z: abs(p(j, z))**2, a, b, c, d, step, contour)

def plotQ(a = -1.5, b = 1.5, c = -1, d = 1, step= 0.01, contour=False):
    plotComplex(Q, a, b, c, d, step, contour)
    
def plotRn(a=-0.5, b=0.5, c=-0.5, d=0.5, step=0.05, contour=True):
    plotComplex(Rn, a, b, c, d, step, contour)
    
def plotR(a=-0.5, b=0.5, c=-0.5, d=0.5, step=0.05, contour=True):
    plotComplex(R, a, b, c, d, step, contour)

def plotRguess(a=-0.5, b=0.5, c=-0.5, d=0.5, step=0.05, contour=True):
    plotComplex(R_guess, a, b, c, d, step, contour)
    
def plotFrac(a=-0.5, b=0.5, c=-0.5, d=0.5, step=0.05, contour=True):
    plotComplex(frac, a, b, c, d, step, contour)

def save(obj, name):
   f = open('{}.pckl'.format(name), 'wb')
   pickle.dump(obj, f)
   f.close()
   
def load(name):
   f = open('{}.pckl'.format(name), 'rb')
   obj = pickle.load(f)
   f.close()
   
   for i in range(len(obj)):
      for j in range(len(obj[i])):
         obj[i][j] = mpf(str(obj[i][j]))
         
   return obj

def mult(a, b, la = -1, lb = -1):
   """
   Takes two n dim vectors representing polynomial coefficients
   and computes scalar product wrt our scalar product
   """
   if la == -1:
      la = len(a)
      lb = len(b)
   
   s = 0
   
   for i in range(la):
      for j in range(lb):
         s += innprods[i, j]*a[i]*b[j]
   
   return s

def divide(vec, dec):
   return array([x/dec for x in vec])

def genPs():
   """
   MGS implementation
   """
   vs = []
   ps = []
   for i in range(n):
      vi = [0]*(i+1)
      vi[i] = 1
      vs.append(vi)

   for i in range(n):
      print('Orthogonalizing polynomial {}'.format(i))
      
      rii = sqrt(mult(vs[i], vs[i], i+1, i+1))
      if rii < 0:
         print('whoops')
      ps.append(divide(vs[i], rii))
      for j in range(i+1, n):
         rij = mult(ps[i], vs[j], i+1, j+1)
         for k in range(i+1):
            vs[j][k] -= rij*ps[i][k]
            
   return ps
      

def checkOrthogonality():
   i = 0
   m = 0
   s = 0
   A = zeros(n)
   for x in range(n):
      for y in range(x, n):
         i += 1
         if i % 100 == 0:
            print('Checked {}%'.format(100*i/(n**2/2)))
         A[x, y] = mult(ps[x], ps[y], x+1, y+1)  
         if x == y:
            A[x, y] -= 1

         if abs(A[x, y]) > m:
             m = abs(A[x, y])

         if abs(A[x, y] < 0.01):
            A[x, y] = 0
         A[y, x] = A[x, y] 
         s += 2*A[x, y]
            
   print(s)
   print(m)
   return A

def plotPjZero():
   xs = list(range(n))
   ys = []
   s = 0
   
   for x in xs:
      ys.append(float(ps[x][0]**2))
      s += ps[x][0]**2
                
   print(s*n**(-1/k_s))
   
   

#innprods = genInnprods()
#ps = genPs()  
