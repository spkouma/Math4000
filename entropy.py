import numpy as np
from math import factorial
import math as math
import matplotlib.pyplot as plt 


#constant γ, provided in Madras et al
gamma_3=1.162
print("gamma_3=",gamma_3)


#constant μ, provided in Madras et al
mu_3=4.683
print("mu_3=",mu_3)


#Amplitude constant, provided by res paper
A_3=1.205
print("A_3=",A_3)


#number of mers on surface, probably a loop param?
j=0
print("j=",j)


# number of sites on surface lattice
a=int((2300*10**((-9)*2))/(0.148*10**((-9)*2)))
print("a=",a)


#number of single-mers in a large-mer
n=6200/27
print("n=",n)


# z_vol arb selected size, say a/3 or a/10
#Z_vol=?
#size of a Mer
M=27
print("M=",M)


#__verify what numbers we shall be usig from the calculations I sent madras__
#dont confuse this with my notation, this is total number of... 
#...Short,Longs in the system mine is the count of their components per chain
#number nM-mers or Large-mers in the system, assume n_L >>a
n_L=2.4092*10**22
print("n_L=",n_L)


#number of M-mers or Short-mers in the system, assume n_s ~ n*n_L??
n_S=2.2887*10**22
print("n_S=",n_S)


#coordination number, number of directions 
q=4
print("q=",q)


#boltzmann constant
k_B=1.38064852*10**(-23)
print("k_B=",k_B)

#Z_vol, arbitrarily set due to lack of volumetrics
Z_vol=(a**2)/10
print("Zvol=",Z_vol)


#Approximation function
def f(j):
    return (Z_vol*A_3*M**(gamma_3-1)*((q-1)**2)*(a-j))/(a*q*n_S)

print("f(0.8)",f(0.8))


#entropy
def S(j):
    return -k_B*np.log(f(j))


# In[] Entropy Graph

x = np.arange(0, (a/(n*M)))
y = S(x)
plt.xlabel('j surface sites occupied by Large')
plt.ylabel('$J\cdot K^{-1}$')
plt.title('Relative Entropy Approximation') 

# Plot the points using matplotlib 
plt.plot(x, y) 
plt.show() 
