import math
import numpy as np
import matplotlib.pyplot as plt



Na = 6.02214e23

rho = 0.0709 #g/cm3

e = 1.602176634e-19*1e9 #nC
e2 = 1.602176634e-19*1e3 #mC

M = 2.016 #g/mol

l = 5 #cm

fact = (2*Na*rho*l/(e*M))*1e-36 #pb-1/nC
fact2 = (2*Na*rho*l/(e2*M))*1e-36 #pb-1/mC



print(fact, "pb-1/nC")
print(fact2, "pb-1/mC")

Q_fall18in = 26.312 + 4.0 + 5.355
Q_fall18out = 11.831 + 20.620
Q_spring19 = 45.994

Q = Q_fall18in  #mC


L_int = fact2*Q #pb-1

print("La luminosité intergré vaut")
print(L_int, "pb-1")




