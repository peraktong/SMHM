import numpy as np
import math

# in B_2013

# use args


class SMHM_relation():


    def f (self,x,alpha,delta,gamma):

        return -math.log10((10**(alpha*x)+1))+delta*(math.log10(1+math.exp(x)))**gamma/(1+math.exp(10**(-x)))



    def SMHM(self,M1,Mh,epsilon,alpha,delta,gamma):

        f0 = self.f(0,alpha,delta,gamma)

        return math.log10(epsilon*M1)+self.f(math.log10(Mh/M1),alpha,delta,gamma)-f0



model = SMHM_relation()

M1 = 1
Mh =2
epsilon =3
alpha = 4
delta=5
gamma =6

args = (M1,Mh,epsilon,alpha,delta,gamma)

result = model.SMHM(*args)

print(result)
