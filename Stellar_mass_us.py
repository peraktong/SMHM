import numpy as np
import math


# here we set the global fb

fb = 0.15


class SMHM_relation():

    def stellar_mass_per_interval(self,f_con,delta_Mh):

        fb =0.15

        return f_con*fb*delta_Mh

    def f_con_z(self,z,f0):
        z0 = 1
        gamma1 = -3
        gamma2 = 0

        if z>z0:
            return f0*((1+z)/(1+z0))**gamma1

        else:
            return f0*((1+z)/(1+z0))**gamma2

    def calculate_mass(self,f0):
        # input array: f_con at zi, deltaMh at zi
        # return an array of M*

        M_stellar = []

        for i in range(0,len(self.z_array)):

            M_stellar.append(self.stellar_mass_per_interval(f_con=self.f_con_z(z=self.z_array[i],f0=f0),delta_Mh=self.delta_mh_array[i]))

        return np.array(M_stellar)


    # I think what we need is to read merger trees and calculate delta_Mh
    # f_con is fit to the data in Behroozi 2013.

    def read_data(self):
        none =1

        f_con = np.array([1,2,3,4])
        self.f_con_array = f_con

        self.delta_mh_array = np.array([1.2*10**12,1.4*10**12,1.5*10**12,1.7*10**12,])

        # Let's generate an z array from 0 to 6
        self.z_array = np.array([0.1,0.4,0.7,0.8])


model = SMHM_relation()

model.read_data()
final = model.calculate_mass(f0=4)
print(final)


