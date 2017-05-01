import numpy as np

import matplotlib.pyplot as plt
import math

ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/2017.4_9_cosmology/data/data_all/"

def log10(x):
    return math.log10(x)




class read_and_plot():


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

    # Attention! This is for one Halo. This means the input z_array should be the z_array for one Halo!
    def calculate_mass(self,f0):
        # input array: f_con at zi, deltaMh at zi
        # return an array of M*

        z_array = self.z_array
        delta_mh_array = self.delta_mh_array

        M_stellar_interval = []

        for i in range(0,len(z_array)):

            M_stellar_interval.append(self.stellar_mass_per_interval(f_con=self.f_con_z(z_array[i],f0=f0),delta_Mh=delta_mh_array[i]))

        # add them togethrm

        M_stellar = []

        for j in range(0,len(M_stellar_interval)):

            # attention! small t small a big z: You Use Z in your calculation!
            M_stellar.append(np.sum(M_stellar_interval[j:]))

        self.M_stellar = np.array(M_stellar)
        print(M_stellar_interval)
        print(M_stellar)
        return M_stellar


    def input_array(self,z_array,delta_mh_array):

        self.z_array = z_array
        self.delta_mh_array = delta_mh_array




    def read_behroozi(self):

        path = ori + "behroozi_lookup.dat"

        # read file:

        result = np.array([0, 0, 0])
        for line in open(path, 'r'):
            data = line.rstrip().split(" ")

            result_i = []
            for i in data:
                if i == "":
                    none = 1
                else:
                    result_i.append(i)
            result_i = np.array(result_i, dtype=float)
            if len(result_i) == 3:
                result = np.vstack((result, result_i))

        print(result.shape)

        self.result_b = result


    def read_ours(self):

        path = ori + "parsed_tree1000.dat"

        # read file:

        result_tree = np.array([0, 0, 0])
        for line in open(path, 'r'):
            data = line.rstrip().split(" ")

            result_i = []
            for i in data:
                if i == "":
                    none = 1
                else:
                    result_i.append(i)

            result_i = np.array(result_i,dtype=float)
            if len(data)==3:
                result_tree = np.vstack((result_tree,result_i))
            else:
                none =1

        print(result_tree.shape)

        self.z = result_tree[1:101, 0]


        # Let's do the first 100 first:

        result_tree =result_tree[1:101,:]


        self.result_tree = result_tree


        # print(self.z)

        Mh = self.result_tree[:, 1]

        delta_Mh = []

        for j in range(0, len(Mh) - 1):
            delta_Mh.append(Mh[j] - Mh[j+1])

        delta_Mh.append(0)
        delta_Mh = np.array(delta_Mh)

        # fusion z Mh delta Mh  + (M_stellar)

        # Then calculate median
        fusion = np.c_[result_tree[:,0],delta_Mh]

        # calculate mass:

        # if u want to fit the function to the data in B_2013

        self.input_array(z_array=fusion[:,0],delta_mh_array=fusion[:,1])

        f0 = 0.01
        M_stellar = self.calculate_mass(f0=f0)

        print(np.array(fusion.shape),np.array(M_stellar).shape)

        fusion = np.c_[fusion,M_stellar]

        self.fusion = fusion

        """
                

        # Calculate Median:

        target = list(set(fusion[:,0]))


        for i in range(0,len(target)):

            print("Doing %.2f %%"%(i/len(target)*100))

            index = np.where(name == target[i])
            index = np.array(index)
            index = index.ravel()

            rms_old_i = np.std(VBARY[index])




        
        """


    def plot(self):


        # calculate delta_Mh

        fusion = self.fusion

        M_stellar = fusion[:,2]

        mask = M_stellar <=0
        M_stellar[mask] = 100





        # [log10(x) for x in M_stellar]
        plt.plot(1/(1+fusion[:,0]),[log10(x) for x in M_stellar],"ro")
        plt.show()



model = read_and_plot()

# model.read_behroozi()
model.read_ours()
model.plot()



"""

# Toy model for the three testing files
path = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/2017.4_9_cosmology/data/test_data/fulltree_100.1"

array = np.loadtxt(path)

print(array.shape)

#Plot them

# indx Mass S step-number desc_indx + scale factor

a = array[:,5]
M_stellar =  array[:,1]

print(a.shape)
print(M_stellar.shape)


plt.plot(a,[log10(x) for x in M_stellar],"ro")
plt.show()
"""

