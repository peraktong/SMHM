import numpy as np
import pickle
import math
import scipy
import matplotlib.pyplot as plt
from termcolor import colored

import pickle

import emcee
import scipy.optimize as op
import pickle
import matplotlib
import matplotlib.pyplot as plt
path = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/Behroozi_2013_data/sfh_stats/stats_a0.334060_absolute_mpeak_10.000_all.dat"




pkl_file = open('z_target.pkl', 'rb')
z_target = pickle.load(pkl_file)
pkl_file.close()




pkl_file = open('Behroozi_revised_M13.0.pkl', 'rb')
B_revised = pickle.load(pkl_file)
pkl_file.close()

# a
us = 1/(1+z_target[::-1])
b = B_revised[:,0]


# attention!! Add this mask to B_revised_data first:





print(1/b-1)
print(len(b))

count = 0

index_us = []
index_b = []

for i in range(0,len(b)):

    index = np.where(abs(us-b[i])<0.004)

    # Only use a>0.2 data:

    if (len(us[index])==1)&(b[i]>0.2):
        index = np.array(index).ravel()
        print(index)

        count+=1
        print(us[index], b[i],us[index]-b[i])
        index_us.append(index)
        index_b.append(i)

index_us = np.array(index_us).ravel()
index_b = np.array(index_b)


#print(index_us)
#print(index_b.ravel())

print(us[index_us])

print(b[index_b])

b_indexed = b[index_b]

delta_distance = []

for k in range(0,len(b_indexed)-1):
    delta_distance.append(b_indexed[k+1]-b_indexed[k])

delta_distance.append(b_indexed[-1]-b_indexed[-2])
delta_distance = np.array(delta_distance)

print("delta_distance_b")
print(delta_distance)

delta_distance_nor = delta_distance*len(delta_distance)/np.nansum(delta_distance)

print(delta_distance_nor)


# save them:


output = open("a_index_us.pkl", 'wb')
pickle.dump(index_us, output)
output.close()

output = open("a_index_Behroozi.pkl", 'wb')
pickle.dump(index_b, output)
output.close()

output = open("a_weight_factor.pkl", 'wb')
pickle.dump(delta_distance_nor, output)
output.close()






