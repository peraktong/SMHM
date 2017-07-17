import numpy as np
from scipy import integrate
import time
from termcolor import colored
import math
import os
from helpers.SimulationAnalysis import SimulationAnalysis, readHlist, iterTrees, getMainBranch
import pickle
import emcee
import scipy.optimize as op
import matplotlib.pyplot as plt


### start:

def log10(x):
    if x > 0:
        return math.log10(x)
    else:
        return -np.inf


def exp(x):
    try:
        return math.exp(x)
    except:
        return np.inf


exp = np.vectorize(exp)
log10 = np.vectorize(log10)

omega_m = 0.272
omega_gamma = 0.728


def E(z):
    return (omega_m * (1 + z) ** 3 + omega_gamma) ** 0.5


def a_to_time_Hogg(a):
    z = 1 / a - 1

    result = integrate.quad(lambda x: ((1 + x) * E(x)) ** (-1), z, np.inf)

    return result[0]


a_to_time_Hogg = np.vectorize(a_to_time_Hogg)

# Quenching fraction:
tinker = np.loadtxt("tinker_fqcen_SDSS_M9.7.dat")

ms_tinker = tinker[:, 0]
fraction_tinker = tinker[:, 1]
error_tinker = tinker[:, 2]


# start from logMh=11 halos:

class Interpolate_median_Bolshoi_tree():
    def update_kwargs(self, kwargs):

        self.kwargs = kwargs

    def interpolate_to_behroozi(self):

        ## Interpolate for M11.0 to 14.0 with 0.1 dex bin

        # read a_target:


        pkl_file = open("a_Behroozi.pkl", 'rb')
        a_target = pickle.load(pkl_file)
        pkl_file.close()

        min = 11.0
        for k in range(0,31):

            mass = min+k*0.1
            mass = "{0:.1f}".format(mass)


            pkl_file = open("Bolshoi_tree_" + str(mass) + ".pkl", 'rb')
            Bolshoi_tree_i = pickle.load(pkl_file)
            pkl_file.close()
            # use numpy.interp

            mh_int_array = []

            for j in range(0,len(Bolshoi_tree_i)):

                tree_i = np.array(Bolshoi_tree_i[j],dtype=float)

                # Remember to reverse them...

                # interpolation

                mh_int_array.append(np.interp(x=a_target,xp=tree_i[:,0][::-1],fp=tree_i[:,1][::-1]))

            mh_int_array = np.array(mh_int_array)

            # print(mh_int_array.shape)


            save_path = "Bolshoi_tree_interpolated_" + str(mass) + ".pkl"

            output = open(save_path, 'wb')
            pickle.dump(mh_int_array, output)
            output.close()

            # save median halos!


            output = open(save_path.replace("interpolated_","interpolated_median_"), 'wb')
            pickle.dump(np.nanmedian(mh_int_array,axis=0), output)
            output.close()

            # check median



model = Interpolate_median_Bolshoi_tree()

# model.read_tree()

model.interpolate_to_behroozi()

