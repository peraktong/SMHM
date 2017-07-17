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

    def read_tree(self):

        # From 11.0 to 15 with bin size=0.1. Your f_con works well for 9.75 to 12.55...


        merger_tree = []

        mi = 11.0

        for i in range(0, 41):
            number = mi + 0.1 * i
            number = "{0:.1f}".format(number)

            merger_tree.append("M" + str(number) + "=None")

        # Split
        merger_tree = dict(s.split("=") for s in merger_tree)
        # print(index_array)

        merger_tree_all = []

        for i in range(0, 41):
            number = mi + 0.1 * i

            pkl_file = open("Bolshoi_tree_" + str(number) + ".pkl", 'rb')
            array_i = pickle.load(pkl_file)
            pkl_file.close()

            merger_tree["M" + str(number)] = array_i

            merger_tree_all.extend(array_i)

        merger_tree_all = np.array(merger_tree_all)

        # prepare merger_tree for the future fitting with HMF weighted.
        self.merger_tree = merger_tree
        self.merger_tree_all = merger_tree_all

    def interpolate_to_behroozi(self):

        ## Only fit for M11.0 to 14.0 with 0.5 dex bin

        min = 11.0
        for k in range(0,7):

            mass = min+k*0.5
            mass = "{0:.1f}".format(mass)

            # READ BOLSHOI

            pkl_file = open("Bolshoi_tree_" + str(mass) + ".pkl", 'rb')
            Bolshoi_tree_i = pickle.load(pkl_file)
            pkl_file.close()


            # READ BEHROOZI

            # fusion: a+s+h+ds+dh+f_con

            path_behroozi = "Behroozi_revised_cen_M_11.0.pkl".replace("11.0",str(mass))


            pkl_file = open(path_behroozi, 'rb')
            Behroozi_i = pickle.load(pkl_file)
            pkl_file.close()


            ## Let's interplate everyhing onto two array: a_behroozi and Mh(N_tree*len(a_nbehroozi))

            a_target = Behroozi_i[:,0]

            output = open("a_Behroozi.pkl", 'wb')
            pickle.dump(a_target, output)
            output.close()



            #plt.plot(a_target, log10(Behroozi_i[:, 2]), "x")

            # use numpy.interp

            for j in range(0,np.min([2,len(Bolshoi_tree_i)])):

                tree_i = np.array(Bolshoi_tree_i[j],dtype=float)

                # Remember to reverse them...
                ai = tree_i[:,0][::-1]

                mi = tree_i[:,1][::-1]

                plt.plot(ai,mi, "r",alpha=0.3)

                #plt.plot(ai,"ro")
                #plt.plot(a_target, "ko")

                # interpolation

                mi_interpolated = np.interp(x=a_target,xp=ai,fp=mi)

                plt.plot(a_target,mi_interpolated, "b",alpha=0.3)

        plt.show()







model = Interpolate_median_Bolshoi_tree()

# model.read_tree()

model.interpolate_to_behroozi()

