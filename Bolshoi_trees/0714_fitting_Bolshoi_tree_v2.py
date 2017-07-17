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

class Fit_Bolshoi_tree():
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


    def read_tree_index(self,index):

        # index from 11.0 to 15.0

        self.index = index

        merger_tree_all = []

        pkl_file = open("Bolshoi_tree_" + str(index) + ".pkl", 'rb')
        array_i = pickle.load(pkl_file)
        pkl_file.close()


        self.merger_tree_all = array_i

    def delta_Ms_to_Ms(self, delta_Ms):

        array_all = np.tile(delta_Ms, (len(delta_Ms), 1))

        array_all = np.tril(array_all, k=0)
        sum = np.sum(array_all, axis=1)

        return sum

    def Mh_to_delta_Mh(self,Mh):

        up = np.append(Mh,[0])
        low = np.append([0],Mh)

        delta_Mh = np.array((up-low)[:-1])

        # Attention!!
        delta_Mh[0] = 0

        return delta_Mh


    def t_to_delta_t(self, t):

        up = np.append(t, [0])
        low = np.append([0], t)

        delta_t = np.array((up - low)[:-1])

        # Attention!!
        delta_t[0] = 0

        return delta_t

    def f_con_Jason(self, a, Mh):

        kwargs = self.kwargs

        f0_log = kwargs["f0_log"]

        f0 = 10 ** (f0_log)

        A1 = kwargs["A1"]
        A2 = kwargs["A2"]
        A3 = kwargs["A3"]
        mht = kwargs["mht"]
        A4 = kwargs["A4"]
        A5 = kwargs["A5"]

        A6 = kwargs["A6"]

        # I think we need to calculate M_critical from reading the data:

        At = A3 * (log10(Mh)) ** 2 + A4 * (log10(Mh)) + A5

        # fix a_s
        As = A6

        f_con_original = f0 * ((Mh / mht) ** A1) * (log10((Mh / mht))) ** A2 * exp(-(At - a) ** 2 / As)

        # Quenched fraction part:

        mhc_log = kwargs["mhc_log"]

        mhc = 10 ** (mhc_log)

        sigmah = kwargs["sigmah"]

        alphah = kwargs["alphah"]

        f_q = (exp((mhc_log - log10(Mh)) / sigmah)) ** alphah

        if Mh > mhc:
            return f_q * f_con_original
        else:
            return f_con_original

            # Let's return:

    def fit_quenched_fraction(self):

        mhc = 10 ** (self.kwargs["mhc_log"])

        f_con_Jason = np.vectorize(self.f_con_Jason)

        merger_tree_all = self.merger_tree_all

        # Let's fit the tree:

        # Fusion = a+Mh+Ms+f_con+SFR
        fusion = []

        judge_array = np.zeros(19)

        length_array = np.zeros(19)

        for num in range(0, len(merger_tree_all)):

            tree_i = np.array(merger_tree_all[num], dtype=float)

            ### Attention!! Inverse!
            a = tree_i[:, 0][::-1]
            Mh = tree_i[:, 1][::-1]

            f_con_array = f_con_Jason(a=a, Mh=Mh)

            delta_Mh_array = self.Mh_to_delta_Mh(Mh=Mh)

            delta_Ms_array = f_con_array * delta_Mh_array

            Ms = self.delta_Ms_to_Ms(delta_Ms=delta_Ms_array)

            # calculate SFR:


            t = a_to_time_Hogg(a=a)

            """

            delta_t = self.t_to_delta_t(t=t)


            sfr_array = delta_Ms_array/delta_t


            """

            index_quenching = np.argmin(abs(Mh - mhc))

            final_log_ms = log10(Ms[-1])

            # index_tinler is between 0 and 18


            index_tinker = np.argmin(abs(final_log_ms - ms_tinker))

            length_array[index_tinker] += 1

            # calculate judge


            bottom = (Ms[index_quenching] - Ms[index_quenching - 1]) / (
                t[index_quenching] - t[index_quenching - 1])
            top = (Ms[-1] - Ms[-2]) / (t[-1] - t[-2])

            if bottom == 0:
                judge = 100
                print("fails")
            else:
                judge = top / bottom

            if judge < 0.1:
                judge_array[index_tinker] += 1
            else:
                # judge_array[index_tinker] += 0
                none = 1

        ## Let's put a mask:

        mask = length_array > 0
        self.mask = mask
        # fitting:


        quenched_fraction_array = judge_array[mask] / length_array[mask]
        self.quenched_fraction_array = quenched_fraction_array

        ivar = (error_tinker[mask]) ** (-2)
        ### Need to normalize ivar!!!

        ivar = ivar * len(ivar) / np.sum(ivar)

        chi_quenched = np.sum((quenched_fraction_array - fraction_tinker[mask]) ** 2 * ivar)
        self.chi_quenched = chi_quenched

        print("chi-squared")
        print(chi_quenched)
        return chi_quenched


## Initialize;

######## Index
# read trees
index = "11.5"


kwargs = {"MCMC": None, "A1": None, "A2": None, "A3": None, "f0_log": None, "A4": None, "A5": None, "A6": None,
          "mht": None,
          "mst": None, "zc": None, "sigmaz": None,
          "alphaz": None, "mhc_log": None, "msc_log": None, "sigmah": None, "sigmag": None, "alphah": None,
          "alphag": None}

# Initial values:
# From scipy best fitting model. Not sure whether to use this or use EMCEE.


kwargs["f0_log"] = log10(0.00251552734235)

kwargs["A1"] = -5.00156921458
kwargs["A2"] = 36.7909234394
kwargs["A3"] = -1.45406047315

kwargs["A4"] = 31.1501956219

kwargs["A5"] = -163.028606123

kwargs["A6"] = 2.83787333037

kwargs["mht"] = 10 ** (8)
kwargs["mst"] = 10 ** (8)

###### Now we add Jeremy's stellar quenching model:

# Here we fit fo, mst,sigmag and alphag.


kwargs["mhc_log"] = 11.66

kwargs["sigmah"] = 0.52

# 3.63
kwargs["alphah"] = 3.6

model = Fit_Bolshoi_tree()

# update kwargs:

model.update_kwargs(kwargs=kwargs)



model.read_tree_index(index=index)
# model.read_tree()
# Let's do the fitting and return the chi-squared(quenched)

# model.fit_quenched_fraction()

# Here comes the fitting:


#### Index

counter = 0


def lnlike(theta, x, y):
    global counter

    # introduce our model here.

    # Theta is a tuple
    # f0_log, mhc_log,sigmah,alphah = theta

    mhc_log, sigmah, alphah = theta

    # kwargs["f0_log"] = f0_log
    kwargs["mhc_log"] = mhc_log

    kwargs["sigmah"] = sigmah
    kwargs["alphah"] = alphah

    model.update_kwargs(kwargs=kwargs)

    # Let's change the model from simple linear to our complex model.

    # y_model is the value from our method

    # Now we define chi as the abs of the fraction of the quenched galaxies


    counter = counter + 1
    print("Doing %d" % counter)
    print(mhc_log, sigmah, alphah)

    chi = model.fit_quenched_fraction()
    # chi = model.calculate_M_stellar()

    # if you want to fit simultaneously. Use this

    print("Tinker fraction")
    print(fraction_tinker[model.mask])

    print("Our fraction")
    print(model.quenched_fraction_array)

    print(colored("chi", "red"))
    print(colored(chi, "red"))

    return -0.5 * (chi)


# Set the range of the parameters
def lnprior(theta):
    mhc_log, sigmah, alphah = theta

    #
    # -3 < f0_log < 0 and
    if 8 < mhc_log < 13 and 0 < sigmah < 3 and -5 < alphah < 5:
        return 0.0
    return -np.inf


# The final ln posterior probability is the sum of the ln prior and ln likelihood



# Then, let's do some MCMC:

start_time = time.time()
print("doing emcee")


# The final ln posterior probability is the sum of the ln prior and ln likelihood


# Define the initial condition of the MCMC chain: Position/initial values

def lnprob(theta, x, y):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf, str(lnlike(theta, x, y) / (-0.5))
    return lp + lnlike(theta, x, y), str(lnlike(theta, x, y) / (-0.5))


# Define the initial condition of the MCMC chain: Position/initial values

ndim, nwalkers = 3, 24

kwargs["MCMC"] = 1

# Or you can replace result["x"] with your initial values
# pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
pos = [
    [kwargs["mhc_log"], kwargs["sigmah"], kwargs["alphah"]] + 1e-3 * np.random.randn(ndim) for i in range(nwalkers)]

# Set up the MCMC chain
x = 1
y = 1
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y), threads=8)

print("running MCMC")
sampler.run_mcmc(pos, 200)

# Now we have an array with dimension of 100*500*3: 100 walkers, 500 steps and 3 parameters

result = sampler.chain
blobs = sampler.blobs

blobs = np.array(blobs)

print(result)

# save it:
# Let's fit M10.0 M10.5 --- M13.0 first

# save_path = "Bolshoi_trees_500_best_fit_v2.pkl"
save_path = "Bolshoi_trees_500_best_fit_v2.pkl".replace("best_fit_","best_fit_"+str(model.index))
save_path_backup = save_path

save_path = save_path.replace("best_fit_", "EMCEE")

output = open(save_path, 'wb')
pickle.dump(result, output)
output.close()

save_path = save_path.replace("EMCEE", "EMCEE_chi_")
output = open(save_path, 'wb')
pickle.dump(blobs, output)
output.close()

######## scipy fitting:


# scipy fitting!!

chi = blobs.T

chi = np.array(chi, dtype=float)
index = np.argwhere(chi == np.nanmin(chi))

length_index = len(index)
# print("length_index")
# print(length_index)


i, j = index[0]

mhc_log, sigmah, alphah = result[i, j, :]
print("chi for different index")
print(mhc_log, sigmah, alphah)
print(chi[i, j])

# kwargs["f0_log"] = f0_log
kwargs["msc_log"] = mhc_log
kwargs["sigmag"] = sigmah
kwargs["alphag"] = alphah

# End doing EMCEE



x = 1
y = 1

nll = lambda *args: -lnlike(*args)

print("doing scipy opt")
# Only fit f0,A1,A2,A3
# I adjust the likelihood function

result = op.minimize(nll,
                     [kwargs["mhc_log"], kwargs["sigmah"], kwargs["alphah"]],
                     args=(x, y))
print("finish doing scipy opt")

print("result x")
print(result["x"])

# Path:


# mst,sigmag,alphag = result["x"]

output = open(save_path_backup, 'wb')
pickle.dump(result["x"], output)
output.close()

stop_time = time.time()

print("Using %.2f seconds" % (stop_time - start_time))













