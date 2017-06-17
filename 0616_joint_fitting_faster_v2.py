import numpy as np
import pickle
import math
import scipy
import matplotlib.pyplot as plt
from termcolor import colored
import pickle
import matplotlib
import emcee
import scipy.optimize as op


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

# Load index for fitting;
fb = 0.17

pkl_file = open('z_target.pkl', 'rb')
z_target = pickle.load(pkl_file)
pkl_file.close()

a_target = 1 / (1 + z_target)

a_target = a_target[::-1]

pkl_file = open('a_index_us.pkl', 'rb')
a_index_us = pickle.load(pkl_file)
pkl_file.close()

pkl_file = open('a_index_Behroozi.pkl', 'rb')
a_index_Behroozi = pickle.load(pkl_file)
pkl_file.close()

# weight factor in a
pkl_file = open('a_weight_factor.pkl', 'rb')
a_weight_factor = pickle.load(pkl_file)
pkl_file.close()


class MCMC():
    def __init__(self, kwargs):
        self.kwargs = kwargs

    def update_kwargs(self, kwargs):
        self.kwargs = kwargs

    def stellar_mass_per_interval(self, f_con, delta_Mh):

        fb = 0.15

        return f_con * fb * delta_Mh

    # f_con for Jeremy

    def f_con_Jeremy(self, z, Mh, Ms):

        kwargs = self.kwargs

        f0 = kwargs["f0"]

        z0 = 1
        gamma1 = -3
        gamma2 = 0

        if z > z0:

            return f0 * ((1 + z) / (1 + z0)) ** gamma1

        else:
            return f0 * ((1 + z) / (1 + z0)) ** gamma2

    # f_con for me

    def f_con_Jason(self, a, Mh, Ms):

        z = 1 / a - 1

        kwargs = self.kwargs

        f0 = kwargs["f0"]

        A1 = kwargs["A1"]
        A2 = kwargs["A2"]
        A3 = kwargs["A3"]
        mht = kwargs["mht"]
        mst = kwargs["mst"]
        A4 = kwargs["A4"]
        A5 = kwargs["A5"]

        A6 = kwargs["A6"]

        # I think we need to calculate M_critical from reading the data:

        # Shall we do this first?

        # Let A4 be a very small negative value

        # at = A4*log10(Mh)+A5

        # Or

        """
        
        
        at = A3*log10(self.median_Mh)+A4

        As = A5*log10(self.median_Mh)+A6

        
        """

        at = A4

        As = A6


        # zt = -0.3164*math.log10(self.M_h_max)+4.2442


        # Remember the critical a/z0
        # Now we use a and a varies from small to big:


        return f0 * ((Mh / mht) ** A1) * ((1. + z) ** A2)* (exp(-(at - a) ** 2 / As))


    ##### Quenching: Ignore it for now

    def f_q(self):

        return 1

    def calculate_M_stellar(self):

        kwargs = self.kwargs

        index = str(self.index)

        # Read_Behroozi's Data

        B_path = "Behroozi_revised_M11.0.pkl"

        B_path = B_path.replace("M11.0", index)

        data_path = "Median_a_Mh_300_M10.0.pkl"

        data_path = data_path.replace("M10.0", index)

        # fusion = np.c_[a, s, h,ds,dh,f_con_array]

        pkl_file = open(B_path, 'rb')
        Behroozi_fusion = pickle.load(pkl_file)
        pkl_file.close()

        pkl_file = open(data_path, 'rb')
        median_fusion = pickle.load(pkl_file)
        pkl_file.close()


        self.median_Mh = np.nanmedian(median_fusion[:,1])

        self.a_target = a_target

        # Let's do it

        a = a_target

        median_Mh = median_fusion[:,1]

        delta_Mh = [median_Mh[0]]

        Ms_now = 1

        Ms = [Ms_now]

        f_con_now = self.f_con_Jason(a[0], median_Mh[0], Ms_now)

        f_con_array = [f_con_now]

        for j in range(1, len(median_Mh)):
            delta_Mh_j = median_Mh[j] - median_Mh[j - 1]

            delta_Mh.append(delta_Mh_j)

            # calculate M_stellar
            f_con_array_j = self.f_con_Jason(a[j], median_Mh[j], Ms_now)

            f_con_array.append(f_con_array_j)

            delta_Ms = f_con_array_j * fb * delta_Mh_j * self.f_q()

            Ms_now = Ms_now + delta_Ms

            Ms.append(Ms_now)

        delta_Mh = np.array(delta_Mh)

        median_Ms = np.array(Ms)

        f_con_array = np.array(f_con_array)

        # Now we already have small a to big a:

        # Our fusion has a+ median_Ms, median_Mh median_f_con_array

        fusion = np.c_[a,median_Ms,median_Mh , f_con_array]


        # No scatter for now:


        log_median_Ms = [log10(x) for x in median_Ms]

        log_Behroozi_Ms = [log10(x) for x in Behroozi_fusion[:, 1]]

        log_median_Ms = np.array(log_median_Ms)
        log_Behroozi_Ms = np.array(log_Behroozi_Ms)

        chi = np.sum((log_median_Ms[a_index_us] - log_Behroozi_Ms[a_index_Behroozi]) ** 2*a_weight_factor)

        self.chi = chi


        return chi


    def fitting_separately(self, index):

        self.index = index

    def fitting_simultaneously(self):

        index_array = ["M11.0", "M11.5", "M12.0", "M12.5"]

        chi_all = []
        # No scatter for now

        for x in range(0, len(index_array)):
            self.fitting_separately(index=index_array[x])
            chi_all.append(self.calculate_M_stellar())

        chi_all = np.array(chi_all)


        return np.nanmean(chi_all)



#### Let's do it!




kwargs = {"MCMC": None, "A1": None, "A2": None, "A3": None, "f0": None, "A4": None, "A5": None, "A6": None, "mht": None,
          "mst": None, "zc": None, "sigmaz": None,
          "alphaz": None, "mhc": None, "msc": None, "sigmah": None, "sigmag": None, "alphah": None, "alphag": None}

# Initial values:
# different f0!

# Here As = zs.
kwargs["f0"] = 0.1

kwargs["A1"] = 0
kwargs["A2"] = -3
kwargs["A3"] = 0

kwargs["A4"] = 0.8

kwargs["A5"] = 0

# I use zs not as because as is a keyword in Python
kwargs["A6"] = 0.1

kwargs["mht"] = 10 ** (8)
kwargs["mst"] = 10 ** (8)


# First, let's fit separately
all_index = ["M11.0", "M11.5", "M12.0", "M12.5", "M13.0"]

model = MCMC(kwargs=kwargs)

#### Index
model.fitting_separately(index=all_index[3])

counter = 0


def lnlike(theta, x, y):
    global counter

    global para_chi_sca

    # introduce our model here.

    # Theta is a tuple
    f0, A1, A2, A3, A4,A5, A6 = theta

    kwargs["f0"] = f0
    kwargs["A1"] = A1
    kwargs["A2"] = A2
    kwargs["A3"] = A3

    kwargs["A4"] = A4
    kwargs["A5"] = A5


    kwargs["A6"] = A6

    # Only fit 4 parameters for now
    # kwargs["mht"] = mht
    # kwargs["mst"] = mst

    model.update_kwargs(kwargs=kwargs)

    # Let's change the model from simple linear to our complex model.

    # y_model is the value from our method


    chi = model.calculate_M_stellar()

    # if you want to fit simultaneously. Use this

    yerr = 0.3

    inv_sigma2 = 1.0 / (yerr ** 2)

    counter = counter + 1
    print("Doing %d" % counter)
    print(chi, f0, A1, A2, A3, A4, A5,A6)

    # Here things become simpler because we have a constant y_err
    # return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

    return -0.5 * (chi * inv_sigma2 - np.sum(np.log(inv_sigma2)))


# Set the range of the parameters
def lnprior(theta):
    f0, A1, A2, A3, A4,A5, A6 = theta

    # Let A4<0!!
    if 0 < f0 < 1 and -1 < A1 < 1 and -5 < A2 < -0.1 and -1.5 < A3 < 1.5 and -1 < A4 < 2 and -1.5 < A5 < 1.5and 0.001 < A6 < 1:
        return 0.0
    return -np.inf


# The final ln posterior probability is the sum of the ln prior and ln likelihood

def lnprob(theta, x, y):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y)


# From scipy fit:


## This part is for scipy opt
## We need this because we need them as initial values!!
## That's what you discussed with Jeremy.


# nll = lambda *args: -lnlike(*args)
nll = lambda *args: -lnprob(*args)
print("doing scipy opt")
# Only fit f0,A1,A2,A3,A4,zs
# I adjust the likelihood function

# x is a and y is logm: But it doesn't matter since lnlike has nothing to do with x and y.
x = 0
y = 0

result = op.minimize(nll, [kwargs["f0"], kwargs["A1"], kwargs["A2"], kwargs["A3"], kwargs["A4"], kwargs["A5"], kwargs["A6"]],
                     args=(x, y))

print("finish doing scipy opt")

print("result x")
print(result["x"])
f0, A1, A2, A3, A4,A5, A6 = result["x"]

kwargs["f0"] = f0
kwargs["A1"] = A1
kwargs["A2"] = A2
kwargs["A3"] = A3

kwargs["A4"] = A4
kwargs["A5"] = A5
kwargs["A6"] = A6

print("final")
print(f0, A1, A2, A3, A4,A5, A6)

save_path = "Fast_Separately_best_fit_M11.0_v2.pkl"
save_path = save_path.replace("M11.0", str(model.index))

output = open(save_path, 'wb')
pickle.dump(result["x"], output)
output.close()



# Then, let's do some MCMC:


print("doing emcee")



# The final ln posterior probability is the sum of the ln prior and ln likelihood

def lnprob(theta, x, y):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf, str(lnlike(theta, x, y)/(-0.5))
    return lp + lnlike(theta, x, y), str(lnlike(theta, x, y)/(-0.5))


# Define the initial condition of the MCMC chain: Position/initial values

ndim, nwalkers = 7, 100



kwargs["MCMC"] = 1

# Or you can replace result["x"] with your initial values
# pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
pos = [
    [kwargs["f0"], kwargs["A1"], kwargs["A2"], kwargs["A3"], kwargs["A4"], kwargs["A5"], kwargs["A6"]] + 1e-3 * np.random.randn(
        ndim) for i in range(nwalkers)]

# Set up the MCMC chain
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y), threads=8)

print("running MCMC")
sampler.run_mcmc(pos, 2000)

# Now we have an array with dimension of 100*500*3: 100 walkers, 500 steps and 3 parameters

result = sampler.chain
blobs = sampler.blobs

blobs = np.array(blobs)

print(result)

# save it:

save_path = save_path.replace("best_fit_", "EMCEE")

output = open(save_path, 'wb')
pickle.dump(result, output)
output.close()

save_path = save_path.replace("EMCEE", "EMCEE_chi_")
output = open(save_path, 'wb')
pickle.dump(blobs, output)
output.close()

