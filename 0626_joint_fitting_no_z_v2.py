import numpy as np

import math

import pickle
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

    # f_con for me

    ##### Quenching: Ignore it for now

    def f_q(self):

        return 1

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

        """
        At = A4 * log10(Mh) + A5
        As = A6
        
        
        At = A3*(log10(Mh))**3+A4 * (log10(Mh))**2 + A5* (log10(Mh))+A6

        # fix a_s
        As = 0.147



        """

        At = A3*(log10(Mh))**2+A4 * (log10(Mh)) + A5

        # fix a_s
        As = A6



        return f0 * ((Mh / mht)**A1)*(log10((Mh / mht)))**A2  * exp(-(At - a) ** 2 / As)

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

        self.median_Mh = np.nanmedian(median_fusion[:, 1])

        self.a_target = a_target

        # Let's do it

        a = a_target

        median_Mh = median_fusion[:, 1]

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

        fusion = np.c_[a, median_Ms, median_Mh, f_con_array]

        # No scatter for now:


        log_median_Ms = [log10(x) for x in median_Ms]

        log_Behroozi_Ms = [log10(x) for x in Behroozi_fusion[:, 1]]

        log_median_Ms = np.array(log_median_Ms)
        log_Behroozi_Ms = np.array(log_Behroozi_Ms)

        chi = np.sum((log_median_Ms[a_index_us] - log_Behroozi_Ms[a_index_Behroozi]) ** 2 * a_weight_factor)

        self.chi = chi

        return chi


    def calculate_scatter0(self):


        kwargs = self.kwargs

        index = str(self.index)

        # Read_Behroozi's Data

        B_path = "Behroozi_revised_M11.0.pkl"

        B_path = B_path.replace("M11.0", index)

        data_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/M10.0/"

        data_path = data_path.replace("M10.0", index)

        pkl_file = open(B_path, 'rb')
        Behroozi_fusion = pickle.load(pkl_file)
        pkl_file.close()

        # Only use 100 halos sometimes. Make it quicker
        n_halo = 100

        # result all is z + Mass_stellar

        result_all = np.array([[10, 10, 10, 10]])

        for halo in range(1, n_halo):

            file_i = "output_fulltree" + str(halo) + ".txt"

            result_i = np.loadtxt(data_path + file_i)

            # calculate M_stellar

            z = result_i[:, 5]

            a = 1 / (1 + z)

            a = np.array(a)

            a = a[::-1]

            M_h = result_i[:, 1]

            M_h = np.array(M_h)

            M_h = M_h[::-1]

            self.median_Mh = np.nanmedian(M_h)

            delta_Mh = [M_h[0]]

            Ms_now = 1

            Ms = [Ms_now]

            f_con_now = self.f_con_Jason(a[0], M_h[0], Ms_now)

            f_con_array = [f_con_now]

            for j in range(1, len(M_h)):
                delta_Mh_j = M_h[j] - M_h[j - 1]
                delta_Mh.append(delta_Mh_j)

                # calculate M_stellar
                f_con_array_j = self.f_con_Jason(a[j], M_h[j], Ms_now)

                f_con_array.append(f_con_array_j)

                delta_Ms = f_con_array_j * fb * delta_Mh_j * self.f_q()

                Ms_now = Ms_now + delta_Ms

                Ms.append(Ms_now)

            delta_Mh = np.array(delta_Mh)

            Ms = np.array(Ms)

            # Plot us

            # plt.plot(a, [log10(x) for x in Ms], "k", alpha=alpha, linewidth=1)


            f_con_array = np.array(f_con_array)

            # Now we already have small a to big a:

            # Our fusion has a+ Ms, Mh f_con_array

            fusion = np.c_[a, Ms, M_h, f_con_array]
            # plot

            # plt.plot(1/(1+fusion[:,0]),[log10(x) for x in stellar_mass],"k",alpha=alpha,linewidth=1)

            try:
                result_all = np.vstack((result_all, fusion))

            except:
                none = 1

        # calculate median

        # result_all : z + Mh
        result_all = np.array(result_all)

        a_all = result_all[:, 0]
        Ms_all = result_all[:, 1]

        median_Ms = []

        scatter = []

        for ii in range(0, len(a_target)):
            index_ii = np.where(a_all == a_target[ii])

            median_Ms.append(np.nanmedian(Ms_all[index_ii]))

            log_Ms_all = np.array([log10(x) for x in Ms_all])

            scatter_i = abs(np.percentile(log_Ms_all[index_ii], [84]) - np.percentile(log_Ms_all[index_ii], [50]))

            scatter.append(scatter_i)

        median_Ms = np.array(median_Ms)
        scatter = np.array(scatter)

        self.a_target = a_target
        self.median_Ms = median_Ms
        self.scatter = scatter
        self.scatter0 = scatter[-1]

        return scatter[-1]



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

kwargs["A1"] = 0.1
kwargs["A2"] = 0.1
kwargs["A3"] = 0

kwargs["A4"] = 0

kwargs["A5"] = 0.8

kwargs["A6"] = 0.147

kwargs["mht"] = 10 ** (8)
kwargs["mst"] = 10 ** (8)

# First, let's fit separately
all_index = ["M11.0", "M11.5", "M12.0", "M12.5", "M13.0"]

model = MCMC(kwargs=kwargs)

### index:


model.fitting_separately(index=all_index[0])


#### Index

counter = 0


def lnlike(theta, x, y):
    global counter

    global para_chi_sca

    # introduce our model here.

    # Theta is a tuple
    f0, A1, A2, A3, A4, A5, A6 = theta

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


    chi = model.fitting_simultaneously()
    # chi = model.calculate_M_stellar()

    # if you want to fit simultaneously. Use this

    yerr = 0.3

    inv_sigma2 = 1.0 / (yerr ** 2)

    counter = counter + 1
    print("Doing %d" % counter)
    print(chi, f0, A1, A2, A3, A4, A5, A6)

    # Here things become simpler because we have a constant y_err
    # return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
    # return -0.5 * (chi * inv_sigma2)
    return -0.5 * (chi * inv_sigma2 - np.sum(np.log(inv_sigma2)))


# Set the range of the parameters
def lnprior(theta):
    f0, A1, A2, A3, A4, A5, A6 = theta

    # Don't simply let A4<0. It's not that simple...
    if 0.001 < f0 < 1 and -5 < A1 < 5 and -5 < A2 < 5 and -3 < A3 < 3 and -3 < A4 < 3 and -1 < A5 < 2 and -1 < A6 < 2:
    #if 0.001 < f0 < 2 and -5 < A1 < 5 and -5 < A2 < 5 and -0.01 < A3 < 0.01 and -5 < A4 < 5 and -5 < A5 < 5 and 0.001 < A6 < 3:
        return 0.0
    return -np.inf


# The final ln posterior probability is the sum of the ln prior and ln likelihood

def lnprob(theta, x, y):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y)



# No scipy for now

# Then, let's do some MCMC:

# Now only use

print("doing emcee")


# The final ln posterior probability is the sum of the ln prior and ln likelihood


# Define the initial condition of the MCMC chain: Position/initial values

def lnprob(theta, x, y):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf, str(lnlike(theta, x, y) / (-0.5))
    return lp + lnlike(theta, x, y), str(lnlike(theta, x, y) / (-0.5))


# Define the initial condition of the MCMC chain: Position/initial values

ndim, nwalkers = 7, 200

kwargs["MCMC"] = 1

# Or you can replace result["x"] with your initial values
# pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
pos = [
    [kwargs["f0"], kwargs["A1"], kwargs["A2"], kwargs["A3"], kwargs["A4"], kwargs["A5"],
     kwargs["A6"]] + 1e-3 * np.random.randn(
        ndim) for i in range(nwalkers)]

# Set up the MCMC chain
x = 1
y = 1
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y), threads=8)

print("running MCMC")
sampler.run_mcmc(pos, 4000)

# Now we have an array with dimension of 100*500*3: 100 walkers, 500 steps and 3 parameters

result = sampler.chain
blobs = sampler.blobs

blobs = np.array(blobs)

print(result)

# save it:
# .replace("EMCEE","EMCEE_"+str(model.index))

save_path = "Simultaneously_no_z_best_fit_v2.pkl"
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



chi = np.array(chi,dtype=float)
index = np.argwhere(chi == np.nanmin(chi))

length_index = len(index)
# print("length_index")
#print(length_index)

#### Choose minimum scatter0
scatter0_array = []
for mk in range(0,length_index):

    i, j = index[mk]

    f0, A1, A2, A3, A4, A5, A6 = result[i, j, :]
    print("chi for different index")
    print(f0, A1, A2, A3, A4, A5, A6)
    print(chi[i, j])

    kwargs["f0"] = f0
    kwargs["A1"] = A1
    kwargs["A2"] = A2
    kwargs["A3"] = A3
    kwargs["A4"] = A4
    kwargs["A5"] = A5

    kwargs["A6"] = A6

    model.update_kwargs(kwargs=kwargs)

    scatter0_i = model.calculate_scatter0()
    print("scatter_i")
    print(scatter0_i)
    print("chi-squared")
    print(chi[i,j])
    scatter0_array.append(scatter0_i)

scatter0_array = np.array(scatter0_array)

min_index = np.argmin(scatter0_array)


# Find the least chi-squared on:

# Simultaneously

i, j = index[min_index]

f0, A1, A2, A3, A4, A5, A6 = result[i, j, :]
print("chi for different index")
print(f0, A1, A2, A3, A4, A5, A6)
print(chi[i, j])

kwargs["f0"] = f0
kwargs["A1"] = A1
kwargs["A2"] = A2
kwargs["A3"] = A3
kwargs["A4"] = A4


kwargs["A5"] = A5
kwargs["A6"] = A6


nll = lambda *args: -lnlike(*args)

print("doing scipy opt")
# Only fit f0,A1,A2,A3
# I adjust the likelihood function

result = op.minimize(nll, [kwargs["f0"], kwargs["A1"], kwargs["A2"], kwargs["A3"], kwargs["A4"], kwargs["A5"], kwargs["A6"]],
                     args=(x, y))
print("finish doing scipy opt")

print("result x")
print(result["x"])
f0, A1, A2, A3,A4,A5,A6 = result["x"]


output = open(save_path_backup, 'wb')
pickle.dump(result["x"], output)
output.close()




