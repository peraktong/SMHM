import numpy as np
import pickle
import math
import scipy
import matplotlib.pyplot as plt
from termcolor import colored

import pickle

import emcee
import scipy.optimize as op


def log10(x):
    if x > 0:
        return math.log10(x)
    else:
        return -np.inf


ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/"

data_path = ori + "M12" + "/"


def log10(x):
    if x > 0:
        return math.log10(x)
    else:
        return -np.inf


"""

    if x>0:
        return math.log10(x)
    else:
        return float("nan")

"""

log10 = np.vectorize(log10)


def exp(x):
    try:
        return math.exp(x)
    except OverflowError:
        if x > 0:
            return np.inf
        else:
            return 0


"""

    try:
        return math.exp(x)
    except OverflowError:
        return float("nan")


"""

exp = np.vectorize(exp)

# use M12, read data


pkl_file = open('B_M12.pkl', 'rb')
B_M12 = pickle.load(pkl_file)
pkl_file.close()


class MCMC():
    def __init__(self, kwargs):
        self.kwargs = kwargs

    # Ours.


    def stellar_mass_per_interval(self, f_con, delta_Mh):

        fb = 0.15

        return f_con * fb * delta_Mh

    # Now we have a different f_con:

    def f_con_z(self, z, Mh, Ms):

        kwargs = self.kwargs

        f0 = kwargs["f0"]

        A1 = kwargs["A1"]
        A2 = kwargs["A2"]
        A3 = kwargs["A3"]

        z0 = 1
        gamma1 = -3
        gamma2 = 0

        """


        if z > z0:

            print("f_con")
            print(f0 * ((1 + z) / (1 + z0)) ** gamma1)
            return f0 * ((1 + z) / (1 + z0)) ** gamma1

        else:
            return f0 * ((1 + z) / (1 + z0)) ** gamma2




        """
        # print("ZA123")
        # print(f0,A1,A2,A3,f0*(Mh**A1)*(Ms**A2)*((1.+z)**A3))


        mht = kwargs["mht"]
        mst = kwargs["mst"]

        if z > z0:
            return f0 * ((log10(Mh / mht)) ** A1) * ((log10(Ms / mst)) ** A2) * ((1. + z) ** A3)
        else:
            return f0 * ((log10(Mh / mht)) ** A1) * ((log10(Ms / mst)) ** A2) * ((1. + z) ** A3)

    ## Here comes f_q

    def f_q(self, z=None, Mh=None, Ms=None):

        # a little tricky..
        kwargs = self.kwargs
        ty = kwargs["ty"]

        zc = kwargs["zc"]
        sigmaz = kwargs["sigmaz"]
        alphaz = kwargs["alphaz"]

        mhc = kwargs["mhc"]
        sigmah = kwargs["sigmah"]
        alphah = kwargs["alphah"]

        msc = kwargs["msc"]
        sigmag = kwargs["sigmag"]
        alphag = kwargs["alphag"]

        # Let's ignore the quenching function for now.

        if ty == "redshift":

            # since we have non-integer exponent, there are complex numbers
            # abs?

            # print(abs(((z-zc)/sigmaz)**alphaz))

            # print(abs(((z-zc)/sigmaz)**alphaz))

            return 1

            """

            if z>zc:
                return 1
            else:

                return (exp((z - zc) / sigmaz)) ** alphaz


            """






        elif ty == "halo":
            if Mh > mhc:
                return exp(((log10(mhc) - log10(Mh)) / sigmah) ** alphah)
            else:
                return 1

        elif ty == "galaxy":
            if Ms > msc:
                return exp(((log10(msc) - log10(Ms)) / sigmag) ** alphag)
            else:
                return 1

    def calculate_mass(self):
        # input array: f_con at zi, deltaMh at zi
        # return an array of M*

        ## kwargs
        kwargs = self.kwargs
        f0 = self.kwargs["f0"]

        z_array = self.z_array
        delta_mh_array = self.delta_mh_array

        # Now M_stellar_interval Mh and M* are just numbers, not arrays.




        # Revert the array

        z_array_inverse = z_array[::-1]
        delta_mh_array_inverse = delta_mh_array[::-1]
        Mh_all_reverse = self.Mh_all[::-1]

        #### Now we need to input Mh and M* for f_con

        # initial values



        # Time start at big z, small

        M_stellar_interval = []

        # intial is from the f_con in Jeremy's  (2)
        Ms_now = Mh_all_reverse[0] * 0.1 * 0.0045

        for i in range(0, len(z_array_inverse)):
            Mh = Mh_all_reverse[i]

            M_stellar_interval_i = self.stellar_mass_per_interval(
                f_con=self.f_con_z(z=z_array_inverse[i], Mh=Mh, Ms=Ms_now), delta_Mh=delta_mh_array_inverse[i])

            Ms_now += M_stellar_interval_i

            M_stellar_interval.append(M_stellar_interval_i)

        M_stellar_interval = np.array(M_stellar_interval)
        #### Add f_q


        cal_fq = np.vectorize(self.f_q)

        f_q_array_temp = np.array(cal_fq(z=z_array_inverse))
        # print(z_array)
        # print(f_q_array_temp)
        # print(f_q_array_temp.shape,M_stellar_interval.shape)

        # multiply them (f_q)

        M_stellar_interval = M_stellar_interval * f_q_array_temp

        M_stellar_interval = np.array(M_stellar_interval)

        # calculate

        M_stellar = []

        # Now the m_stellar_interval is inversed!!

        for j in range(0, len(M_stellar_interval)):
            # attention! small t small a big z: You Use Z in your calculation!
            # And now they are inversed!

            M_stellar.append(np.sum(M_stellar_interval[0:j]))

        self.M_stellar = np.array(M_stellar)
        # print(M_stellar_interval)
        # print(M_stellar)
        return M_stellar

    def calculate_median_us(self):

        # Let's call it B_13 for now
        pkl_file = open('B_M12.pkl', 'rb')
        B_M12 = pickle.load(pkl_file)
        pkl_file.close()

        self.B_M12 = B_M12

        # a little tricky here
        kwargs = self.kwargs

        trans = 0.05

        # Only use 100 halos for now. Make it quicker
        n_halo = 300

        # result all is z + Mass_stellar

        result_all = np.array([[0, 0]])

        # What we need is 1 and 5
        # 1 Halo Mass
        # 5 red shift z

        for halo in range(1, n_halo):

            file_i = "output_fulltree" + str(halo) + ".txt"

            result_i = np.loadtxt(data_path + file_i)

            # calculate M_stellar

            M_h = result_i[:, 1]

            self.Mh_all = M_h
            # print(M_h)

            delta_Mh = []

            for j in range(0, len(M_h) - 1):
                delta_Mh.append(M_h[j] - M_h[j + 1])

            delta_Mh.append(M_h[len(M_h) - 1])

            delta_Mh = np.array(delta_Mh)

            # calculate M_stellar

            z = result_i[:, 5]

            if halo == 1:
                z_target = z

            else:
                none = 1

            # read parameters

            self.z_array = z
            self.delta_mh_array = delta_Mh

            # The calculate mass module is different:
            M_stellar = self.calculate_mass()

            # Now they are from big z to small z:

            fusion = np.c_[z[::-1], M_stellar]
            # plot

            # plt.plot(1/(1+fusion[:,0]),[log10(x) for x in stellar_mass],"k",alpha=alpha,linewidth=1)

            try:
                result_all = np.vstack((result_all, fusion))

            except:
                none = 1

        # add median

        # result_all : z + Mh
        result_all = np.array(result_all)

        # calculate median:


        # result_all : z + Mh
        result_all = np.array(result_all)

        # print("result shape")
        # print(result_all.shape)

        # since all of them are from Main branch, it's okay to do this
        target = z_target[::-1]

        target = list(set(target))

        M_stellar_all = []

        index_z0 = np.where(abs(result_all[:, 0] - target[-1]) < 0.01)

        index_z0 = np.array(index_z0).ravel()

        M_stellar_z0 = result_all[:, 1][index_z0]

        self.M_stellar_z0 = M_stellar_z0

        residual = []

        for i in range(0, len(target)):
            # print("Doing %.2f %%" % (i / len(target) * 100))

            index = np.where(abs(result_all[:, 0] - target[i]) < 0.01)
            index = np.array(index)
            index = index.ravel()

            # Mass
            array = result_all[:, 1][index]

            # in dex

            array_dex = np.array([log10(x) for x in array])


            scatter_i = abs(np.percentile(array_dex,[84])-np.percentile(array_dex,[50]))
            residual.append(scatter_i)


            # print(target[i])
            # print(array)
            M_stellar_all.append(np.nanmedian(array))

        M_stellar_all = np.array(M_stellar_all)

        residual = np.array(residual).ravel()

        scatter = (np.nansum(residual ** 2) / len(residual)) ** 0.5

        # We limit scatter0 to be smaller than 0.2
        scatter_0 = residual[-1]

        # print(M_stellar_all)

        # Now we have z_target + M_stellar_all


        target = np.array(target)
        # print(target)
        # print(M_stellar_all)

        a_us = 1 / (1 + target)

        m_us = [log10(x) for x in M_stellar_all]
        m_us = np.array(m_us, dtype=float)
        # print(m_us)

        # sort:

        index = np.argsort(a_us)

        a_us = a_us[index]
        m_us = m_us[index]

        a_us = np.array(a_us, dtype=float)
        m_us = np.array(m_us, dtype=float)

        target = self.B_M12[:, 0]
        index_us = []

        for mk in range(0, len(target)):
            index = np.where(a_us == target[mk])

            index_us.append(index)

        index_us = np.array(index_us, dtype=int)

        # This is the part where we calculate the chi!

        a_us_chi = a_us[index_us]
        m_us_chi = m_us[index_us]

        # give a chi-squared


        m_Behroozi = np.array(self.B_M12[:, 1]).ravel()
        m_us_chi = np.array(m_us_chi).ravel()

        self.a_us = a_us
        self.m_us = m_us
        self.a_Behroozi = np.array(self.B_M12[:, 0]).ravel()

        self.m_Behroozi = m_Behroozi
        self.m_us_chi = m_us_chi


        # limit scatter0<0.2

        if (scatter_0 > 0.2)&(np.sum((m_Behroozi - m_us_chi) ** 2)<1):
            # &(self.kwargs["MCMC"]!=1)
            # self.chi_temp = np.inf

            # I think we can make it tricky here:
            # if chi_temp <1, set it to inf. Else set to normal!!
            # When running MCMC, you'd better cancel this term!!

            self.chi_temp = np.sum((m_Behroozi - m_us_chi) ** 2)

        else:

            self.chi_temp = np.sum((m_Behroozi - m_us_chi) ** 2)
        print("check chi")
        print(self.chi_temp)



        return self.chi_temp


# Now let's come to the main topic:


## The initial values of these parameters are in kwargs

# Let's look into f_con and have a new f_con, which are related to z, Mh and M_*. The parameters we use i A1, A2 and A3


kwargs = {"MCMC":None,"A1": None, "A2": None, "A3": None, "f0": None, "mht": None, "mst": None, "zc": None, "sigmaz": None,
          "alphaz": None, "mhc": None, "msc": None, "sigmah": None, "sigmag": None, "alphah": None, "alphag": None}

# Let's do the redshift first:

kwargs["ty"] = "redshift"

# initial values of kwargs:
# Now we have f0, A1, A2 and A3

# log f_con = log(f0) + A1*log(M_halo[z]) + A2*log(M_star[z]) + A3*log(1+z)
#


kwargs["f0"] = 0.1

kwargs["A1"] = 0
kwargs["A2"] = 0
kwargs["A3"] = -3

kwargs["mht"] = 10 ** (7)
kwargs["mst"] = 10 ** (7)

# Threshold




model = MCMC(kwargs=kwargs)
model.calculate_median_us()

# print(kwargs)


# define x, y and y_err

# x is the selected redshift from Behroozi 2013, which is 26-27 long.

# y is the data from Behroozi 2013, which is about 26-27 long. You can read them from B_M12.pkl

# y_err is the scatter of the model, which is set manually


## Thus, you need an effective way to convert 101(x) in 26(y):

## Read data:


# use M12, read data


pkl_file = open('B_M12.pkl', 'rb')
B_M12 = pickle.load(pkl_file)
pkl_file.close()

# x is a, y is logm:

x = B_M12[:, 0]
y = B_M12[:, 1]

# Let's define ln PDF functions:

# Set y_err =  the scatter of the model

yerr = 0.3
counter = 0


def lnlike(theta, x, y):
    global counter

    # introduce our model here.

    # Theta is a tuple
    f0, A1, A2, A3 = theta

    kwargs["f0"] = f0
    kwargs["A1"] = A1
    kwargs["A2"] = A2
    kwargs["A3"] = A3

    # Only fit 4 parameters for now
    #kwargs["mht"] = mht
    #kwargs["mst"] = mst

    counter = counter + 1
    print("Doing %d" % counter)

    model = MCMC(kwargs=kwargs)

    # Let's change the model from simple linear to our complex model.

    # y_model is the value from our method
    chi = model.calculate_median_us()

    # print(kwargs)


    y_model = model.m_us_chi
    inv_sigma2 = 1.0 / (yerr ** 2)

    # Here things become simpler because we have a constant y_err
    # return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))


    return -0.5 * (chi * inv_sigma2 - np.sum(np.log(inv_sigma2)))


# Set the range of the parameters
def lnprior(theta):
    f0, A1, A2, A3 = theta

    if 0 < f0 < 3 and -5 < A1 < 5 and -5 < A2 < 5 and -15 < A3 < 15 :
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


nll = lambda *args: -lnlike(*args)

print("doing scipy opt")
# Only fit f0,A1,A2,A3
# I adjust the likelihood function

result = op.minimize(nll, [kwargs["f0"], kwargs["A1"], kwargs["A2"], kwargs["A3"]],
                     args=(x, y))
print("finish doing scipy opt")

print("result x")
print(result["x"])
f0, A1, A2, A3 = result["x"]

kwargs["f0"] = f0
kwargs["A1"] = A1
kwargs["A2"] = A2
kwargs["A3"] = A3

print("final")
print(f0, A1, A2, A3)

output = open('EMCEE_new_fcon_M12_scipy_log.pkl', 'wb')
pickle.dump(result["x"], output)
output.close()

output = open('kwargs_new_fcon_log.pkl', 'wb')
pickle.dump(kwargs, output)
output.close()

###



print("doing emcee")

# Define the initial condition of the MCMC chain: Position/initial values

ndim, nwalkers = 4, 40

kwargs["MCMC"] = 1

# Or you can replace result["x"] with your initial values
# pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
pos = [
    [kwargs["f0"], kwargs["A1"], kwargs["A2"], kwargs["A3"]] + 1e-4 * np.random.randn(
        ndim) for i in range(nwalkers)]

# Set up the MCMC chain
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y), threads=8)

print("running MCMC")
sampler.run_mcmc(pos, 300)

# Now we have an array with dimension of 100*500*3: 100 walkers, 500 steps and 3 parameters

result = sampler.chain
print(result)

# save it:


output = open('EMCEE_new_fcon_M12_log.pkl', 'wb')
pickle.dump(result, output)
output.close()







