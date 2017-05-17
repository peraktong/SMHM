import numpy as np
import pickle
import math
import scipy
import matplotlib.pyplot as plt
from termcolor import colored

import pickle

import emcee
import scipy.optimize as op



ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/"

data_path = ori + "M13" + "/"

plot_path = "/Users/caojunzhi/Downloads/upload_2017.5.1_Jeremy/"


def log10(x):
    if x > 0:
        return math.log10(x)
    else:
        return -np.inf


ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/"

data_path = ori + "M13" + "/"

plot_path = "/Users/caojunzhi/Downloads/upload_2017.5.1_Jeremy/"


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

# use M13, read data


pkl_file = open('B_M13.pkl', 'rb')
B_M13 = pickle.load(pkl_file)
pkl_file.close()


class MCMC():
    def __init__(self, kwargs):
        self.kwargs = kwargs

        self.kwargs_temp = kwargs

    # Ours.


    def stellar_mass_per_interval(self, f_con, delta_Mh):

        fb = 0.15

        return f_con * fb * delta_Mh

    def f_con_z(self, z, f0):

        z0 = 1
        gamma1 = -3
        gamma2 = 0

        if z > z0:
            return f0 * ((1 + z) / (1 + z0)) ** gamma1

        else:
            return f0 * ((1 + z) / (1 + z0)) ** gamma2

    def input_array(self, z_array, delta_mh_array):

        self.z_array = z_array
        self.delta_mh_array = delta_mh_array

    ## Here comes f_q

    def f_q(self, z=None, Mh=None, Ms=None):

        # a little tricky..
        kwargs = self.kwargs_temp
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

        if ty == "redshift":

            # since we have non-integer exponent, there are complex numbers
            # abs?

            # print(abs(((z-zc)/sigmaz)**alphaz))

            #print(abs(((z-zc)/sigmaz)**alphaz))
            if z>zc:
                return 1
            else:
                return (exp((z - zc) / sigmaz)) ** alphaz




        elif ty == "halo":
            if Mh>mhc:
                return exp(((log10(mhc) - log10(Mh)) / sigmah) ** alphah)
            else:
                return 1

        elif ty == "galaxy":
            if Ms>msc:
                return exp(((log10(msc) - log10(Ms)) / sigmag) ** alphag)
            else:
                return 1

    def calculate_mass(self):
        # input array: f_con at zi, deltaMh at zi
        # return an array of M*

        ## kwargs
        kwargs = self.kwargs_temp
        f0 = self.kwargs_temp["f0"]

        z_array = self.z_array
        delta_mh_array = self.delta_mh_array

        M_stellar_interval = []

        #### Add f_q

        for i in range(0, len(z_array)):
            result = self.stellar_mass_per_interval(f_con=self.f_con_z(z_array[i], f0=f0), delta_Mh=delta_mh_array[i])

            M_stellar_interval.append(result)

        M_stellar_interval = np.array(M_stellar_interval)
        #### Add f_q


        cal_fq = np.vectorize(self.f_q)

        f_q_array_temp = np.array(cal_fq(z=z_array))
        # print(z_array)
        # print(f_q_array_temp)
        # print(f_q_array_temp.shape,M_stellar_interval.shape)

        # multiply them (f_q)

        M_stellar_interval = M_stellar_interval * f_q_array_temp

        M_stellar_interval = np.array(M_stellar_interval)

        # calculate

        M_stellar = []

        for j in range(0, len(M_stellar_interval)):
            # attention! small t small a big z: You Use Z in your calculation!
            M_stellar.append(np.sum(M_stellar_interval[j:]))

        self.M_stellar = np.array(M_stellar)
        # print(M_stellar_interval)
        # print(M_stellar)
        return M_stellar

    def calculate_median_us(self):


        pkl_file = open('B_M13.pkl', 'rb')
        B_M13 = pickle.load(pkl_file)
        pkl_file.close()

        self.B_M13 = B_M13


        # a little tricky here
        kwargs = self.kwargs_temp

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

            self.input_array(z_array=z, delta_mh_array=delta_Mh)

            # The calculate mass module is different:
            M_stellar = self.calculate_mass()

            fusion = np.c_[z, M_stellar]
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
        target = z_target

        target = list(set(target))

        M_stellar_all = []

        for i in range(0, len(target)):
            # print("Doing %.2f %%" % (i / len(target) * 100))

            index = np.where(abs(result_all[:, 0] - target[i]) < 0.01)
            index = np.array(index)
            index = index.ravel()

            # Mass
            array = result_all[:, 1][index]

            # in dex

            array_dex = np.array([log10(x) for x in array])

            # print(target[i])
            # print(array)
            M_stellar_all.append(np.nanmedian(array))

        M_stellar_all = np.array(M_stellar_all)

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

        target = self.B_M13[:, 0]
        index_us = []

        for mk in range(0, len(target)):
            index = np.where(a_us == target[mk])

            index_us.append(index)

        index_us = np.array(index_us, dtype=int)

        # This is the part where we calculate the chi!

        a_us_chi = a_us[index_us]
        m_us_chi = m_us[index_us]

        # give a chi-squared


        m_Behroozi = np.array(self.B_M13[:, 1]).ravel()
        m_us_chi = np.array(m_us_chi).ravel()

        self.m_Behroozi = m_Behroozi
        self.m_us_chi = m_us_chi

        self.chi_temp = np.sum((m_Behroozi - m_us_chi) ** 2)

        return self.chi_temp

    def do_MCMC(self, step):

        pkl_file = open('B_M13.pkl', 'rb')
        B_M13 = pickle.load(pkl_file)
        pkl_file.close()

        self.B_M13 = B_M13

        chi = self.calculate_median_us()
        temp = np.array([])

        for i in range(0, step):

            print("doing %d" % i)

            """

            # Tricky : Let's do it periodically

            residual = i%4

            if residual==0:
                self.kwargs_temp["f0"] = np.random.normal(self.kwargs["f0"], 1)
            elif residual==1:
                self.kwargs_temp["zc"] = np.random.normal(self.kwargs["zc"], 1)

            elif residual ==2:

                self.kwargs_temp["sigmaz"] = np.random.normal(self.kwargs["sigmaz"], 1)
            elif residual ==3:
                self.kwargs_temp["alphaz"] = np.random.normal(self.kwargs["alphaz"], 1)

            """

            self.kwargs_temp["f0"] = np.random.normal(self.kwargs["f0"], 1)
            self.kwargs_temp["zc"] = np.random.normal(self.kwargs["zc"], 1)
            self.kwargs_temp["sigmaz"] = np.random.normal(self.kwargs["sigmaz"], 1)
            self.kwargs_temp["alphaz"] = np.random.normal(self.kwargs["alphaz"], 1)

            # calculate chi-squared_temp

            chi_temp = self.calculate_median_us()

            # if (ratio>r)&(np.isfinite(chi-chi_temp))&(self.kwargs_temp["f0"]>0)&(self.kwargs_temp["f0"]<2)&(self.kwargs_temp["zc"]>0)&(self.kwargs_temp["zc"]<6)&(self.kwargs_temp["sigmaz"]>0)&(self.kwargs_temp["sigmaz"]<5)&(self.kwargs_temp["alphaz"]>0)&(self.kwargs_temp["alphaz"]<5):

            # from model scatter we derive sig
            sig = 0.215
            ratio = exp((-chi_temp + chi) / (2 * sig ** 2))

            r = np.random.rand()
            print(ratio, chi, chi_temp)

            # Let's save the variables for all steps:
            temp = np.append(temp, self.kwargs)

            output = open('redshift_result_all_M13.pkl', 'wb')
            pickle.dump(temp, output)
            output.close()

            #

            if (ratio > r) & (np.isfinite(chi - chi_temp)):
                self.kwargs = self.kwargs_temp
                chi = chi_temp

                print(colored("better", "red"))
                print((np.nanmean(chi)) ** 0.5)
                print(self.kwargs)
                # Let's save it:

                output = open('redshift_result_M13.pkl', 'wb')
                pickle.dump(self.kwargs, output)
                output.close()



            else:
                self.kwargs_temp = self.kwargs

            ## The initial values of these parameters are in kwargs


            kwargs = {"f0": None, "zc": None, "sigmaz": None, "alphaz": None, "mhc": None, "msc": None, "sigmah": None,
                      "sigmag": None, "alphah": None, "alphag": None}

            # Let's do the redshift first:

            kwargs["ty"] = "redshift"

            kwargs["f0"] = 0.34
            kwargs["zc"] = 3.14
            kwargs["sigmaz"] = 1.1
            kwargs["alphaz"] = 2.14

            print(kwargs)




# Now let's come to the main topic:


## The initial values of these parameters are in kwargs


kwargs = {"f0":None,"zc":None, "sigmaz": None ,"alphaz":None,"mhc":None,"msc":None,"sigmah":None,"sigmag":None,"alphah":None,"alphag":None}


# Let's do the redshift first:

kwargs["ty"] = "redshift"

# initial values of kwargs:

kwargs["f0"] = 0.34
kwargs["zc"] = 3.14
kwargs["sigmaz"] = 1.1
kwargs["alphaz"] = 2.14

# print(kwargs)


# define x, y and y_err

# x is the selected redshift from Behroozi 2013, which is 26-27 long.

# y is the data from Behroozi 2013, which is about 26-27 long. You can read them from B_M13.pkl

# y_err is the scatter of the model, which is set manually


## Thus, you need an effective way to convert 101(x) in 26(y):

## Read data:


# use M13, read data


pkl_file = open('B_M13.pkl', 'rb')
B_M13 = pickle.load(pkl_file)
pkl_file.close()


# x is a, y is logm:

x = B_M13[:,0]
y = B_M13[:,1]


# Let's define ln PDF functions:

# Set y_err =  the scatter of the model

yerr = 0.4
counter = 0

def lnlike(theta, x, y):

    global  counter

    # introduce our model here.

    # Theta is a tuple
    f0,zc,sigmaz,alphaz = theta

    kwargs["f0"] = f0
    kwargs["zc"] = zc
    kwargs["sigmaz"] = sigmaz
    kwargs["alphaz"] = alphaz

    counter = counter+1
    print("Doing %d"%counter)



    model = MCMC(kwargs=kwargs)

    # Let's change the model from simple linear to our complex model.

    # y_model is the value from our method
    chi = model.calculate_median_us()

    # print(kwargs)
    print("check chi")
    print(chi)


    y_model = model.m_us_chi
    inv_sigma2 = 1.0/(yerr**2)


    # Here things become simpler because we have a constant y_err
    #return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))


    return -0.5 * (np.sum((y_model - y) ** 2 * inv_sigma2 - np.log(inv_sigma2)))


# Set the range of the parameters
def lnprior(theta):

    f0, zc, sigmaz, alphaz = theta

    if 0.1 < f0 < 0.5 and 2.5 < zc < 5 and 0.8 < sigmaz < 3 and 2 < alphaz < 5.5:
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
result = op.minimize(nll, [kwargs["f0"],kwargs["zc"],kwargs["sigmaz"],kwargs["alphaz"]], args=(x, y))
print("finish doing scipy opt")

print("result x")
print(result["x"])
f0,zc,sigmaz,alphaz = result["x"]

kwargs["f0"] = f0
kwargs["zc"] = zc
kwargs["sigmaz"] = sigmaz
kwargs["alphaz"] = alphaz

output = open('EMCEE_redshift_M13_scipy.pkl', 'wb')
pickle.dump(result["x"], output)
output.close()


output = open('kwargs_redshift.pkl', 'wb')
pickle.dump(kwargs, output)
output.close()



###

print("doing emcee")


# Define the initial condition of the MCMC chain: Position/initial values

ndim, nwalkers = 4, 40

# Or you can replace result["x"] with your initial values
# pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
pos = [
    [kwargs["f0"],kwargs["zc"],kwargs["sigmaz"],kwargs["alphaz"]]+ 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# Set up the MCMC chain
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y),threads=8)

"""


f = open("EMCEE_redshift_M13.dat", "w")
f.close()

for result in sampler.sample(pos, iterations=500, storechain=False):
    position = result[0]
    f = open("EMCEE_redshift_M13.dat", "a")
    for k in range(position.shape[0]):
        f.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
    f.close()


"""
print("running MCMC")
sampler.run_mcmc(pos, 200)

# Now we have an array with dimension of 100*500*3: 100 walkers, 500 steps and 3 parameters

result = sampler.chain
print(result)

# save it:


output = open('EMCEE_redshift_M13.pkl', 'wb')
pickle.dump(result, output)
output.close()





