import numpy as np
import pickle
import math
import scipy
import matplotlib.pyplot as plt
import matplotlib
from termcolor import colored
import corner
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


class diagnostic():

    def __init__(self,chain,name):

        self.chain = chain
        self.name = name



    def walker_plot(self):

        name =self.name

        chain = self.chain
        n_parameter =len(chain[0,0,:])

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 15}

        matplotlib.rc('font', **font)



        for i in range(0,n_parameter):


            plt.subplot(n_parameter, 1, i+1)

            alpha = .3

            for j in range(0,len(chain[:,0,0])):

                plt.plot(chain[j,:,i],"k",linewidth=0.7, alpha=alpha)


            plt.ylabel("%s" % name[i], fontsize=20)

        plt.xlabel("Step")
        plt.legend()

        plt.suptitle("The values of walkers vs step", fontsize=20)

        # share x
        plt.subplots_adjust(hspace=.5)

        # plt.show()
        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(16.5, 11.5)

        # select some plots:


        save_path = "/Users/caojunzhi/Downloads/upload_20170519_Jeremy/" + "redshift_walker" + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()

    def corner_plot(self):
        name = self.name

        # Only choose resultf form steps after 50
        samples = self.chain.reshape((-1, 4))

        samples[:, 2] = np.exp(samples[:, 2])


        f0, A1, A2, A3 = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                     zip(*np.percentile(samples, [16, 50, 84],
                                                        axis=0)))

        print(f0,A1,A2,A3)

        #fig = corner.corner(samples, labels=name,truths=[f0,A1,A2,A3,mht,mst])

        fig = corner.corner(samples, labels=name)


        fig.savefig("/Users/caojunzhi/Downloads/upload_20170519_Jeremy/"+"redshift_triangle.png")



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

    def save_plot(self):


        mass_name = "M12"

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 20}

        matplotlib.rc('font', **font)


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

            slope = 1

            plt.plot(1 / (1 + z[::-1]), [slope * log10(x) for x in M_stellar], "k", alpha=trans, linewidth=1)

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

        self.chi_temp = np.sum((m_Behroozi - m_us_chi) ** 2)



        slope = 1



        plt.plot(a_us, slope * m_us, "b", linewidth=3, label="$Median \quad value$")



        std = (np.sum((m_Behroozi - slope * m_us_chi) ** 2) / len(m_Behroozi))**0.5

        # print(residual.shape)

        scatter = (np.nansum(residual**2)/len(residual))**0.5
        scatter_0 = residual[-1]

        plt.plot([], [], "ko", label="$\sigma\quad log[M_{*}]_{z=0} = %.3f$" % (scatter_0))

        plt.plot([], [], "ko", label="$STD = %.3f$" % (std))


        # Plot Behroozi 2013

        a_Behroozi = B_M12[:,0]
        m_Behroozi = B_M12[:,1]

        plt.plot(a_Behroozi, m_Behroozi, "ro", linewidth=3, label="$Behroozi \quad 2013$")

        plt.plot([], [], "k", label="$Subsample \quad of \quad Haloes$", alpha=trans)
        plt.plot([], [], "ko", label="$log[M_h] \quad = %.f (dex)$" % (np.nanmax([log10(x) for x in M_h])))



        plt.legend()

        # save it:
        plt.suptitle("$log[M_*]\quad vs\quad scale \quad factor\quad new \quad f_{con}$")

        plt.xlabel("$Scale\quad factor\quad a$")
        plt.ylabel("$log[M_*]\quad (dex)$", fontsize=20)

        axes = plt.gca()
        axes.set_xlim([0.2, 1])
        #axes.set_ylim([9, 12])

        # share x axis



        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(18.5, 12.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170519_Jeremy/" + "brain_storm_emcee_" + str(mass_name) + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


# Now let's come to the main topic:


## The initial values of these parameters are in kwargs

# Let's look into f_con and have a new f_con, which are related to z, Mh and M_*. The parameters we use i A1, A2 and A3



kwargs = {"A1": None, "A2": None, "A3": None, "f0": None, "zc": None, "sigmaz": None, "alphaz": None, "mhc": None,
          "msc": None, "sigmah": None, "sigmag": None, "alphah": None, "alphag": None}

# Let's do the redshift first:

kwargs["ty"] = "redshift"

# initial values of kwargs:
# Now we have f0, A1, A2 and A3

# log f_con = log(f0) + A1*log(M_halo[z]) + A2*log(M_star[z]) + A3*log(1+z)
#


# read parameters for f0, A1 A2 and A3


pkl_file = open('EMCEE_new_fcon_M12_scipy_log.pkl', 'rb')
scipy_result = pickle.load(pkl_file)
pkl_file.close()

f0,A1,A2,A3 = scipy_result

print(scipy_result)



kwargs["f0"] = f0

kwargs["A1"] = A1
kwargs["A2"] = A2
kwargs["A3"] = A3

kwargs["mht"] = 10e7
kwargs["mst"] = 10e7

model = MCMC(kwargs=kwargs)
model.save_plot()




# Some diagnostic plots:

pkl_file = open('EMCEE_new_fcon_M12_log.pkl', 'rb')
emcee = pickle.load(pkl_file)
pkl_file.close()


name = ["f0","A1","A2","A3"]

model = diagnostic(chain=emcee[:,50:,:],name=name)

model.walker_plot()
model.corner_plot()




