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
import corner


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

plot_path = "/Users/caojunzhi/Downloads/upload_201706_Jeremy/"


class MCMC():
    def __init__(self, kwargs):
        self.kwargs = kwargs

    def update_kwargs(self, kwargs):
        self.kwargs = kwargs

    def stellar_mass_per_interval(self, f_con, delta_Mh):

        fb = 0.15

        return f_con * fb * delta_Mh

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

        At = A3 * (log10(Mh)) ** 2 + A4 * (log10(Mh)) + A5

        # fix a_s
        As = A6

        return f0 * ((Mh / mht) ** A1) * (log10((Mh / mht))) ** A2 * exp(-(At - a) ** 2 / As)

    def f_q(self):

        return 1

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




    def Plot_M_stellar(self):

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 15}

        matplotlib.rc('font', **font)

        alpha = 0.05

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
        n_halo = 300

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

            # Attention!

            log_Ms_all = np.array([log10(x) for x in Ms_all])

            scatter_i = abs(np.percentile(log_Ms_all[index_ii], [84]) - np.percentile(log_Ms_all[index_ii], [50]))

            scatter.append(scatter_i)
            #print(np.percentile(log_Ms_all[index_ii], [84]),np.percentile(log_Ms_all[index_ii], [50]))
            #print(scatter_i)

        median_Ms = np.array(median_Ms)
        scatter = np.array(scatter)

        self.a_target = a_target
        self.median_Ms = median_Ms
        self.scatter = scatter

        # let's calculate chi-squared between us and Behroozi

        log_median_Ms = [log10(x) for x in median_Ms]

        log_Behroozi_Ms = [log10(x) for x in Behroozi_fusion[:, 1]]

        log_median_Ms = np.array(log_median_Ms)
        log_Behroozi_Ms = np.array(log_Behroozi_Ms)

        chi = np.sum((log_median_Ms[a_index_us] - log_Behroozi_Ms[a_index_Behroozi]) ** 2)

        slope = np.sum(log_Behroozi_Ms[a_index_Behroozi] * log_median_Ms[a_index_us]) / np.sum(
            log_median_Ms[a_index_us] ** 2)

        self.slope = slope

        self.chi = chi

        # Let's plot median_us vs Behroozi:

        plt.plot([], [], "ko", label="$\sigma\quad log[M_{*}]_{z=0} = %.3f$" % (scatter[-1]*slope))

        plt.plot([], [], "ko", label="$STD = %.3f$" % ((chi / len(scatter)) ** 0.5))

        plt.plot(Behroozi_fusion[:, 0][a_index_Behroozi], log_Behroozi_Ms[a_index_Behroozi], "ro", markersize=10,
                 label="$Behroozi \quad 2013$")

        plt.plot([], [], "k", label="$Subsample \quad of \quad Haloes$", alpha=alpha)

        plt.plot([], [], "ko", label="$log[M_h] \quad = %s (dex)$" % (self.index.replace("M", "")))

        # Do remember to update the value of f0

        # us:



        plt.plot(a_target, slope * log_median_Ms, "b", linewidth=3, label="$Median \quad value$")

        # Plot sub_halo, which has slope term:
        # This is for simultaneously fitting!


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

            log_us = np.array([log10(x) for x in Ms])

            plt.plot(a, slope * log_us, "k", alpha=alpha, linewidth=1)

        plt.legend()

        # save it:


        # save it:
        plt.suptitle("$log[M_*]\quad vs\quad scale \quad factor\quad new \quad f_{con}$")

        plt.xlabel("$Scale\quad factor\quad a$")
        plt.ylabel("$log[M_*]\quad (dex)$", fontsize=20)

        axes = plt.gca()
        axes.set_xlim([0.2, 1])
        axes.set_ylim([7, 12.5])
        # axes.set_ylim([6, 12])

        # share x axis


        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(18.5, 12.5)

        save_path = plot_path + "Fitting_sim_v2" + str(self.index).replace(".","") + ".png"
        fig.savefig(save_path, dpi=200)

        plt.close()

    def Plot_20_Halo_f_con(self):

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 15}

        matplotlib.rc('font', **font)

        alpha = 0.05

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


            f_con_array = np.array(f_con_array)

            if halo<20:
                plt.plot(a, f_con_array, "k", alpha=alpha, linewidth=1)

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
        f_con_all = result_all[:, 3]

        median_Ms = []

        median_f_con = []

        scatter = []

        for ii in range(0, len(a_target)):
            index_ii = np.where(a_all == a_target[ii])

            median_Ms.append(np.nanmedian(Ms_all[index_ii]))

            median_f_con.append(np.nanmedian(f_con_all[index_ii]))

            log_Ms_all = np.array([log10(x) for x in Ms_all])

            scatter_i = abs(np.percentile(log_Ms_all[index_ii], [84]) - np.percentile(log_Ms_all[index_ii], [50]))

            scatter.append(scatter_i)

        median_Ms = np.array(median_Ms)
        median_f_con = np.array(median_f_con)
        scatter = np.array(scatter)

        self.a_target = a_target
        self.median_Ms = median_Ms
        self.scatter = scatter

        # let's calculate chi-squared between us and Behroozi

        log_median_Ms = [log10(x) for x in median_Ms]

        log_Behroozi_Ms = [log10(x) for x in Behroozi_fusion[:, 1]]

        log_median_Ms = np.array(log_median_Ms)
        log_Behroozi_Ms = np.array(log_Behroozi_Ms)

        chi = np.sum((log_median_Ms[a_index_us] - log_Behroozi_Ms[a_index_Behroozi]) ** 2)

        self.chi = chi

        # Median of our f_con
        plt.plot(a_target, median_f_con, "bo", label = "$Median\quad of\quad Our\quad f_{con}$")


        # Let's plot median_us vs Behroozi:

        plt.plot(Behroozi_fusion[:, 0][a_index_Behroozi], Behroozi_fusion[:, 5][a_index_Behroozi], "ro", markersize=10,
                 label="$Behroozi \quad 2013$")

        plt.plot([], [], "k", label="$Our\quad f_{con} \quad of \quad some\quad halos$", alpha=alpha)

        plt.plot([], [], "ko", label="$log[M_h] \quad = %s (dex)$" % (self.index.replace("M", "")))

        plt.legend()

        # save it:


        # save it:

        plt.suptitle("$f_{con}\quad vs\quad scale \quad factor$")

        plt.xlabel("$Scale\quad factor\quad a$")
        plt.ylabel("$f_{con}$", fontsize=20)

        axes = plt.gca()
        axes.set_xlim([0.2, 1])
        # axes.set_ylim([0, 0.1+np.nanmax(Behroozi_fusion[:, 5][a_index_Behroozi])])
        axes.set_ylim([0, 0.1 + np.max([np.nanmax(Behroozi_fusion[:, 5][a_index_Behroozi]),np.nanmax(median_f_con)])])

        # share x axis


        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(18.5, 12.5)

        save_path = plot_path + "f_con_Fitting_sim_v2" + str(self.index).replace(".","") + ".png"
        fig.savefig(save_path, dpi=200)

        plt.close()

    def fitting_separately(self, index):

        self.index = index


class diagnostic():
    def __init__(self, chain):

        name = ["f0", "A1", "A2", "A3", "A4","A5", "A6"]

        self.chain = chain
        self.name = name

    def walker_plot(self):

        name = self.name

        chain = self.chain
        n_parameter = len(chain[0, 0, :])

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 15}

        matplotlib.rc('font', **font)

        for i in range(0, n_parameter):

            plt.subplot(n_parameter, 1, i + 1)

            alpha = .3

            data = np.array(chain[:, :, i]).ravel()


            for j in range(0, len(chain[:, 0, 0])):
                plt.plot(chain[j, :, i], "k", linewidth=0.7, alpha=alpha)

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


        save_path = "/Users/caojunzhi/Downloads/upload_201706_Jeremy/" + "Walker_plot_sim_v1_" + ".png"
        fig.savefig(save_path, dpi=200)

        plt.close()

    def corner_plot(self, fusion):

        name = ["f0", "A1", "A2", "A3", "A4", "A5","A6"]

        # Only choose resultf form steps after 50


        samples = fusion.reshape(-1, 7)

        print(samples.shape)

        f0, A1, A2, A3, A4, A5,A6 = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                                    zip(*np.percentile(samples, [16, 50, 84],
                                                                       axis=0)))

        print(f0[0], A1[0], A2[0], A3[0], A4[0], A5[0], A6[0])

        self.f0 = f0
        self.A1 = A1
        self.A2 = A2
        self.A3 = A3
        self.A4 = A4
        self.A5 = A5
        self.A6 = A6

        print("f0")
        print(f0)

        print("A1")
        print(A1)
        print("A2")
        print(A2)
        print("A3")
        print(A3)
        print("A4")
        print(A4)
        print("A5")
        print(A5)
        print("A6")
        print(A6)

        fig = corner.corner(samples, labels=name,
                            truths=[f0[0], A1[0], A2[0], A3[0], A4[0], A5[0], A6[0]])

        # fig = corner.corner(samples, labels=name)


        fig.savefig("/Users/caojunzhi/Downloads/upload_201706_Jeremy/" + "triangle_sim_v1.png")


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


# Threshold
all_index = ["M11.0", "M11.5", "M12.0", "M12.5", "M13.0"]

model = MCMC(kwargs=kwargs)

index_i = all_index[0]
model.fitting_separately(index=index_i)


# Path

save_path = "Simultaneously_no_z_best_fit_v2.pkl"

save_path_backup = save_path


save_path = save_path.replace("best_fit_", "EMCEE")


pkl_file = open(save_path, 'rb')
result = pickle.load(pkl_file)
pkl_file.close()


pkl_file = open(save_path.replace("EMCEE", "EMCEE_chi_"), 'rb')
chi = pickle.load(pkl_file)
pkl_file.close()

chi = chi.T


# print(result.shape)

# Load result:

## The least chi-squared result


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





# Or use scipy best fitting result:

## load scipy best fitting result:

pkl_file = open(save_path_backup, 'rb')
best_fit_result = pickle.load(pkl_file)
pkl_file.close()


f0, A1, A2, A3, A4, A5, A6 = best_fit_result
print("scipy best fitting result")
print(f0, A1, A2, A3, A4, A5, A6)


kwargs["f0"] = f0
kwargs["A1"] = A1
kwargs["A2"] = A2
kwargs["A3"] = A3
kwargs["A4"] = A4


kwargs["A5"] = A5
kwargs["A6"] = A6



# Fix A_s for different halos and try it.

model.update_kwargs(kwargs=kwargs)


model.Plot_M_stellar()
model.Plot_20_Halo_f_con()

# diagnostic plot:


pkl_file = open(save_path, 'rb')
result = pickle.load(pkl_file)
pkl_file.close()



fusion =result
model = diagnostic(chain=result)

model.walker_plot()

model.corner_plot(fusion=fusion)














