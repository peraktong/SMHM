import numpy as np
import pickle
import matplotlib.pyplot as plt
import math
import matplotlib
import corner

ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/"


pkl_file = open('B_M13.pkl', 'rb')
B_M13 = pickle.load(pkl_file)
pkl_file.close()



pkl_file = open('redshift_result_M13.pkl', 'rb')
kwargs = pickle.load(pkl_file)
pkl_file.close()

# load emcee:

pkl_file = open('EMCEE_redshift_M13.pkl', 'rb')
emcee = pickle.load(pkl_file)
pkl_file.close()

# Only choose resultf form steps after 50
samples = emcee[:, 50:, :].reshape((-1, 4))

samples[:, 2] = np.exp(samples[:, 2])
f0,zc,sigmaz,alphaz = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))

f0 = np.median(f0)
zc = np.median(zc)
sigmaz = np.median(sigmaz)
alphaz = np.median(alphaz)

"""

# [ 0.88009791  3.34350292  1.63356693  1.00601424]


kwargs["f0"] = 0.88009791
kwargs["zc"] = 3.34350292
kwargs["sigmaz"] = 1.63356693
kwargs["alphaz"] = 1.00601424



kwargs["f0"] = f0
kwargs["zc"] = zc
kwargs["sigmaz"] = sigmaz
kwargs["alphaz"] = alphaz



"""



kwargs["f0"] = f0
kwargs["zc"] = zc
kwargs["sigmaz"] = sigmaz
kwargs["alphaz"] = alphaz




def log10(x):
    if x > 0:
        return math.log10(x)
    else:
        return float("nan")


log10 = np.vectorize(log10)


def exp(x):
    try:
        return math.exp(x)
    except OverflowError:
        return float("nan")


exp = np.vectorize(exp)


class read_and_plot():
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



    ## Here comes f_q

    def f_q(self,z=None,Mh=None,Ms=None):

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

            #print(abs(((z-zc)/sigmaz)**alphaz))
            if z>zc:
                return 1
            else:
                return (exp((z - zc) / sigmaz)) ** alphaz



        elif ty == "halo":
            return exp(((log10(mhc)-log10(Mh))/sigmah)**alphah)

        elif ty=="galaxy":
            return exp(((log10(msc)-log10(Ms))/sigmag)**alphag)

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
        #print(z_array)
        #print(f_q_array_temp)
        #print(f_q_array_temp.shape,M_stellar_interval.shape)

        #multiply them (f_q)

        M_stellar_interval = M_stellar_interval*f_q_array_temp

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



    def input_array(self, z_array, delta_mh_array):

        self.z_array = z_array
        self.delta_mh_array = delta_mh_array

    # save each of them in small files

    # Index is M11, M12, M13 or M14
    # Now we have 300 halos for our tree

    # Let's calculate Behroozi 2013 results:

    def f(self, x, alpha, delta, gamma):
        return -log10(10 ** (alpha * x) + 1) + delta * ((log10(1 + exp(x))) ** gamma / (1 + exp(10 ** (-x))))

    def log_M_stellar(self, epsilon, M1, Mh, delta, alpha, gamma):

        return log10(epsilon * M1) + self.f(x=log10(Mh / M1), alpha=alpha, delta=delta, gamma=gamma) - self.f(x=0,
                                                                                                              alpha=alpha,
                                                                                                              delta=delta,
                                                                                                              gamma=gamma)

    def read_kwargs(self,kwargs):

        self.kwargs = kwargs
        self.kwargs_temp = kwargs

    def read_save_plot(self, index):

        mass_name = index

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 20}

        matplotlib.rc('font', **font)

        trans = 0.05
        n_halo = 100

        data_path = ori + index + "/"

        # result all is z + Mass

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

            # set f0
            f0 = 0.1
            self.input_array(z_array=z, delta_mh_array=delta_Mh)
            M_stellar = self.calculate_mass()

            fusion = np.c_[z, M_stellar]

            stellar_mass = fusion[:, 1]

            # plot

            # plt.plot(1/(1+fusion[:,0]),[log10(x) for x in stellar_mass],"k",alpha=alpha,linewidth=1)

            try:
                result_all = np.vstack((result_all, fusion))

            except:
                none = 1

        # add median

        # result_all : z + Mh
        result_all = np.array(result_all)

        # print("result shape")
        # print(result_all.shape)

        # since all of them are from Main branch, it's okay to do this

        target = z_target

        target = list(set(target))

        M_stellar_all = []

        ## ?? Use 68% percent range of values?
        residual = []


        for i in range(0, len(target)):
            # print("Doing %.2f %%" % (i / len(target) * 100))

            index = np.where(abs(result_all[:, 0] - target[i]) < 0.01)
            index = np.array(index)
            index = index.ravel()

            # Mass
            array = result_all[:, 1][index]



            # in dex
            # 68% range of values

            array_dex = np.array([log10(x) for x in array])
            scatter_i = np.std(array_dex)
            mask_i = (array_dex<np.nanmedian(array_dex)+scatter_i)&(array_dex>np.nanmedian(array_dex)-scatter_i)

            residual.extend(array_dex[mask_i]-np.nanmedian(array_dex))


            # print(target[i])
            # print(array)
            M_stellar_all.append(np.nanmedian(array))

        M_stellar_all = np.array(M_stellar_all)
        residual = np.array(residual).ravel()

        # print(scatter)

        # Let's do it!

        # Plot Behroozi 2013

        a_Behroozi = B_M13[:,0]
        m_Behroozi = B_M13[:,1]

        plt.plot(a_Behroozi, m_Behroozi, "ro", linewidth=3, label="$Behroozi \quad 2013$")

        plt.plot([], [], "k", label="$Subsample \quad of \quad Haloes$", alpha=trans)
        plt.plot([], [], "ko", label="$log[M_h] \quad = %.f (dex)$" % (np.nanmax([log10(x) for x in M_h])))

        M_stellar_all = np.array(M_stellar_all)

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

        # Also save this for convenient:


        # Let's fit us to B_2013, since our result is proportional to logf0, we can use a simple trick
        # Only fit the B_2013[mask] point


        # another mask for Behroozi for M13 and M14

        """

        mask2 = a_Behroozi<0.35

        a_Behroozi = a_Behroozi[mask2]
        m_Behroozi = m_Behroozi[mask2]


        """

        target = a_Behroozi
        index_us = []

        for mk in range(0, len(target)):
            index = np.where(a_us == target[mk])

            index_us.append(index)

        index_us = np.array(index_us, dtype=int)

        a_us_chi = a_us[index_us]
        m_us_chi = m_us[index_us]

        # y_b = slope*x_us, find a via minimizing chi-squared

        # Hogg's linear regression document:


        m_Behroozi = np.array(m_Behroozi).ravel()
        m_us_chi = np.array(m_us_chi).ravel()

        slope = 1



        plt.plot(a_us, slope * m_us, "b", linewidth=3, label="$Median \quad value$")



        std = (np.sum((m_Behroozi - slope * m_us_chi) ** 2) / len(m_Behroozi))**0.5

        # print(residual.shape)
        scatter = (np.sum(residual**2)/len(residual))**0.5

        plt.plot([], [], "ko", label="$\sigma\quad log[M_{*}] = %.3f$" % (scatter))

        plt.plot([], [], "ko", label="$STD = %.3f$" % (std))

        # plot sub-halos


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

            self.input_array(z_array=z, delta_mh_array=delta_Mh)
            M_stellar = self.calculate_mass()

            # plot
            # Add slope



            plt.plot(1 / (1 + z), [slope * log10(x) for x in M_stellar], "k", alpha=trans, linewidth=1)

            try:
                result_all = np.vstack((result_all, fusion))

            except:
                none = 1

        plt.legend()

        # save it:
        plt.suptitle("$log[M_*]\quad vs\quad scale \quad factor$")

        plt.xlabel("$Scale\quad factor\quad a$")
        plt.ylabel("$log[M_*]\quad (dex)$", fontsize=20)


        """
           
        
        axes = plt.gca()
        axes.set_xlim([0.2, 1])
        axes.set_ylim([9, 12])
        
   
        
        """

        # share x axis



        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(18.5, 12.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170516_Jeremy/" + "redshift_emcee_" + str(mass_name) + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


# Let's define some diagnostic plots:


class diagnostic():

    def __init__(self,chain,name,save_path):

        self.chain = chain
        self.name = name

        self.save_path = save_path

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

            plt.xlabel("Step")
            plt.ylabel("%s" % name[i], fontsize=20)


            
        axes = plt.gca()
        axes.set_xlim([15660, 15780])

            


        plt.legend()

        plt.suptitle("The values of walkers vs step", fontsize=20)

        # share x
        plt.subplots_adjust(hspace=.0)

        # plt.show()
        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(16.5, 11.5)

        # select some plots:


        save_path = "/Users/caojunzhi/Downloads/upload_20170516_Jeremy/" + "redshift_walker" + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()

    def corner_plot(self):
        name = self.name

        fig = corner.corner(samples, labels=name,
                            truths=[f0,zc,sigmaz,alphaz])


        fig.savefig(self.save_path+"triangle.png")

## Output:

model = read_and_plot()

# model.read_save_plot(index="M11")

# Do use M13

index = "M13"
model.read_kwargs(kwargs=kwargs)

print(kwargs)

model.read_save_plot(index=index)



name = ["f0","zc","sigma_z","alpha_z"]

save_path = "/Users/caojunzhi/Downloads/upload_20170516_Jeremy/"

model = diagnostic(chain=emcee[:, 50:, :],name=name,save_path=save_path)

model.walker_plot()
model.corner_plot()









