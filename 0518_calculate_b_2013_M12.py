import numpy as np
import pickle
import matplotlib.pyplot as plt
import math
import matplotlib

ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/"


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

    # Attention! This is for one Halo. This means the input z_array should be the z_array for one Halo!
    def calculate_mass(self, f0):
        # input array: f_con at zi, deltaMh at zi
        # return an array of M*

        z_array = self.z_array
        delta_mh_array = self.delta_mh_array

        M_stellar_interval = []

        for i in range(0, len(z_array)):
            M_stellar_interval.append(
                self.stellar_mass_per_interval(f_con=self.f_con_z(z_array[i], f0=f0), delta_Mh=delta_mh_array[i]))

        # add them togethrm

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

    def calculate_SM_a_Behroozi(self, index):

        data_path = ori + index + "/"

        n_halo = 300

        # The structure is a and log M_stellar
        result_all_Behroozi = np.array([[0, 0]])

        # What we need is 1 and 5
        # 1 Halo Mass
        # 5 red shift z

        for halo in range(1, n_halo):

            file_i = "output_fulltree" + str(halo) + ".txt"

            result_i = np.loadtxt(data_path + file_i)

            # calculate M_stellar

            M_h = result_i[:, 1]
            # print(M_h)


            z = result_i[:, 5]

            if halo == 1:
                z_target = z

            else:
                none = 1

            # We only need to plot the final matrix:

            #  a = 1/(1+z)

            #### parameters from Behroozi 2013

            a = 1 / (1 + z)

            v = exp(-4 * a ** 2)

            epsilon = 10 ** (-1.777 + (-0.006 * (a - 1) * v - 0.119 * (a - 1)))

            M1 = 10 ** (11.514 + (-1.793 * (a - 1) - 0.251 * z) * v)

            alpha = -1.412 + (0.731 * (a - 1)) * v

            delta = 3.508 + (2.608 * (a - 1) - 0.043 * z) * v

            gamma = 0.316 + (1.319 * (a - 1) + 0.279 * z) * v

            log_M_stellar = self.log_M_stellar(epsilon=epsilon, M1=M1, Mh=M_h, delta=delta, alpha=alpha, gamma=gamma)

            fusion = np.c_[z, log_M_stellar]

            # z, Mh and log M_stellar

            try:
                result_all_Behroozi = np.vstack((result_all_Behroozi, fusion))

            except:
                none = 1

        result_all_Behroozi = np.array(result_all_Behroozi)

        # print(result_all_Behroozi.shape)

        # Combine results:


        # since all of them are from Main branch, it's okay to do this
        target = z_target

        target = list(set(target))

        M_stellar_all = []

        stellar_mass_b = []

        for i in range(0, len(target)):
            # print("Doing %.2f %%" % (i / len(target) * 100))

            index = np.where(abs(result_all_Behroozi[:, 0] - target[i]) < 0.01)
            index = np.array(index)
            index = index.ravel()

            # stellar Mass in dex
            array = result_all_Behroozi[:, 1][index]

            # in dex

            stellar_mass_b.append(np.median(array))

        median_Behroozi = np.c_[target, stellar_mass_b]
        median_Behroozi = np.array(median_Behroozi)

        self.median_Behroozi = median_Behroozi

    def read_save_plot(self, index):

        mass_name = index

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 20}

        matplotlib.rc('font', **font)

        trans = 0.05
        n_halo = 300

        data_path = ori + index + "/"

        plot_path = "/Users/caojunzhi/Downloads/upload_2017.5.1_Jeremy/"

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
            M_stellar = self.calculate_mass(f0=f0)

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

        # Let's do it!

        # Plot Behroozi 2013

        median_Behroozi = self.median_Behroozi

        a_Behroozi = 1 / (1 + median_Behroozi[:, 0])
        m_Behroozi = median_Behroozi[:, 1]

        index = np.argsort(a_Behroozi)

        a_Behroozi = a_Behroozi[index]
        m_Behroozi = m_Behroozi[index]

        # important! Fit these points!

        # Put mask here! For M13...

        mask = (m_Behroozi > 0) & (a_Behroozi > 0.2)

        a_Behroozi = np.array(a_Behroozi[mask]).ravel()
        m_Behroozi = np.array(m_Behroozi[mask]).ravel()

        #### save Behroozi:
        # in fact, they are a and logm!
        B_M13 = np.c_[a_Behroozi,m_Behroozi]


        output = open('B_M12.pkl', 'wb')
        pickle.dump(B_M13, output)
        output.close()


model = read_and_plot()

# model.read_save_plot(index="M11")


index = "M12"
model.calculate_SM_a_Behroozi(index=index)
model.read_save_plot(index=index)






