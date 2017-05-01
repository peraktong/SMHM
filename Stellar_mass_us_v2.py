import numpy as np
import pickle
import matplotlib.pyplot as plt
import math
import matplotlib
ori = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/2017.4_9_cosmology/data/data_all/"
save = "/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/result_z_vs_m/"

def log10(x):
    return math.log10(x)


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
        #print(M_stellar_interval)
        #print(M_stellar)
        return M_stellar

    def input_array(self, z_array, delta_mh_array):

        self.z_array = z_array
        self.delta_mh_array = delta_mh_array

    def read_behroozi(self):

        path = ori + "behroozi_lookup.dat"

        # read file:

        result = np.array([0, 0, 0])
        for line in open(path, 'r'):
            data = line.rstrip().split(" ")

            result_i = []
            for i in data:
                if i == "":
                    none = 1
                else:
                    result_i.append(i)
            result_i = np.array(result_i, dtype=float)
            if len(result_i) == 3:
                result = np.vstack((result, result_i))

        print(result.shape)

        self.result_b = result



    # save each of them in small files
    def read_ours(self,index):

        pkl_file = open('/Users/caojunzhi/Desktop/NYU/Laboratory/My_code/Cosmology_2017.4-8_Jason/tree_1halo_1.pkl', 'rb')
        redshift = pickle.load(pkl_file)
        pkl_file.close()

        self.redshift = redshift[:,0]

        self.index = index

        path = ori + "parsed_tree.dat".replace(".dat",str(index)+".dat")

        # read file:

        result_tree = np.array([0, 0, 0])
        halo =1


        for line in open(path, 'r'):
            data = line.rstrip().split(" ")
            # print(data)

            result_i = []
            for i in data:
                if i == "":
                    none = 1
                else:
                    result_i.append(i)

            result_i = np.array(result_i, dtype=float)
            #print(result_i)



            if len(result_i) == 3:
                result_tree = np.vstack((result_tree, result_i))

            elif (len(result_i)>5)&(len(np.atleast_2d(result_tree)[:,0])>1):

                # save it

                length = len(result_tree[:, 0])
                result_tree = result_tree[1:, :]
                #print(halo)
                #print(result_tree.shape)

                # calculate fusion matrix: Z vs M_stellar


                z = result_tree[:, 0]

                Mh = result_tree[:, 1]

                delta_Mh = []

                for j in range(0, len(Mh) - 1):
                    delta_Mh.append(Mh[j] - Mh[j + 1])

                delta_Mh.append(0)
                delta_Mh = np.array(delta_Mh)

                # fusion z Mh delta Mh  + (M_stellar)

                # Then calculate median

                # calculate mass:

                # if u want to fit the function to the data in B_2013

                self.input_array(z_array=z, delta_mh_array=delta_Mh)


                # set f0 from paper https://arxiv.org/pdf/1207.6105.pdf

                f0 = 0.01
                M_stellar = self.calculate_mass(f0=f0)

                fusion = np.c_[z, M_stellar]

                # The first index_halo =1 + tree_index=1 file has the correct z scale

                self.fusion = fusion

                # tree_index_halo_index

                save_path = save + "tree_" + str(index) + "halo_" + str(halo) + ".pkl"

                self.save_path = save + "tree_" + str(index) + "halo_"


                output = open(save_path, 'wb')
                pickle.dump(fusion, output)
                output.close()

                halo += 1
                self.halo_index_max = halo

                result_tree = np.array([0, 0, 0])

    def read_save_plot(self):

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 20}

        matplotlib.rc('font', **font)

        alpha = 0.1

        plot_path = "/Users/caojunzhi/Downloads/upload_2017.5.1_Jeremy/"

        result_all = np.array([[0, 0]])
        for halo in range(1,self.halo_index_max):


            pkl_file = open(self.save_path+str(halo)+".pkl", 'rb')
            result_i  = pickle.load(pkl_file)
            pkl_file.close()


            try:
                result_all = np.vstack((result_all, result_i))

            except :
                none=1


            M_stellar = result_i[:,1]

            mask = M_stellar <= 0
            M_stellar[mask] = None

            # plot

            plt.plot(1/(1+result_i[:,0]),[log10(x) for x in M_stellar],"k",alpha=alpha,linewidth=1)


        # add median

        # z + Mh
        result_all = np.array(result_all)

        print("result shape")
        print(result_all.shape)

        target = self.redshift

        target = list(set(target))

        M_stellar_all = []
        scatter = []

        for i in range(0, len(target)):
            # print("Doing %.2f %%" % (i / len(target) * 100))

            index = np.where(abs(result_all[:,0] -target[i])<0.01)
            index = np.array(index)
            index = index.ravel()

            # Mass
            array = result_all[:,1][index]

            mask = array>0

            scatter.append(np.std([log10(x) for x in array[mask]]))

            #print(target[i])
            #print(array)
            M_stellar_all.append(np.nanmedian(array))

        # Let's do it!

        scatter = np.array(scatter)

        scatter = (np.nansum(scatter**2)/len(scatter))**2
        print(scatter)

        M_stellar_all = np.array(M_stellar_all)

        mask = M_stellar_all<=0

        M_stellar_all[mask] = None

        target = np.array(target)
        #print(target)
        #print(M_stellar_all)


        plt.plot(1/(1+target),[log10(x) for x in M_stellar_all],"bo",linewidth =3,label="$Median \quad value$")
        plt.plot([],[],"ko",label = "$Subsample \quad of \quad Halo$",alpha=alpha)
        plt.plot([], [], "ko", label="$log[M_h] \quad = %.3f$"%(np.nanmax([log10(x) for x in M_stellar_all])))




        # calculate scatter
        scatter = scatter
        plt.plot([], [], "ko", label="$\sigma log[M_*]\quad = %.3f$" % (scatter))

        plt.legend()



        target = 1

        # save it:
        plt.suptitle("$log[M_*]\quad vs\quad scale \quad factor$")



        plt.xlabel("$Scale\quad factor\quad a$")
        plt.ylabel("$Flux$", fontsize=20)
        #axes = plt.gca()
        #axes.set_xlim([0,1])
        # share x axis



        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(18.5, 12.5)

        save_path = plot_path + "tree_"+str(self.index) + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


model = read_and_plot()

choose = [12,123,225,367,590,774,982,889]

for i in range(990,1000):
    model.read_ours(index=i)
    model.read_save_plot()


for j in choose:
    model.read_ours(index=j)
    model.read_save_plot()



