import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle

x = [11,11.5,12,12.5]

index_array = ["M11.0","M11.5","M12.0","M12.5"]

name = ["f0","A1","A2","A3","A4","A5","A6"]



kwargs = {"MCMC": None, "A1": None, "A2": None, "A3": None, "f0": None, "A4": None, "A5": None, "A6": None, "mht": None,
          "mst": None, "zc": None, "sigmaz": None,
          "alphaz": None, "mhc": None, "msc": None, "sigmah": None, "sigmag": None, "alphah": None, "alphag": None}



f0 = []

A1 = []

A2 = []

A3 = []

A4 = []

A5 = []

A6 = []

for i in range(0,4):
    # Path
    save_path = "Fast_Separately_best_fit_M11.0_v2.pkl"
    save_path = save_path.replace("M11.0", str(index_array[i]))

    save_path = save_path.replace("best_fit_", "EMCEE")

    pkl_file = open(save_path, 'rb')
    result = pickle.load(pkl_file)
    pkl_file.close()

    pkl_file = open(save_path.replace("EMCEE", "EMCEE_chi_"), 'rb')
    chi = pickle.load(pkl_file)
    pkl_file.close()

    chi = chi.T

    ## The least chi-squared result




    chi = np.array(chi, dtype=float)
    index = np.argwhere(chi == np.nanmin(chi))
    i, j = index[0]

    f0_i, A1_i, A2_i, A3_i, A4_i, A5_i, A6_i = result[i, j, :]

    print(f0_i, A1_i, A2_i, A3_i, A4_i, A5_i, A6_i)
    print(chi[i, j])

    f0.append(f0_i)
    A1.append(A1_i)
    A2.append(A2_i)
    A3.append(A3_i)
    A4.append(A4_i)
    A5.append(A5_i)
    A6.append(A6_i)



font = {'family': 'normal',
        'weight': 'bold',
        'size': 20}

matplotlib.rc('font', **font)


f, axarr = plt.subplots(2, 4,sharex=True)

axarr[0, 0].plot(x,f0,"o",label=name[0])
axarr[0, 0].set_title(name[0])

axarr[0,0].set_xlabel("log[M_h](z=0)")
axarr[0,0].set_ylabel(name[0])

axarr[0, 1].plot(x,A1,"o",label=name[1])
axarr[0, 1].set_title(name[1])

axarr[0,1].set_xlabel("log[M_h](z=0)")
axarr[0,1].set_ylabel(name[1])


axarr[0, 2].plot(x,A2,"o",label=name[2])
axarr[0, 2].set_title(name[2])
axarr[0,2].set_xlabel("log[M_h](z=0)")
axarr[0,2].set_ylabel(name[2])


axarr[0, 3].plot(x,A3,"o",label=name[3])
axarr[0, 3].set_title(name[3])
axarr[0,3].set_xlabel("log[M_h](z=0)")
axarr[0,3].set_ylabel(name[3])





axarr[1, 0].plot(x,A4,"o",label=name[4])
axarr[1, 0].set_title(name[4])
axarr[1,0].set_xlabel("log[M_h](z=0)")
axarr[1,0].set_ylabel(name[4])


axarr[1, 1].plot(x,A5,"o",label=name[5])
axarr[1, 1].set_title(name[5])
axarr[1,1].set_xlabel("log[M_h](z=0)")
axarr[1,1].set_ylabel(name[5])



axarr[1, 2].plot(x,A6,"o",label=name[6])
axarr[1, 2].set_title(name[6])
axarr[1,2].set_xlabel("log[M_h](z=0)")
axarr[1,2].set_ylabel(name[6])




plt.suptitle("Values of parameters vs log[Mh]")


fig = matplotlib.pyplot.gcf()

# adjust the size based on the number of visit

fig.set_size_inches(18.5, 12.5)

save_path = "/Users/caojunzhi/Downloads/upload_201706_Jeremy/" + "parameters_vs_halo_mass_faster_v2.png"

fig.savefig(save_path, dpi=500)

plt.close()