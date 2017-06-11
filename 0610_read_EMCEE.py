import numpy as np
import pickle

save_path = "Separately_EMCEEM11.0_v6.pkl"

pkl_file = open(save_path, 'rb')
result = pickle.load(pkl_file)
pkl_file.close()

save_path = save_path.replace("EMCEE", "EMCEE_chi_")

pkl_file = open(save_path, 'rb')
chi = pickle.load(pkl_file)
pkl_file.close()

chi = chi.T

print(result.shape,chi.shape)



result = result.reshape(-1,7)
f0, A1, A2, A3, A4, A5, zs = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                 zip(*np.percentile(result[50:,:], [16, 50, 84],
                                                    axis=0)))

print(f0[0], A1[0], A2[0], A3[0], A4[0], A5[0], zs[0])




"""

chi = np.array(chi,dtype=float)
index = np.argwhere(chi == np.nanmin(chi))
i,j = index[0]

f0, A1, A2, A3, A4, A5, zs = result[i,j,:]
print(f0,A1,A2,A3,A4,A5,zs)
print(chi[i,j])


"""