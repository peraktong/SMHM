import numpy as np
import pickle
import matplotlib.pyplot as plt
from termcolor import colored
import emcee
import scipy.optimize as op

# example:
# For our circumstance, we can still improve


# Let's define some functions

count=0

def lnlike(theta, x, y, yerr):
    global count
    count = count+1
    print(count)
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

# The final ln posterior probability is the sum of the ln prior and ln likelihood

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)



# Let's create some fake data:
# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)


## This part is for scipy opt
## We need this because we need them as initial values!!
## That's what you discussed with Jeremy.

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))

print(result["x"])
m_ml, b_ml, lnf_ml = result["x"]


# Now let's come to the main topic:

# Define the initial condition of the MCMC chain: Position/initial values

ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# Set up the MCMC chain
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))

sampler.run_mcmc(pos, 500)

# Now we have an array with dimension of 100*500*3: 100 walkers, 500 steps and 3 parameters

result = sampler.chain


save_path = "/Users/caojunzhi/Downloads/upload_20170516_Jeremy/standard_examples"
# plots












