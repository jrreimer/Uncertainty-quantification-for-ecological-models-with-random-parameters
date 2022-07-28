# Code accompanying "Uncertainty quantification for ecological models with random parameters"
# Reproduces partial data for Figure 4: saves this data into "scenario1.mat"

import sys

import numpy as np
from scipy import integrate
from scipy.io import savemat

from UncertainSCI.distributions import BetaDistribution
from UncertainSCI.indexing import TotalDegreeSet
from UncertainSCI.pce import PolynomialChaosExpansion
from UncertainSCI.utils.linalg import weighted_lsq

from np_algae_model import np_solver

a = 0
b = 1
c = 0.8
d = 0.1
e = 0

x0 = np.zeros(2) # [N(0) P(0)]
x0[0] = 0.22
x0[1] = 0.001

# Time integration
tmax = 160.; # max. time value
t_eval = np.arange(0, tmax, 0.05)
N = 2*t_eval.size

solver = np_solver((0, tmax), t_eval)

order = 6

quantile_levels = [0.05, 0.5, 0.95]

# Number of quad points
M = order+1

### UncertainSCI setup
dim = 1
bounds = np.zeros([2,dim])
bounds[0], bounds[1] = 0.8, 1.2 # bounds for b

alpha, beta = 1, 1
bdist = BetaDistribution(domain=np.reshape(bounds, [2, 1]), alpha=alpha, beta=beta)

iset = TotalDegreeSet(dim=dim,order=order)

pce = PolynomialChaosExpansion(iset, bdist)

# Generate samples
p, w = bdist.polys.polys1d[0].gauss_quadrature(M)
samples = bdist.transform_to_standard.mapinv(bdist.transform_standard_dist_to_poly.mapinv(p))
pce.set_samples(np.reshape(samples, [M, dim]))

# Collect model output
model_output = np.zeros([pce.samples.shape[0], N])
args = np.zeros(5)
args[0] = a
args[2:] = [c, d, e]
for ind, xval in enumerate(pce.samples):
    args[1] = xval
    model_output[ind,:] = solver(x0, args)['y'].flatten()

# More complicated code to circumvent LOOCV computation
pce.model_output = model_output
p_standard = pce.map_to_standard_space(pce.samples)
V = pce.distribution.polys.eval(p_standard, pce.index_set.get_indices())
pce.coefficients, pce.accuracy_metrics['residuals'] = weighted_lsq(V, pce.model_output, pce.weights)

#pce.build(model_output=model_output)

pce_mean = np.reshape(pce.mean(), [2, t_eval.size])
pce_stdev = np.reshape(pce.stdev(), [2, t_eval.size])
temp = pce.quantile(quantile_levels, M=int(1e4))
pce_quantiles = [None,]*len(quantile_levels)
for qind in range(len(quantile_levels)):
    pce_quantiles[qind] = np.reshape(temp[qind,:], [2, t_eval.size])

### Monte Carlo runs
NMC = [11, 10000]

mc_mean = [None,]*len(NMC)
mc_stdev = [None,]*len(NMC)
mc_samples = [None,]*len(NMC)
mc_quantiles = []

for mcind, mcsize in enumerate(NMC):

    mc_samples[mcind] = bdist.MC_samples(mcsize)

    model_output = np.zeros([mcsize, N])
    args = np.zeros(5)
    args[0] = a
    args[2:] = [c, d, e]
    for ind, xval in enumerate(mc_samples[mcind]):
        args[1] = xval

        model_output[ind,:] = solver(x0, args)['y'].flatten()

    mc_mean[mcind] = np.reshape(np.mean(model_output, axis=0), [2, t_eval.size])
    mc_stdev[mcind] = np.reshape(np.std(model_output, axis=0), [2, t_eval.size])
    temp = []
    for qind in range(len(quantile_levels)):
        qlevel = quantile_levels[qind]
        temp.append(np.reshape(np.quantile(model_output, qlevel, axis=0), [2, t_eval.size]))
    mc_quantiles.append(temp)

### Saving data in Matlab format
mdic = {'pce_mean':pce_mean, 
        'pce_stdev':pce_stdev, 
        'mc_mean1':mc_mean[0], 
        'mc_mean2':mc_mean[1], 
        'mc_stdev1':mc_stdev[0], 
        'mc_stdev2':mc_stdev[1], 
        'mc_samples1':mc_samples[0], 
        'mc_samples2':mc_samples[1], 
        'mc_quantiles00': mc_quantiles[0][0], 
        'mc_quantiles01':mc_quantiles[0][1], 
        'mc_quantiles02':mc_quantiles[0][2], 
        'mc_quantiles10':mc_quantiles[1][0], 
        'mc_quantiles11':mc_quantiles[1][1], 
        'mc_quantiles12':mc_quantiles[1][2]}

savemat("scenario1.mat", mdic)
print("Saved data to {0:s}".format("scenario1.mat"))
