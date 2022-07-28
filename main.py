# Code accompanying "Uncertainty quantification for ecological models with random parameters"
# Reproduces Figures 5, S3, S4, and S5

from itertools import chain, combinations, product

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

from UncertainSCI.distributions import BetaDistribution, NormalDistribution, TensorialDistribution
from UncertainSCI.pce import PolynomialChaosExpansion

from np_algae_model import np_solver, plot_labels
from utils import normal_parameters_from_lognormal, plot_spread, compute_batch_skewness

# Change visualization to serif fonts
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times']
rcParams['pdf.fonttype'] = 42

def get_tensorized_grid(pce, M1d):
    """
    Returns tensorized Gauss quadrature using M1d points per dimension.
    """

    dist = pce.distribution

    ps, ws = [], []
    for dind in range(dim):
        tempp, tempw = dist.polys.polys1d[dind].gauss_quadrature(M1d)
        ps.append(tempp)
        ws.append(tempw)

    # Tensorize
    pall = np.meshgrid(*ps)
    wall = np.meshgrid(*ws)
    w = np.ones(M1d**dim)
    for dind in range(dim):
        pall[dind] = pall[dind].flatten()
        w *= wall[dind].flatten()
    p = np.vstack(pall).T

    p = dist.transform_to_standard.mapinv(dist.transform_standard_dist_to_poly.mapinv(p))
    p = np.reshape(p, [M1d**dim, dim])

    return p, w

# Problem setup
a = 0.00075
b = 1 # Nominal value of B
c = 0.8
d = 0.1
e = 0.0005

x0 = np.zeros(2) # [N(0) P(0)]
x0[0] = 0.05
x0[1] = 0.001

# Time integration
t0 = 0.  # time of interest
tmax = 150.;  # max. time value

t_eval = np.arange(t0, tmax, 0.5)
Nt = t_eval.size
solver = np_solver((t0, tmax), t_eval)

# UncertainSCI setup
dim = 3
order = 10

# Number of quad points
M = order+1

# Number of MC realizations to draw for plots
simNMC = 20

alpha, beta = 1, 1  # For a uniform distribution
bounds = np.zeros([2,1])
bounds[0,0], bounds[1,0] = 0.6, 1.4 # bounds for b
bdist = BetaDistribution(domain=np.reshape(bounds[:,0], [2, 1]), alpha=alpha, beta=beta)

NP_mu = np.array([0.22, 0.001])

N_cvs = [0.5, 0.3, 0.1]
P_cvs = [0.5, 0.3, 0.1]

comb = product(N_cvs, P_cvs)

# Manually arrange the 3x3 plots
rbg_green = [118/255., 176/255., 65/255.]
rbg_purple = [125/255., 91/255., 166/255.]

# Nominal + Mean/stdev plot
plt.figure(1, figsize=(12,9))
# CV + skewness plot
plt.figure(2, figsize=(12,9))
# Sensitivity plot
plt.figure(3, figsize=(12,9))

cv_axes = []
skew_axes = []

for (N0_cv, P0_cv) in comb:

    print("Running N0_cv={0:1.1e}, P0_cv={1:1.1e}".format(N0_cv, P0_cv))

    # Construct lognormal distributions + full 3-parameter distribution
    normal_mu, normal_cov = normal_parameters_from_lognormal(NP_mu, np.array([N0_cv, P0_cv]), 0.)
    GN = NormalDistribution(mean=[normal_mu[0]], cov=normal_cov[0,0]*np.eye(1))
    GP = NormalDistribution(mean=[normal_mu[1]], cov=normal_cov[1,1]*np.eye(1))
    dist = TensorialDistribution(distributions=(bdist, GN, GP))

    # Instantiate PC
    pce = PolynomialChaosExpansion(distribution=dist, order=order)

    # Assign samples
    pce.samples, w = get_tensorized_grid(pce, M)

    # Loop through samples for PC data
    model_output = np.zeros([pce.samples.shape[0], Nt])
    args = np.zeros(5)
    args[0] = a
    args[2:] = [c, d, e]
    for ind, xval in enumerate(pce.samples):
        args[1] = xval[0]
        x0[0] = np.exp(xval[1])
        x0[1] = np.exp(xval[2])

        model_output[ind,:] = solver(x0, args)['y'][1,:].flatten()

    # Generate Monte Carlo realizations for plotting
    args = np.zeros(5)
    args[0] = a
    args[2:] = [c, d, e]
    mc_outputs = np.zeros([simNMC, Nt])
    for mcind in range(simNMC):
        temp = dist.MC_samples(1)[0,:]
        args[1] = temp[0]
        x0[0] = np.exp(temp[1])
        x0[1] = np.exp(temp[2])
        mc_outputs[mcind,:] = solver(x0, args)['y'][1,:].flatten()

    # Gather nominal curve data
    args[1] = np.mean(bounds[:,0])
    x0[0] = NP_mu[0]
    x0[1] = NP_mu[1]
    nominal = solver(x0, args)['y'][1,:].flatten()

    # Build PC coefficients
    p_standard = dist.transform_standard_dist_to_poly.map(
                        dist.transform_to_standard.map(pce.samples))

    V = pce.distribution.polys.eval(p_standard,
                                    pce.index_set.get_indices())

    V = np.multiply(V.T, np.sqrt(w)).T
    model_output = np.multiply(model_output.T, np.sqrt(w)).T
    pce.coefficients, residuals = np.linalg.lstsq(V, model_output, rcond=None)[:2]

    # Construct global sensitivities
    variable_interactions = list(chain.from_iterable(combinations(range(dim), r) for r in range(1, dim+1)))
    global_sensitivity = pce.global_sensitivity(variable_interactions).T
    stdev = pce.stdev()
    mean = pce.mean()

    # Skewness: compute how many Gauss quadrature points we need to compute cubic terms
    M1d = int(np.ceil((3*order+1)/2))
    quadp, quadw = get_tensorized_grid(pce, M1d)
    skewness = quadw @ ((pce.pce_eval(quadp) - mean)/stdev)**3

    variable_texts = ["$B$", "$N_0$", "$P_0$"]
    labels = plot_labels(variable_interactions, variable_texts)

    ### Visualization

    # Manually arrange the 3x3 plots
    N0_index = N_cvs.index(N0_cv)
    P0_index = P_cvs.index(P0_cv)
    plot_index = 3*N0_index + P0_index + 1 

    ## Nominal + mean/stdev plots
    plt.figure(1)
    plt.subplot(3, 3, plot_index)
    plt.plot(t_eval, nominal, color=rbg_green, label='Nominal')
    plt.plot(t_eval, mean, color=rbg_purple, linestyle='--', label='Mean')
    plt.fill_between(t_eval, mean-stdev, mean+stdev, interpolate=True, facecolor=rbg_purple, alpha=0.2, label='$\pm 1$ stdev range')
    ymax = max(plot_spread(nominal)[1], plot_spread(mean)[1])
    for mcind in range(simNMC):
        if mcind==0:
            plt.plot(t_eval, mc_outputs[mcind,:], 'k', linewidth=1, alpha=0.1, label='Realizations')
        else:
            plt.plot(t_eval, mc_outputs[mcind,:], 'k', linewidth=1, alpha=0.1)
        ymax = max(ymax, plot_spread(mc_outputs[mcind,:])[1])
    plt.xlim([t0, tmax])
    plt.ylim([0, ymax])
    if (N0_index, P0_index) == (0,0):
        plt.legend(frameon=False,bbox_to_anchor=(0.5, 1))
    if P0_index == 0:
        plt.ylabel('algae $P$')
    if plot_index < 7:
        plt.tick_params(labelbottom=False)
    else:
        plt.xlabel('Time $t$ [days]')
    if N0_index == 0: # Then label columns with P0 CV
        plt.title(('CV($P_0$)=%f' % P0_cv).rstrip('0').rstrip('.'))
    if P0_index == 2: # Then label rows on right with N0 CV
        plt.text(1.1, 0.5, ('CV($N_0$)=%f' % N0_cv).rstrip('0').rstrip('.'), transform=plt.gca().transAxes )
    if P0_index > 0:
        plt.gca().tick_params(labelleft=False)

    ## CV + skewness plots
    plt.figure(2)
    plt.subplot(3, 3, plot_index)
    cv_line = plt.plot(t_eval, stdev/mean, color='k', label='CV')
    plt.xlim((t0, tmax))
    plt.ylim(plot_spread(stdev/mean))
    plt.gca().tick_params(axis='y', labelcolor='k')
    if plot_index < 7:
        plt.tick_params(labelbottom=False)
    else:
        plt.xlabel('Time $t$ [days]')
    cv_axes.append(plt.gca())

    # N_0 changes along rows
    # P_0 changes along columns
    if N0_index == 0: # Then label columns with P0 CV
        plt.title(('CV($P_0$)=%f' % P0_cv).rstrip('0').rstrip('.'))
    if P0_index == 2: # Then label rows on right with N0 CV
        plt.text(1.2, 0.5, ('CV($N_0$)=%f' % N0_cv).rstrip('0').rstrip('.'), transform=plt.gca().transAxes )
    if P0_index > 0:
        plt.tick_params(labelleft=False)
    else:
        plt.ylabel('CV', color='k')

    ax2 = plt.gca().twinx()
    skew_line = ax2.plot(t_eval, skewness, 'k', linestyle='--', label='skewness')
    ax2.set_ylim(plot_spread(skewness))
    ax2.tick_params(axis='y', labelcolor='k')
    skew_axes.append(ax2)
    if P0_index < 2:
        plt.tick_params(labelright=False)
    else:
        ax2.set_ylabel('Skewness', color='k')

    lines = cv_line + skew_line
    labs = [l.get_label() for l in lines]
    if (N0_index, P0_index) == (0,0):
        ax2.legend(lines, labs, frameon=False)

    ## Sensitivity plots
    plt.figure(3)
    plt.subplot(3, 3, plot_index)
    plt.stackplot(t_eval, global_sensitivity.T, labels=labels)
    plt.xlim((t0, tmax))
    plt.ylim((0, 1))
    if (N0_index, P0_index) == (0,2):
        plt.legend(frameon=False, bbox_to_anchor=(0.1, 1))
    if plot_index < 7:
        plt.tick_params(labelbottom=False)
    else:
        plt.xlabel('Time $t$ [days]')
    if P0_index > 0:
        plt.tick_params(labelleft=False)
    else:
        plt.ylabel('Global sensitivity $S_I$')
    if N0_index == 0: # Then label columns with P0 CV
        plt.title(('CV($P_0$)=%f' % P0_cv).rstrip('0').rstrip('.'))
    if P0_index == 2: # Then label rows on right with N0 CV
        plt.text(1.1, 0.5, ('CV($N_0$)=%f' % N0_cv).rstrip('0').rstrip('.'), transform=plt.gca().transAxes )


    ## Figure 5 plot
    if (P0_cv, N0_cv) == (0.3, 0.3):
        plt.figure(4, figsize=(10,6.5))
        plt.subplot(311)

        # Use trapezoid rule to compute averages
        dt = t_eval[1] - t_eval[0]
        nominal_avg = (np.sum(dt*nominal) - dt/2*nominal[0] - dt/2*nominal[-1])/(tmax-t0)
        mean_avg = (np.sum(dt*mean) - dt/2*mean[0] - dt/2*mean[-1])/(tmax-t0)

        label_str = r"Nominal, $\overline{P}" + "={0:1.2e}$".format(nominal_avg)
        plt.plot(t_eval, nominal, color=rbg_green, label=label_str)
        label_str = r"Mean, $\overline{P}" + "={0:1.2e}$".format(mean_avg)
        plt.plot(t_eval, mean, color=rbg_purple, linestyle='--', label=label_str)
        plt.fill_between(t_eval, mean-stdev, mean+stdev, interpolate=True, facecolor=rbg_purple, alpha=0.2, label='$\pm 1$ stdev range')
        ymax = max(plot_spread(nominal)[1], plot_spread(mean)[1])
        for mcind in range(simNMC):
            if mcind==0:
                plt.plot(t_eval, mc_outputs[mcind,:], 'k', linewidth=1, alpha=0.1, label='Realizations')
            else:
                plt.plot(t_eval, mc_outputs[mcind,:], 'k', linewidth=1, alpha=0.1)
            ymax = max(ymax, plot_spread(mc_outputs[mcind,:])[1])
        plt.xlim([t0, tmax])
        plt.ylim([0, ymax])
        plt.legend(frameon=False, bbox_to_anchor=(2.3, 1))
        plt.ylabel('algae $P$')
        plt.tick_params(labelbottom=False)
        plt.title('(a)')

        plt.subplot(312)
        cv_line = plt.plot(t_eval, stdev/mean, color='k', label='CV')
        plt.ylabel('CV', color='k')
        plt.xlim((t0, tmax))
        plt.ylim(plot_spread(stdev/mean))
        plt.gca().tick_params(axis='y', labelcolor='k')
        plt.tick_params(labelbottom=False)

        ax2 = plt.gca().twinx()
        skew_line = ax2.plot(t_eval, skewness, 'k', linestyle='--', label='skewness')
        ax2.set_ylabel('Skewness', color='k')
        ax2.set_ylim(plot_spread(skewness))
        ax2.tick_params(axis='y', labelcolor='k')
        plt.title('(b)')

        lines = cv_line + skew_line
        labs = [l.get_label() for l in lines]
        ax2.legend(lines, labs, frameon=False)

        plt.subplot(313)
        plt.stackplot(t_eval, global_sensitivity.T, labels=labels)
        plt.xlim((t0, tmax))
        plt.ylim((0, 1))
        plt.legend(frameon=False, bbox_to_anchor=(1.3, 1))
        plt.ylabel('Global sensitivity $S_I$')
        plt.xlabel('Time $t$ [days]')
        plt.title('(c)')

        plt.subplots_adjust(right=0.802, hspace=0.28)

# Now row-standardize y axes on CV + skewness plots
plt.figure(2)
for row in range(len(N_cvs)):
    skew_ymin, cv_ymin = np.Inf, np.Inf
    skew_ymax, cv_ymax = -np.Inf, -np.Inf
    for col in range(len(P_cvs)):
        plot_index = 3*row + col
        #plt.subplot(3, 3, plot_index)
        cv_ymin = min(cv_ymin, cv_axes[plot_index].get_ylim()[0])
        cv_ymax = max(cv_ymax, cv_axes[plot_index].get_ylim()[1])
        skew_ymin = min(skew_ymin, skew_axes[plot_index].get_ylim()[0])
        skew_ymax = max(skew_ymax, skew_axes[plot_index].get_ylim()[1])
    for col in range(len(P_cvs)):
        plot_index = 3*row + col
        #plt.subplot(3, 3, plot_index)
        cv_axes[plot_index].set_ylim([cv_ymin, cv_ymax])
        skew_axes[plot_index].set_ylim([skew_ymin, skew_ymax])

# Now row+col-standardize y axes on mean + stdev plot
plt.figure(1)
ymin = np.Inf
ymax = -np.Inf
for row in range(3):
    for col in range(3):
        plot_index = 3*row + col + 1 
        plt.subplot(3, 3, plot_index)
        ymin = min(ymin, plt.ylim()[0])
        ymax = max(ymax, plt.ylim()[1])

for row in range(3):
    for col in range(3):
        plot_index = 3*row + col + 1 
        plt.subplot(3, 3, plot_index)
        plt.ylim([ymin, ymax])


plt.figure(1)
plt.suptitle('Figure S3')
plt.subplots_adjust(right=0.802, hspace=0.06, wspace=0.06)

plt.figure(2)
plt.suptitle('Figure S4')
plt.subplots_adjust(right=0.802, hspace=0.06, wspace=0.06)

plt.figure(3)
plt.suptitle('Figure S5')
plt.subplots_adjust(right=0.802, hspace=0.06, wspace=0.06)

plt.figure(4)
plt.suptitle('Figure 5')

# plt.interactive(True)
# plt.show(block=True)