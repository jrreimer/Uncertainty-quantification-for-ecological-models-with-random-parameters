import numpy as np

def rho_range(m, cv):
    """
    Returns the lower- and upper- bounds for the correlation coefficient
    of a bivariate lognormal random variable with bivariate mean m and
    bivariate coefficient of variation cv.
    """

    # Bounds for rho
    L = np.exp(np.sqrt(np.log(1 + cv[0]**2) * np.log(1 + cv[1]**2)))
    rho_upper = (L-1)/(cv[0]*cv[1])
    rho_lower = (1/L-1)/(cv[0]*cv[1])

    return rho_lower, rho_upper

def normal_parameters_from_lognormal(m, cv, rho):
    """
    Returns the mean mu and covariance cov of a bivariate normal
    distribution for random variables (Y1, Y2) such that (X1,X2) =
    (exp(Y1), exp(Y2)) is a bivariate lognormal random variable whose
    mean is ms, coefficient of variation is cvs, and correlation(X1,X2)
    is rho.
    """

    # Bounds for rho
    rho_lower, rho_upper = rho_range(m, cv)

    assert rho_lower <= rho, "Correlation rho is too small"
    assert rho <= rho_upper, "Correlation rho is too large"

    mu = np.log(m/np.sqrt(1+cv**2))
    cov = np.zeros([2,2])

    cov[0,0] = np.log(1 + cv[0]**2)
    cov[1,1] = np.log(1 + cv[1]**2)
    cov[0,1] = np.log(1 + rho*cv[0]*cv[1])
    cov[1,0] = cov[0,1]

    return mu, cov

def plot_spread(data):

    minval, maxval = np.min(data), np.max(data)
    spread = maxval-minval
    return (minval-spread/10, maxval+spread/10)

def compute_batch_skewness(generator, mean, stdev, ensemble_size=int(1e7), batch_size=int(1e4), axis=0):
    """
    Computes sample skewness along a particular axis using batch sampling.
    Assumes the generator generates an array with the ensemble along axis 0.
    """

    # For now we'll do this the simple and probably dumb way with cumulative
    # sums
    cumsum = 0
    m = 0

    while m < ensemble_size:
        batch_size = min(batch_size, ensemble_size-m)
        samples = (generator(batch_size)-mean)/stdev

        cumsum = (m*cumsum + np.sum(samples**3, axis=0))/(m+batch_size)

        m += batch_size

    return cumsum

    #N = array.shape[axis]
    #mean = np.mean(array, axis=axis)
    #stdev = np.std(array, axis=axis)
    #centered_array = array - mean

    #skewness = (1/N * np.sum(centered_array**3))/(stdev**3)
    #return skewness
