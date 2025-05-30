from scipy.stats import f, ncf
import numpy as np

def main(stat_viewpoint, power, effect_size, k, alpha, eps=1, rho=0, pop_size=None):
    if stat_viewpoint == "Enumerative":
        if pop_size is None:
            raise ValueError("pop_size must be provided for the Enumerative viewpoint.")
        result = fpc(pop_size, sample_size(power, effect_size, k, alpha, eps, rho))
        print(f"Enumerative sample size (with FPC): {result}")
    elif stat_viewpoint == "Generalised":
        result = sample_size(power, effect_size, k, alpha, eps, rho)
        print(f"Generalised sample size: {result}")
    else:
        raise ValueError("stat_viewpoint must be either 'Enumerative' or 'Generalised'.")
    return result


def power_calc_RMANOVA(n, effect_size, k, alpha, eps = 1, rho = 0):
    F_c = f.ppf(1-alpha, k-1, (k)*(n-1))
    #print(F_c)
    #print(f.cdf(3.002, k-1, (k)*(n-1)))
    #lambda_ = effect_size**2*n*k*eps/(1-rho)
    #print(lambda_)
    power = 1-ncf.cdf(F_c, k-1, (k)*(n-1), effect_size**2*n*k*eps/(1-rho))
    return power


def sample_size(power, effect_size, k, alpha, eps=1, rho=0):
    n = 2
    while True:
        power_result = power_calc_RMANOVA(n, effect_size, k, alpha, eps, rho)
        if power_result >= power:
            return n
        n += 1


def fpc(pop_size, samp_size):
    return int(np.ceil(samp_size / (1 + (samp_size - 1) / pop_size)))
