from scipy.stats import f, ncf
import numpy as np

def main():
    print(fpc(298, sample_size(0.9, 0.1, 3, 0.05, 1, 0)))

def power_calc_RMANOVA(n, effect_size, k, alpha, eps = 1, rho = 0):
    F_c = f.ppf(1-alpha, k-1, (k)*(n-1))
    #print(F_c)
    #print(f.cdf(3.002, k-1, (k)*(n-1)))
    #lambda_ = effect_size**2*n*k*eps/(1-rho)
    #print(lambda_)
    power = 1-ncf.cdf(F_c, k-1, (k)*(n-1), effect_size**2*n*k*eps/(1-rho))
    return power

def sample_size(power, effect_size, k, alpha, eps=1, rho=0):
    n = 1
    while True:
        power_result = power_calc_RMANOVA(n, effect_size, k, alpha, eps, rho)
        if power_result >= power:
            return n
        n += 1

def fpc(pop_size, samp_size):
    return int(np.ceil(samp_size/(1+(samp_size-1)/pop_size)))

main()