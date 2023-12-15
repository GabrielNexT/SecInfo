from collections import Counter
from scipy.special import gammainc
import math
import numpy as np
import scipy.integrate as integrate
from functools import cache
from scipy.stats import norm

def bytes_to_bits(b):
    return ''.join(format(byte, '08b') for byte in b)


def approximate_entropy(e, m=None):
    n = len(e)

    if m is None:
        m = int(math.floor(math.log(n, 2) - 5))

    def phi(c):
        return sum(map(lambda x: x * math.log(x), c.values()))

    def get_counts(e, m):
        augmented_e = e + e[0:m-1]
        an = len(augmented_e)
        patterns = (augmented_e[i:i+m] for i in range(an-m+1))
        counts = Counter(patterns)
        countsv = dict(map(lambda x: (x[0],x[1]/n), counts.items()))
        return countsv

    ap_en = phi(get_counts(e, m)) - phi(get_counts(e, m+1))
    chi_squared = 2 * n * (math.log(2) - ap_en)
    pval = 1 - gammainc(2**(m-1), chi_squared/2)
    #return ap_en, chi_squared, pval
    return pval

def keystream_ap_en(ks):
    return approximate_entropy(bytes_to_bits(ks))


def cumsum(e, mode=0):
    def nrange(start, stop):
        if start is not int:
            start = math.floor(start)
        if stop is not int:
            stop = math.floor(stop)
        return range(start, stop+1)
    n = len(e)
    if mode == 1:
        e = reversed(e)
    vals = list(map(lambda x: 1 if int(x) else -1, e))
    sums = np.cumsum(vals)
    z = max(np.abs(sums))
    sum1 = 0
    for k in nrange((-n/z+1)/4, (n/z-1)/4):
        sum1 += normal_cdf((4*k+1)*z/np.sqrt(n)) - normal_cdf((4*k-1)*z/np.sqrt(n))
    sum2 = 0
    for k in nrange((-n/z-3)/4, (n/z-1)/4):
        sum2 += normal_cdf((4*k+3)*z/np.sqrt(n)) - normal_cdf((4*k+1)*z/np.sqrt(n))
    pval = 1 - sum1 + sum2
    #return z, pval
    return pval
        
def keystream_cumsum(ks):
    return cumsum(bytes_to_bits(ks))

@cache
def normal_cdf(z):
    t = 1/np.sqrt(2*np.pi) 
    return norm.cdf(z)
#   val, err = integrate.quad(lambda x: np.exp(-x**2/2), -np.inf, z)
#   print((z, val, err))
    return t * val

def keystream_cumsum(ks):
    return cumsum(bytes_to_bits(ks))

def monobit_test(e):
    n = len(e)
    s = sum(map(lambda x: 1 if int(x) else -1, e))
    sobs = abs(s)/np.sqrt(n)
    pval = math.erfc(sobs/np.sqrt(2))
    #return sobs, sobs, pval
    return pval

def run_test(e):
    n = len(e)
    tao = 2 / np.sqrt(n)
    prop = sum(map(lambda x: 1 if int(x) else 0, e))/n
    if (np.abs(prop - 0.5) >= tao):
        return prop, 0, 0
    vobs = 1
    for k in range(0, n-1):
        v = 1 if e[k] != e[k+1] else 0
        vobs += v
    pval = math.erfc(abs(vobs - 2*n*prop*(1-prop))/(2*np.sqrt(2*n)*prop*(1-prop)))
    #return prop, vobs, pval
    return pval

def keystream_run_test(ks):
    return run_test(bytes_to_bits(ks))

def fisher_method(pvals):
    if pvals is np.array:
        return -2 * np.sum(np.log(pvals))
    return -2 * sum(map(lambda x: math.log(x), pvals))
