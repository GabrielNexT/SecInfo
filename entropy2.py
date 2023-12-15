from scipy.stats import entropy, combine_pvalues
import numpy as np
import pandas as pd

from rc4 import (
    ksa,
    meksa,
    prga,
    meprga,
    ApEn,
    rec4a_prga,
    ideal_rc4_prga,
    vmpc_prga,
    rc4plus_prga,
)
from metrics import keystream_ap_en, keystream_run_test, keystream_cumsum, fisher_method

KEYS = list(
    map(
        lambda x: x.encode("utf8"),
        [
            "ETHAR546",
            "METHAQTG13",
            "THEKRAABASS3345",
            "AD56790877HNKUNNHG6",
            "FVHHGFHG64877HYNG9TG",
        ],
    )
)



SAMPLES = 2**10


def test_keys():
    df = pd.DataFrame(columns=['key', 'cipher', 'ApEn', 'run_test', 'cumsum', 'entropy', 'random'])
    for key in KEYS:
        print("Key: {}".format(key))

        gen_rc4 = np.array([x for (_, x) in zip(range(SAMPLES), prga(ksa(key)))])
        gen_rc4a = np.array([x for (_, x) in zip(range(SAMPLES), rec4a_prga(key, key))])
        gen_merc4 = np.array([x for (_, x) in zip(range(SAMPLES), meprga(meksa(key)))])
        gen_ideal_rc4 = np.array([x for (_, x) in zip(range(SAMPLES), ideal_rc4_prga(key))])
        gen_vmpc = np.array([x for (_, x) in zip(range(SAMPLES), vmpc_prga(key))])
        gen_rc4plus = np.array([x for (_, x) in zip(range(SAMPLES), rc4plus_prga(key))])


        tests = dict(ApEn = keystream_ap_en, run_test = keystream_run_test,
                     cumsum = keystream_cumsum)
        ciphers = dict(rc4 = gen_rc4, rc4a = gen_rc4a, merc4 = gen_merc4,
                       ideal_rc4 = gen_ideal_rc4, vmpc = gen_vmpc, rc4plus =
                       gen_rc4plus)

        
        for (cname, cipher) in ciphers.items():
            mask = (df['key'] == key) & (df['cipher'] == cname)
            h = entropy(cipher)
            if df[mask].empty:
                df.loc[len(df)] = {'key': key, 'cipher': cname, 'entropy': h, 'random': 0}
                mask = (df['key'] == key) & (df['cipher'] == cname)
            rval = 0
            for (tname, test) in tests.items():
                res = test(cipher)
                accept = 1 if res > 0.01 else 0
                acceptstr = 'random' if accept else '~random'
                if accept:
                    rval += 1
                df.loc[mask, [tname, 'entropy', 'random']] = [res, h, rval]
    return df


def main():
    mcols = ['ApEn', 'run_test', 'cumsum', 'entropy']
    df = test_keys()
    print(df)
    #df.to_csv('results.csv', index=False)
    #print(df.groupby('cipher')[mcols].agg(['mean', 'std']))
    #df.to_latex('results.tex', index=False, float_format="%.3f")
    #print(df.groupby('cipher')[mcols[:-1]].agg(fisher_method))
    fisher = df.groupby('cipher')[mcols[:-1]].agg(lambda x: combine_pvalues(x, method='fisher')[1])
    fisher.to_latex('fisher.tex', float_format="%.3f")

if __name__ == "__main__":
    main()
