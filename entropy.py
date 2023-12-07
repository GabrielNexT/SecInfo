from scipy.stats import entropy
import numpy as np

from rc4 import ksa, meksa, prga, meprga, ApEn, rec4a_prga, ideal_rc4_prga, vmpc_prga, rc4plus_prga
KEYS = list(map(lambda x: x.encode('utf8'),
                ['ETHAR546', 'METHAQTG13', 'THEKRAABASS3345',
                 'AD56790877HNKUNNHG6', 'FVHHGFHG64877HYNG9TG']))
SAMPLES = 2**10

def test_key(key: str):
    print('Key: {}'.format(key))

    gen_rc4 = np.array([x for (_,x) in zip(range(SAMPLES), prga(ksa(key)))])
    gen_rc4a = np.array([x for (_,x) in zip(range(SAMPLES), rec4a_prga(key, key))])
    gen_merc4 = np.array([x for (_,x) in zip(range(SAMPLES), meprga(meksa(key)))])
    gen_ideal_rc4 = np.array([x for (_,x) in zip(range(SAMPLES), ideal_rc4_prga(key))])
    gen_vmpc = np.array([x for (_,x) in zip(range(SAMPLES), vmpc_prga(key))])
    gen_rc4plus = np.array([x for (_,x) in zip(range(SAMPLES), rc4plus_prga(key))])

    print('Entropy of RC4: {}'.format(entropy(gen_rc4, base=2)))
    print('Entropy of MERC4: {}'.format(entropy(gen_merc4, base=2)))
    print('Entropy of RC4A: {}'.format(entropy(gen_rc4a, base=2)))
    print('Entropy of Ideal RC4: {}'.format(entropy(gen_ideal_rc4, base=2)))
    print('Entropy of VMPC: {}'.format(entropy(gen_vmpc, base=2)))
    print('Entropy of RC4+: {}'.format(entropy(gen_rc4plus, base=2)))
    print()
    print('ApEn of RC4: {}'.format(ApEn(gen_rc4, 2, 3)))
    print('ApEn of MERC4: {}'.format(ApEn(gen_merc4, 2, 3)))
    print('ApEn of RC4A: {}'.format(ApEn(gen_rc4a, 2, 3)))
    print('ApEn of Ideal RC4: {}'.format(ApEn(gen_ideal_rc4, 2, 3)))
    print('ApEn of VMPC: {}'.format(ApEn(gen_vmpc, 2, 3)))
    print('ApEn of RC4+: {}'.format(ApEn(gen_rc4plus, 2, 3)))


def main():
    for key in KEYS:
        test_key(key)
        print()

if __name__ == '__main__':
    main()
