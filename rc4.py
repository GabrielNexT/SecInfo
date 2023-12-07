import numpy as np
import random


def rec4a_prga(key1: bytes, key2: bytes):
    s1 = list(range(256))
    s2 = list(range(256))

    j = 0

    for i in range(256):
        j = (j + s1[i] + key1[i % len(key1)]) % 256
        s1[i], s1[j] = s1[j], s1[i]

    for i in range(256):
        j = (j + s1[i] + key2[i % len(key2)]) % 256
        s2[i], s2[j] = s2[j], s2[i]

    i = j1 = j2 = 0
    while True:
        i = (i + 1) % 256
        j1 = (j1 + s1[i]) % 256
        s1[i], s1[j1] = s1[j1], s1[i]
        ks1 = s2[(s1[i] + s1[j1]) % 256]

        j2 = (j2 + s2[i]) % 256
        s2[i], s2[j2] = s2[j2], s2[i]
        ks2 = s1[(s2[i] + s2[j2]) % 256]

        yield ks1 ^ ks2


def compute_hash(seed: bytes):
    p = 31
    mod = int(1e9 + 9)
    hash = 0
    pow = 1
    for c in seed:
        hash = (hash + (c - ord("a") * pow)) % mod
        pow = (pow * p) % mod
    return hash


def ideal_rc4_prga(key: bytes):
    s = list(range(256))

    rand = random.Random(compute_hash(key))

    i = 0
    while True:
        j = rand.randint(0, 255)
        s[i], s[j] = s[j], s[i]
        old_i = i
        i = (i + 1) % 256
        yield (s[old_i] + s[j]) % 256


def vmpc_prga(key: bytes):
    s = list(range(256))
    l = len(key)
    j = 0

    for m in range(767):
        n = m % 256
        j = s[(j + s[n] + key[m % l]) % 256]
        s[n], s[j] = s[j], s[n]

    i = j = 0
    while True:
        a = s[i % 256]
        j = s[(j + a) % 256]
        ks = s[(s[s[j]] + 1) % 256]
        s[i], s[j] = s[j], s[i]
        i = (i + 1) % 256
        yield ks


def rc4plus_prga(key: bytes):
    s = list(range(256))
    l = len(key)
    j = 0

    for i in range(256):
        j = (j + s[i] + key[i % l]) % 256
        s[i], s[j] = s[j], s[i]

    i = j = 0
    while True:
        i = (i + 1) % 256
        a = s[i]
        j = (j + a) % 256
        b = s[j]
        s[i], s[j] = s[j], s[i]

        c = (s[(i << 5 ^ j >> 3) % 256] + s[(j << 5 ^ i >> 3) % 256]) % 256
        yield ((s[(a + b) % 256] + s[(c ^ 0xAA) % 256]) ^ s[(j + b) % 256]) % 256


def ksa(key):
    kl = len(key)
    S = np.arange(256, dtype=np.int32)
    j = 0
    for i in range(256):
        j = (j + S[i] + key[i % kl]) % 256
        S[i], S[j] = S[j], S[i]
    return S


def prga(s):
    i = 0
    j = 0
    while True:
        i = (i + 1) % 256
        j = (j + s[i]) % 256
        s[i], s[j] = s[j], s[i]
        # print(f"{s[i]=} + {s[j]=} = {s[i] + s[j]}")
        yield s[(s[i] + s[j]) % 256]


def meksa(key):
    S = np.arange(256, dtype=np.int32)
    j = 0
    for i in range(256):
        ic = 255 - S[i]
        j = (j + S[ic] + key[i % len(key)]) % 256
        S[i], S[j] = S[j], S[i]
    return S


def meprga(s):
    i = j1 = j2 = 0
    while True:
        i = (i + 1) % 256
        j1 = (j1 + s[i]) % 256
        s[i], s[j1] = s[j1], s[i]
        j2 = (j2 + s[i]) % 256
        s[i], s[j2] = s[j2], s[i]
        yield s[(s[i] + s[j1] + s[j2]) % 256] ^ j2


def rc4(key, plaintext):
    keystream = prga(ksa(key))
    key_char = zip(keystream, plaintext)
    return bytes([ks ^ c for (ks, c) in key_char])


def rc4a(key, plaintext):
    key2 = "".join(random.sample(str(key), len(key))).encode("utf8")
    keystream = rec4a_prga(key, key2)
    key_char = zip(keystream, plaintext)
    return bytes([ks ^ c for (ks, c) in key_char])


def merc4(key, plaintext):
    keystream = meprga(meksa(key))
    key_char = zip(keystream, plaintext)
    return bytes([ks ^ c for (ks, c) in key_char])


def vmpc(key, plaintext):
    keystream = vmpc_prga(key)
    key_char = zip(keystream, plaintext)
    return bytes([ks ^ c for (ks, c) in key_char])


def ideal_rc4(key, plaintext):
    keystream = ideal_rc4_prga(key)
    key_char = zip(keystream, plaintext)
    return bytes([ks ^ c for (ks, c) in key_char])


def rc4_plus(key, plaintext):
    keystream = rc4plus_prga(key)
    key_char = zip(keystream, plaintext)
    return bytes([ks ^ c for (ks, c) in key_char])


def ApEn(U, m, r) -> float:
    """Approximate_entropy."""

    def _maxdist(x_i, x_j):
        return max([abs(ua - va) for ua, va in zip(x_i, x_j)])

    def _phi(m):
        x = [[U[j] for j in range(i, i + m - 1 + 1)] for i in range(N - m + 1)]
        C = [
            len([1 for x_j in x if _maxdist(x_i, x_j) <= r]) / (N - m + 1.0)
            for x_i in x
        ]
        return (N - m + 1.0) ** (-1) * sum(np.log(C))

    N = len(U)

    return abs(_phi(m + 1) - _phi(m))
