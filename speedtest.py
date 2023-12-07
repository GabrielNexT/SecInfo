from csv import DictReader
import time

from rc4 import ideal_rc4, rc4, rc4_plus, rc4a, merc4, vmpc

FILE_NAME = "chat_dataset.csv"

KEYS = list(map(lambda x: x.encode('utf8'),
                ['ETHAR546', 'METHAQTG13', 'THEKRAABASS3345',
                 'AD56790877HNKUNNHG6', 'FVHHGFHG64877HYNG9TG']))
SAMPLES = 2**10

def load_messages_list() -> list[bytes]:
    messages = []
    with open(FILE_NAME) as file:
        reader = DictReader(file)
        for row in reader:
            messages.append(row["message"].encode('utf8'))
    print(f"{len(messages)} loaded from {FILE_NAME}")
    return messages

def test_key(key: str, messages: list[bytes]):
    print(f"{key=}")

    start = time.time()
    for m in messages:
        rc4(key, m)
    end = time.time()
    print(f"Time elapse for RC4 {end-start:.6f} second(s)")

    start = time.time()
    for m in messages:
        merc4(key, m)
    end = time.time()
    print(f"Time elapse for MERC4 {end-start:.6f} second(s)")

    start = time.time()
    for m in messages:
        rc4a(key, m)
    end = time.time()
    print(f"Time elapse for RC4A {end-start:.6f} second(s)")

    start = time.time()
    for m in messages:
        vmpc(key, m)
    end = time.time()
    print(f"Time elapse for VMPC {end-start:.6f} second(s)")

    start = time.time()
    for m in messages:
        ideal_rc4(key, m)
    end = time.time()
    print(f"Time elapse for Ideal RC4 {end-start:.6f} second(s)")
    
    start = time.time()
    for m in messages:
        rc4_plus(key, m)
    end = time.time()
    print(f"Time elapse for RC4+ {end-start:.6f} second(s)")

def main():
    messages = load_messages_list()
    for key in KEYS:
        print()
        test_key(key, messages)
        print()

if __name__ == '__main__':
    main()
