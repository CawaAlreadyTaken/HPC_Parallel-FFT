from random import randint

LOG2_DEGREE = 18
MAX_COEFFICIENT_MODULO = 100

OUTPUT_FILE = "../data/dataset_1_1.txt"

degree = pow(2, LOG2_DEGREE)
print(f"Creating input type 1 of degree {degree}...")
with open(OUTPUT_FILE, 'w') as f:
    polynomial = [randint(-MAX_COEFFICIENT_MODULO, MAX_COEFFICIENT_MODULO) for _ in range(degree)]
    f.write(str(degree) + "\n")
    f.write(" ".join([str(x) for x in polynomial]))
    f.write("\n")
    polynomial = [randint(-MAX_COEFFICIENT_MODULO, MAX_COEFFICIENT_MODULO) for _ in range(degree)]
    f.write(str(degree) + "\n")
    f.write(" ".join([str(x) for x in polynomial]))
    f.write("\n")

print(f"[*] Done. Dataset saved in {OUTPUT_FILE}")
