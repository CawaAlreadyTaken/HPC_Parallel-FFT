from random import randint

LOG2_DEGREE = 29
MAX_COEFFICIENT_MODULO = 100

OUTPUT_FILE = "../data/dataset_1_5.txt"

degree = pow(2, LOG2_DEGREE)
print(f"Creating input type 1 of degree {degree}...")

with open(OUTPUT_FILE, 'w') as f:
    f.write(str(degree) + "\n")
    for _ in range(degree):
        coefficient = randint(-MAX_COEFFICIENT_MODULO, MAX_COEFFICIENT_MODULO)
        f.write(str(coefficient) + " ")
    f.write("\n")

    # Genera e scrivi il secondo polinomio
    f.write(str(degree) + "\n")
    for _ in range(degree):
        coefficient = randint(-MAX_COEFFICIENT_MODULO, MAX_COEFFICIENT_MODULO)
        f.write(str(coefficient) + " ")
    f.write("\n")

print(f"[*] Done. Dataset saved in {OUTPUT_FILE}")
