import sys
import math
import numpy as np

print("Float case: \n\n")
epsilon = np.float32(1.0)
zero = np.float32(0.0)
one = np.float32(1.0)
two = np.float32(2.0)

# print("Double case: \n\n")
# epsilon = np.float64(1.0)
# zero = np.float64(0.0)
# one = np.float64(1.0)
# two = np.float64(2.0)

while one + epsilon/two > one:
    epsilon /= two

bits = -int(math.log2(epsilon))

print(f"Epsilon: {epsilon}\n")
print(f"Bits in mantissa: {bits} bits + one bit for sign\n")
print(f"Bits in exponent: {64 - two - bits} bits + one bit for sign\n")

value = one
max_exponent = zero
while not math.isinf(value):
    value *= two
    max_exponent += one

value = one
min_exponent = zero
while value != zero:
    value /= two
    min_exponent -= one

print(f"Max exponent: {max_exponent}\n")
print(f"Min exponent: {min_exponent}\n")

if (one == one + epsilon/two):
    print("one == one + epsilon/2")
if (one < one + epsilon):
    print("one < one + epsilon")
if (one < one + epsilon + epsilon/two):
    print("one < one + epsilon + epsilon/2")
if (one + epsilon/two < one + epsilon):
    print("one + epsilon/2 < one + epsilon")
if (one + epsilon/two < one + epsilon + epsilon/two):
    print("one + epsilon/2 < one + epsilon  + epsilon/2")
if (one + epsilon < one + epsilon + epsilon/two):
    print("one + epsilon < one + epsilon  + epsilon/2")