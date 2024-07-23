import time
import numpy as np

N = 40000
a = 100

A = np.random.rand(N, a, a)
B = np.random.rand(N, a, a)
C = np.zeros((N, a, a))

start = time.time()
C = A@B
end = time.time()
print("Time taken for numpy matrix multiplication: ", end-start)
start = time.time()
C *= 1
end = time.time()
print("Time taken for numpy matrix multiplication: ", end-start)


