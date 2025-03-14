import math
import numpy as np
import time

N = 20
pi = 3.14159265358979323846
x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
X = [0] * N

def DFT():
    for k in range(0,N):
        for n in range(0,N):
            angle = (-2*pi*k*n)/N
            X[k] += x[n]*complex(math.cos(angle), math.sin(angle))

def FFT_NUMPY():
    X = np.fft.fft(x)

start_time = time.time()
DFT()
end_time = time.time()
execution_time1 = end_time - start_time
print(f"Execution time: {execution_time1:.10f} seconds")

start_time = time.time()
FFT_NUMPY()
end_time = time.time()
execution_time2 = end_time - start_time
print(f"Execution time: {execution_time2:.10f} seconds")

speedup = execution_time1 / execution_time2
print(f"Percentage improvement: {((speedup - 1) * 100):.2f}%")

"""
inesh@Mcw:~/Examples/Threading$ /usr/local/bin/python /home/inesh/Examples/Threading/FFTM/DFT.py
Execution time: 0.0003230572 seconds
Execution time: 0.0000326633 seconds
speed difference is :889.051094890511%
"""
