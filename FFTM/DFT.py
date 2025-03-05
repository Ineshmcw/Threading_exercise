import math

N = 4
pi = 3.14159265358979323846
x = [1,2,3,4]
X = [0] * N

for k in range(0,N):
    for n in range(0,N):
        angle = (-2*pi*k*n)/N
        X[k] += x[n]*complex(math.cos(angle), math.sin(angle))

print(X)


