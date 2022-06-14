import lamberthub
import numpy as np
import timeit

NUM_TRIALS = 10_000
print([solver.__name__ for solver in lamberthub.ALL_SOLVERS])
print(lamberthub.__version__)

r1 = (np.random.random(3 * NUM_TRIALS).reshape(NUM_TRIALS, 3) - 0.5) * 1.0e11
r2 = (np.random.random(3 * NUM_TRIALS).reshape(NUM_TRIALS, 3) - 0.5) * 1.0e11
tof = np.random.random(NUM_TRIALS) * 1.0e8
v1s = np.zeros(3 * NUM_TRIALS).reshape(NUM_TRIALS, 3)
v2s = np.zeros(3 * NUM_TRIALS).reshape(NUM_TRIALS, 3)

def benchmark_gooding():
    global v1s, v2s
    for i in range(NUM_TRIALS):
        v1, v2 = lamberthub.gooding1990(1.327e20, r1[i], r2[i], tof[i], prograde=False, low_path=True, maxiter=100, atol=1e-5, rtol=1e-7, full_output=False)
        v1s[i] = v1
        v2s[i] = v2

def benchmark_avanzini():
    global v1s, v2s
    for i in range(NUM_TRIALS):
        v1, v2 = lamberthub.avanzini2008(1.327e20, r1[i], r2[i], tof[i], prograde=False, low_path=True, maxiter=100, atol=1e-5, rtol=1e-7, full_output=False)
        v1s[i] = v1
        v2s[i] = v2

def benchmark_arora():
    global v1s, v2s
    for i in range(NUM_TRIALS):
        v1, v2 = lamberthub.arora2013(1.327e20, r1[i], r2[i], tof[i], prograde=False, low_path=True, maxiter=100, atol=1e-5, rtol=1e-7, full_output=False)
        v1s[i] = v1
        v2s[i] = v2

def benchmark_izzo():
    global v1s, v2s
    for i in range(NUM_TRIALS):
        v1, v2 = lamberthub.izzo2015(1.327e20, r1[i], r2[i], tof[i], prograde=False, low_path=True, maxiter=100, atol=1e-5, rtol=1e-7, full_output=False)
        v1s[i] = v1
        v2s[i] = v2

time = timeit.Timer(benchmark_gooding).timeit(number=1)
print(time, time / NUM_TRIALS)
time = timeit.Timer(benchmark_avanzini).timeit(number=1)
print(time, time / NUM_TRIALS)
#time = timeit.Timer(benchmark_arora).timeit(number=1)
#print(time, time / NUM_TRIALS)
time = timeit.Timer(benchmark_izzo).timeit(number=1)
print(time, time / NUM_TRIALS)