from scipy.stats import hypergeom
import time
import numpy as np

#[M, n, N] = [60e6, 30000, 30e6]
[M, n, N] = [60e6, 5, 30e6]
x = np.arange(0, n+1)

rv = hypergeom(M, n, N)
p_hyper = rv.pmf(x)

list_prob = []
for i in range(0, len(x)):
    list_prob.append(float(sum(p_hyper[i:])))

print (x)
list_prob = np.array(list_prob)
print (list_prob);
expect_depth = int(max(x[list_prob >= 0.8]))
print (expect_depth)
prob_expect_depth = round(float(list_prob[x == expect_depth]), 2)
print (prob_expect_depth)
