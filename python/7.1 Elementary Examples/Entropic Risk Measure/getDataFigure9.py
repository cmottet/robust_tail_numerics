# Clear the sytem of all paste variables
def clear_all():
    all = [var for var in globals() if "__" not in (var[:2], var[-2:])]
    for var in all:
        del globals()[var]

clear_all()


from distribution_pty import expon
import robust_tail
import numpy as np
import pandas as pd
import pickle
from math import (exp, log)



#
# Parameters of program 5
#
scale = 1
a = expon.ppf(q=0.7, scale=scale)
beta = 1 - expon.dcdf(x=a, d=0, scale=scale)
eta = expon.dcdf(x=a, d=1, scale=scale)
nu = -expon.dcdf(x=a, d=2, scale=scale)

mu = eta/nu
sigma = 2*beta/nu

#
# Compute the optimal upper bound for different values of theta
#

Theta = np.linspace(start=0.04, stop=2, num=50)


def run_func(theta):
    # Compute the first term in Equation (15)
    term1 = 1.0/(theta + 1)*(1 - exp(-(theta + 1)*a))

    # Compute the second term in Equation (15)
    H = lambda x: exp(-theta*a) / theta**2 * (theta*x + exp(-theta*x) - 1)
    bound_term2 = robust_tail.compute_bound(H=H, mu=mu, sigma=sigma, limsup=0, nu=nu)['bound'][0]

    # Compute the value of Equation (15)
    bound = 1.0/theta*log(term1 + bound_term2)

    output = pd.DataFrame({'theta': theta, 'term1': term1, 'bound_term2': bound_term2, 'bound': bound}, index=[0])
    return output


optim_bound = pd.DataFrame({'theta': [], 'term1': [], 'bound_term2': [], 'bound': []})
for theta in Theta:
    optim_bound = optim_bound.append(run_func(theta), ignore_index=True)


file = open("data/runEntropyExpDist.data", "w")
pickle.dump(optim_bound, file)
file.close()

