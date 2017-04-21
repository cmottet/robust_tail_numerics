
#
# The goal of this function is to fit a mixture of Shifted Pareto and 2-PLT
# that is feasible for program (1), i.e. we find two distribution functions
# F1 and F2 such that
#
# w1 + w2 = 1
# w1 f1(a) + w2 f2(a) = eta
# - w1 f1'(a) - w2 f2'(a) = nu
# w1 (1 - F1(a)) + w2 (1- F2(a)) = beta
#
# where
#  * F1(x) is a two point mass distribution function (2-PLT)
#  * F2 is a pareto distribution, i.e. F2(x) = 1 - (scale/x)^shape for all x => scale
#  and the shape is a number between 0 and 2
#
# The problem of finding such distribution functions is equivalent to finding
# F2(x) such that
#
# eta1 = f1(a) = (eta - w2 f2(a))/w1 = (eta - w2 eta2(a))/w1  => 0
# nu1 = -f1'(a) = (nu + w2 f2'(a))/w1 = (nu - w2 nu2(a))/w1 => 0
# beta1 =  (beta - w2 (1 - F2(a)))/w1  = (beta - w2 beta2)/w1 => 0
#
# and
#
# mu1**2 = (eta1/nu1)**2 <= sigma1 = 2*beta1/nu1
#
# If there exist a pareto distribution F2(x) such that the 3 previous
# equalities and the inequality hold, then there exists a 2-PLT distribution whose
# key parameters (x1,x2) and (p1,p2), as given in Equation (6), can
# be obtained using the function "getDistribution" available in the package "RobustTail"
#
# The process is not deterministic, it requires sampling some parameters
# from uniform distributions. In a reproducibility perspective,
# a "seed" parameter can be passed to the function to set the seed when sampling
# from the uniform distributions
#
# we now test our function
#
# from scipy.stats import expon
# a = expon.ppf(0.7)
# eta = expon.pdf(a)
# nu = expon.pdf(a)
# beta = 1 - expon.cdf(a)
#
# mixture = fitMixturePareto2PLT(a,nu,eta,beta)
#
# # Check the system
# with(mixture, data.frame( nu = with(Pareto, w*nu) + with(PLT, w*nu),
#                           eta = with(Pareto, w*eta) + with(PLT, w*eta),
#                           beta = with(Pareto, w*beta) + with(PLT, w*beta)))
#
# # Check that eta1, beta1, nu1, eta2, beta2, nu2 >= 0
# with(mixture, data.frame( nu = with(Pareto, nu>= 0) & with(PLT, nu >= 0),
#                           eta =with(Pareto, eta >= 0) & with(PLT, eta >= 0),
#                           beta = with(Pareto, beta>= 0) & with(PLT, beta >= 0)))
#
import numpy as np
import robust_tail
import distribution_pty
from scipy.stats import uniform
from math import sqrt


def fit_mixture_pareto_2PLT(a, nu, eta, beta, iter_max=10**4, seed=None):
    mu = eta/nu
    sigma = 2*beta/nu

    for i in range(iter_max):
        # Assign randomly w2 and the shape of F2(x)
        if seed is not None:
            np.random.seed(seed)

        w2 = uniform.rvs(loc=0, scale=1)

        if seed is not None:
            np.random.seed(seed)

        shape = uniform.rvs(loc=0, scale=2)
        w1 = (1 - w2)

        # Compute the necessary (but not sufficient) bounds on the X =
        # (x/scale)**(-shape) of F2(x) for there existence of F1(x)

        # x_ub is an upper bound on X that ensures that
        # eta1, nu1, and beta1 => 0
        x_ub = min([1, beta/w2, eta*a/(w2*shape), nu*a**2/(w2*shape*shape + 1)])

        # we now compute XlB, a necessary lower bound for mu1^1 <= sigma1
        b = mu*a*shape - (shape + 1)*shape*sigma/2 - a**2
        fourac = shape*(shape + 2)*a**2*(sigma - mu**2)
        delta = b**2 - fourac

        x_lb = np.where(delta >= 0, max([nu/(4*w2)*(b + sqrt(delta)), 0]), 0)

        # Get a random value of X respecting the bounds
        # and derive the value of the scale parameter in F2(x)
        if seed is not None:
            np.random.seed(seed)

        x = uniform.rvs(loc=x_lb, scale=x_ub)
        scale = x**(1/shape)*a

        # We now need to compute eta1, beta1, and nu1
        # to check if there exists a 2-PLT F1(x) for
        # this specific F2(x)
        eta2 = distribution_pty.pareto.dcdf(x=a, d=1, scale=scale, b=shape)
        nu2 = -distribution_pty.pareto.dcdf(x=a, d=2, scale=scale, b=shape)
        beta2 = 1 - distribution_pty.pareto.dcdf(x=a, d=0, scale=scale, b=shape)

        eta1 = (eta - w2*eta2)/w1
        nu1 = (nu - w2*nu2)/w1
        beta1 = (beta - w2*beta2)/w1

        mu1 = eta1/nu1
        sigma1 = 2*beta1/nu1

        # Check the necesserary condition for the existence
        # of F1(x)
        if mu1**2 <= sigma1:
            break

        if i > iter_max:
           return "No solution"

    P = robust_tail.get_distribution(mu1, sigma1)

    pareto = {'scale': scale, 'shape': shape, 'eta': eta2, 'nu': nu2, 'beta': beta2, 'w': w2}
    PLT = {'x': P.x, 'p': P.p, 'w': w1, 'nu': nu1, 'eta': eta1, 'beta': beta1}

    output = {'Pareto': pareto, 'PLT': PLT}

    return output
