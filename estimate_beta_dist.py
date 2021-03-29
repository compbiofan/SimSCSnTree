# Given x_0 (which is the furthest point on a Lorenz curve to the diagonal), compute Alpha and Beta for the underlying Beta distribution of the number of reads in a bin. 

import math
from scipy.stats import beta
from scipy.optimize import newton_krylov
from scipy.optimize.nonlin import NoConvergence
import numpy as np

#IxAlpha = 0.4
#IxAlphap1 = 0.3
#total_it = 100

def beta_function(Alpha, Beta):
    #print Alpha, Beta
    return math.gamma(Alpha)*math.gamma(Beta)/math.gamma(Alpha + Beta)

def get_series(Alpha, Beta, x, total_it, it):
    if it < total_it:
        if it % 2 == 1:
            # odd 
            m = (it - 1)/2
            d_this = -(Alpha + m) * (Alpha + Beta + m) * x / ((Alpha + 2 * m ) * (Alpha + 2 * m + 1)) 
        else:
            # even
            m = it / 2
            d_this = m * (Beta - m ) * x / ((Alpha + 2*m - 1)*(Alpha + 2*m))
        return d_this / 1 + get_series(Alpha, Beta, x, total_it, it + 1) 
    else:
        return 0

def get_alpha_beta(x0, y0):
    def F(P):
        Alpha = P[0]
        Beta = P[1]
        #x0 = 0.45
        #y0 = 0.35
        #X = Alpha**2 + Beta**2 - x0
        #Y = Alpha**2 - Beta**2 - y0
        #x = Alpha / (Alpha + Beta)
        #series = 1 / (1 + get_series(Alpha, Beta, x, total_it, 1))
        #print series
        # approximate so that X and Y are zeros
        #X = (x**Alpha) * ((1 - x)**Beta) / (Alpha * beta_function(Alpha, Beta)) * series - x0
        #Y = (x**(Alpha + 1)) * ((1 - x)**Beta) / ((Alpha + 1) * beta_function((Alpha + 1), Beta)) * series - y0
        #X = ((Alpha/(Alpha+Beta))**Alpha) * ((1 - (Alpha/(Alpha+Beta)))**Beta) / (Alpha * beta_function(Alpha, Beta)) * (1 / (1 + get_series(Alpha, Beta, Alpha/(Alpha+Beta), total_it, 1))) - x0
        #Y = ((Alpha/(Alpha+Beta))**(Alpha + 1)) * ((1 - (Alpha/(Alpha+Beta)))**Beta) / ((Alpha + 1) * beta_function((Alpha + 1), Beta)) * (1 / (1 + get_series(Alpha, Beta, Alpha/(Alpha+Beta), total_it, 1))) - y0
        X = beta.cdf(float(Alpha)/(Alpha+Beta), Alpha, Beta) - x0
        Y = beta.cdf(float(Alpha)/(Alpha+Beta), Alpha+1, Beta) - y0
        return [X, Y]

    guess = [10, 5]
    max_n = 1000
    n = 1
    while n < max_n:
        try:
            sol = newton_krylov(F, guess, method = 'lgmres', verbose = 0, rdiff = 0.1, maxiter=50)
            break
        except NoConvergence as e:
            guess = np.random.rand(2) * 10 + 0.1
            print(guess)
        except ValueError as e:
            guess = np.random.rand(2) * 10 + 0.1
            print(guess)
        n = n + 1
    if n == max_n:
        print("Error: the chosen value is not within the computation range. Please re-choose. Suggested value: 0.3 <= x, y <= 0.5, and x >= y")
        return [0, 0]
    else:
        return sol

def get_beta(u, v):
    return get_alpha_beta(u, v)
    #Alpha = guess[0]
    #Beta = guess[1]
    #X = ((Alpha/(Alpha+Beta))**Alpha) * ((1 - (Alpha/(Alpha+Beta)))**Beta) / (Alpha * beta_function(Alpha, Beta)) * (1 / (1 + get_series(Alpha, Beta, Alpha/(Alpha+Beta), total_it, 1)))
    #Y = ((Alpha/(Alpha+Beta))**(Alpha + 1)) * ((1 - (Alpha/(Alpha+Beta)))**Beta) / ((Alpha + 1) * beta_function((Alpha + 1), Beta)) * (1 / (1 + get_series((Alpha+1), Beta, Alpha/(Alpha+Beta), total_it, 1)))
    #print X, Y
