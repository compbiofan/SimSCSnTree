# Given Lorenz curve (approximated by Beta distribution given by Alpha and Beta), a recommended region size, coverage, a reference, acceptance threshold, sample reads from it so that the distribution is accorindg to the Lorenz curve. 
from scipy.stats import beta, norm
from estimate_beta_dist import get_beta
import numpy as np

# Lorenz curve
#Alpha = 0.2
#Beta = 0.5

# coverage
#coverage = 0.25

# reference fasta
#fasta = "test.fa"

# region size
#region_size = 50000

# acceptance threshold
#u = 0.5

#def random_sampling:
#    return np.random.rand(1)

def proposal_probability(a, b, sigma):
    # define Gaussian
    return norm.pdf(abs(b - a)/sigma)

def proposal(x_p):
    # sampling from the same Gaussian distribution
    new_x = np.random.normal(x_p)
    return new_x

# given the point that is furthest from the diagonal on a lorenz curve, estimate the Alpha and Beta for the underlying distribution (Beta)
def get_beta_dist(x0, y0):
    return get_beta(x0, y0)

def beta_probability(x, Alpha, Beta):
    return beta.pdf(x, Alpha, Beta)

def metropolis_hasting(prob_ratio, prop_ratio, u):
    if min(1, prob_ratio * prop_ratio) > u:
        return True
    else:
        return False

# given the previous window, decide how many reads (in terms of the parameter on Beta distribution) to sample in this window
def Markov_sampling(x_p, prob_x_p, Alpha, Beta, u, sigma):
    ini = False
    new_x = 0
    while ini is False:
        new_x = proposal(x_p)
        new_x_p = beta_probability(new_x, Alpha, Beta)
        prob_ratio = new_x_p / prob_x_p
        #prop_ratio = proposal_probability(x_p, new_x, sigma) / proposal_probability(new_x, x_p, sigma)
        # Gaussian, it is symmetric
        prop_ratio = 1
        ini = metropolis_hasting(prob_ratio, prop_ratio, u)
    return [new_x, new_x_p]

# get the # of reads for the next window
def get_next(x_p, prob_x_p, Alpha, Beta, u, sigma):
    [new_x, new_x_p] = Markov_sampling(x_p, prob_x_p, Alpha, Beta, u, sigma)
    return [new_x, new_x_p]

# get the mean read number based on given coverage and window size
# this number is corresponding to Alpha / (Alpha + Beta) = x0 
# this number / x0 * sampled_percentage is the read number at the new window
def get_mean(cov, window_size, readlen):
    return int(float(cov * window_size) / float(readlen))

# starting from here
x0 = float(raw_input("x axis furthest point from diagonal (0, 1):"))
y0 = float(raw_input("y axis furthest point from diagonal (0, 1), x >= y:"))
cov = float(raw_input("coverage (usually < 1):"))
readlen = float(raw_input("read length (35-250):"))
window_size = float(raw_input("window size (25000000):"))
# acceptance rate
u = 0.5
# sigma for Gaussian
sigma = 0.5

[Alpha, Beta] = get_beta_dist(x0, y0)
x0 = Alpha/ (Alpha + Beta)
# TODO, chr len according to actual situation
chr_len = 30000000
n = int(chr_len / window_size)
ini = x0
print ini
mean_read = get_mean(cov, window_size, readlen)
# output: a vector of read numbers for each window
read_numbers = []
f = open("out_x_ps", "w")
for i in range(n):
    if i == 0:
        x_p = ini
        prob_x_p = beta_probability(x_p, Alpha, Beta)
    else:
        [x_p, prob_x_p] = get_next(x_p, prob_x_p, Alpha, Beta, u, sigma)
    f.write(str(x_p) + "\n")
    # transfer to read number
    read_p = x_p / x0 * mean_read
    read_numbers.append(read_p)

print read_numbers
print [Alpha, Beta]
f.close()




