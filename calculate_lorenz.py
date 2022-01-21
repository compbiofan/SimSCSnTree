# given a file with one column composed of the number of reads in a ceratin region, give the Lorenz curve. 
#infile = "test.bin"
import matplotlib.pyplot as plt
from scipy.stats import beta
infile = "out_x_ps"
f = open(infile, "r")
data = []
for line in f:
    l = line.strip()
    data.append(float(l))

data.sort()
total = sum(data)
total_elements = 100
x_element = int(len(data) / total_elements)


# plot and compare the sampled one's Lorenz and the true one
fig, ax = plt.subplots()
# num is the x-axis, the bins
# tmp-sum is the y-axis
num = 0
tmp_sum = 0
# perc_data records the y axis for each percentage of the x axis chance, assuming there are totally 100 points
perc_data = []
for i in data:
    tmp_sum = tmp_sum + i
    if num % x_element == 0:
        # record
        perc_data.append(tmp_sum/total)
    num = num + 1

ax.plot(range(total_elements), perc_data, 'k--', label='Sampled')
x_ = range(1000)
Alpha = 7.7
Beta = 7.7
x = []
y = []
for i in x_:
    i = float(i)/1000
    x.append((beta.cdf(i, Alpha, Beta))*100)
    y.append(beta.cdf(i, Alpha + 1, Beta))

ax.plot(x, y, 'k:', label='Theoretic')
ax.plot(50, 0.4, 'ro', label='Input Point')
plt.title('Comparison of Sampled and Theoretic (sample size = 10k)')
legend = ax.legend(loc='upper left', fontsize='x-large')

plt.show()
