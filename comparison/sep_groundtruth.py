import sys
if len(sys.argv) <= 1:
    print("""
    For the ground truth bed file that were generated in simulator, separate it into different samples, each with a bed file, with the given prefix. The last column is the leaf ID. 
    Usage: python 
    """ + sys.argv[0] + """ [input] [output_prefix] """) 
    sys.exit(0)

input_f = sys.argv[1]
output_prefix = sys.argv[2]
# get all the leaf ids out, each leaf id contains all the lines corresponding to it
h={}
with open(input_f, "r") as f:
    for l in f:
        values = l.rstrip().split()
        key = values[len(values)-1]
        l_ = str(values[0])
        for j in range(len(values)-2):
            l_ = l_ + "\t" + str(values[j + 1])
        l_ = l_ + "\n"
        
        if key in h.keys():
            h[key].append(l_)
        else:
            h[key] = []
            h[key].append(l_)

# write each file
for i in h:
    f_ = output_prefix + str(i) + ".bed"
    f = open(f_, "w")
    for j in range(len(h[i])):
        f.write(h[i][j])
    f.close()

