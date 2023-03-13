# generate a reference array (of chrs) given CNs from an old reference array (of chrs)
# all amplifications are tandem
# return both alleles

from CN import CN
import os

# does not generate a fa file, but return the variable
def gen_ref_from_tree(ID, tree, ref):
    trace = [ID]
    visit = ID
    while visit != 0:
        visit = tree[visit].parentID
        trace.append(visit)
    # now reverse the trace so that the CNV and SNVs can be applied from root to leaf (leaf exclded)
    AB = []
    for i in range(len(trace) - 1):
        j = len(trace) - i - 1
        # gather all CNs together
        for ab in tree[trace[j]].aberrations:
            AB.append(ab)
    return gen_ref_wAberration(ref, AB)


def make_fa(ID, tree, ref, chr_name_array, fa_prefix):
    # make it generalizable when not to dump to a fa
    fa_f_prefix = fa_prefix + str(ID) + "_"
    # record all the nodes that this route visited (from leaf to root)
    trace = [ID]
    visit = ID
    while visit != 0:
        visit = tree[visit].parentID
        trace.append(visit)
    # now reverse the trace so that the CNV can be applied from root to leaf
    CN = []
    for i in range(len(trace)):
        j = len(trace) - i - 1
        # gather all CNs together
        for cn in tree[trace[j]].cn:
            CN.append(cn)
    new_ref = gen_ref(ref, CN)
    write_ref(new_ref, chr_name_array, fa_f_prefix)

def make_fa_wABs(ID, tree, ref, chr_name_array, fa_prefix):
    fa_f_prefix = fa_prefix + str(ID) + "_"
    # record all the nodes that this route visited (from leaf to root)
    trace = [ID]
    visit = ID
    while visit != 0:
        visit = tree[visit].parentID
        trace.append(visit)
    # now reverse the trace so that SNV and CNV can be applied from root to leaf
    AB = []
    for i in range(len(trace)):
        j = len(trace) - i - 1
        # gather all aberrations
        for k in range(len(tree[trace[j]].aberrations)):
            AB.append(tree[trace[j]].aberrations[k])
        #AB.append(tree[trace[j]].aberrations)
    new_ref = gen_ref_wAberration(ref, AB)
    write_ref(new_ref, chr_name_array, fa_f_prefix)


def getlen_ref(template):
    # get only the chromosome name and length, to save space when it's not necessary
    chr_name = []
    len_chr = []
    file = open(template,"r")
    # initialize string
    len_chr_ = 0
    line = file.readline().rstrip('\n')
    while(line != ""):
        if line[0] == '>':
            chr_name.append(line[1:])
            if len_chr_ != 0:
                len_chr.append(len_chr_)
            len_chr_ = 0
        else:
            len_chr_ = len_chr_ + len(line)
        line = file.readline().rstrip('\n')
    len_chr.append(len_chr_)
    file.close()
    return chr_name, len_chr


def init_ref(template):
    # read the template fasta file and save to ref array
    ref = []
    chr_name = []
    len_chr = []
    file = open(template,"r")
    line = "tmp"
    # initialize string
    str_ = ""
    line = file.readline().rstrip('\n')
    while(line != ""):
        if line[0] == '>':
            chr_name.append(line[1:].replace(" ", "_").replace(":", "_").replace("-", "_"))
            if str_ != "":
                ref.append(str_)
                len_chr.append(len(str_))
            str_ = ""
        else:
            str_ = str_ + line
        line = file.readline().rstrip('\n')
    ref.append(str_)
    len_chr.append(len(str_))
    file.close()
    ref_diploid = [ref, ref]
    return ref_diploid, chr_name, len_chr

def gen_ref(ref, CNs):
    # return this reference
    ret_ref = [row[:] for row in ref]
    for i in range(len(CNs)):
        ale = CNs[i].get_CN_Ale()
        pos1, pos2 = CNs[i].get_CN_position()
        chr_ = CNs[i].get_CN_chromosome()
        if CNs[i].CN_Del == 0:
            # amplification
            amp_num = CNs[i].get_CN_amp_num()
            #print amp_num
            #print ret_ref[ale][chr_][pos1:pos2]
            str_amp = amp_num * ret_ref[ale][chr_][pos1:pos2]
            # a bug found here. If str_amp = 0, the first string should be up to pos2, so that nothing is lost.
            ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos2] + str_amp + ret_ref[ale][chr_][pos2:]
        else:
            # deletion
            ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos1] + ret_ref[ale][chr_][pos2:]
    return ret_ref

# generate the reference with both SNVs and CNs
def gen_ref_wAberration(ref, aberrations):
    # return this reference
    ret_ref = [row[:] for row in ref]
    # in the order of the aberrations
    for i in range(len(aberrations)):
        print("Dealing with the " + str(i) + "th aberrations in gen_ref_wAberration")
        if aberrations[i].ab_type == "SNV":
            SNV = aberrations[i].SNV
            ale = SNV.ale
            chr_ = SNV.chr
            pos = SNV.pos
            pos1 = pos + 1
            new_nuc = SNV.new_nuc
            ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos] + new_nuc + ret_ref[ale][chr_][pos1:]

        elif aberrations[i].ab_type == "CNA": 
            CN = aberrations[i].CN
            ale = CN.get_CN_Ale()
            pos1, pos2 = CN.get_CN_position()
            chr_ = CN.get_CN_chromosome()
            if CN.CN_Del == 0:
                # amplification
                amp_num = CN.get_CN_amp_num()
                #print amp_num
                #print ret_ref[ale][chr_][pos1:pos2]
                str_amp = amp_num * ret_ref[ale][chr_][pos1:pos2]
                # a bug found here. If str_amp = 0, the first string should be up to pos2, so that nothing is lost.
                ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos2] + str_amp + ret_ref[ale][chr_][pos2:]
            else:
                # deletion
                ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos1] + ret_ref[ale][chr_][pos2:]
    return ret_ref


# read reference from a file to an array
def read_ref(ref):
    ref_a = []
    for i in [1, 2]:
        file_ = ref + str(i) + ".fa" 
        if os.path.isfile(file_):
            file = open(file_, "r")
            ref_ = []
            str_ = ""
            for line in file:
                line = line.strip()
                if line[0] == ">":
                    if str_ != "":
                        # if there is something in string
                        # add it to ref_, in the same order of the fasta file
                        ref_.append(str_)
                        str_ = ""
                else:
                    str_ = str_ + line
            if str_ != "":
                ref_.append(str_)
            file.close()
            ref_a.append(ref_)
        else:
            print(file_ + " does not exist and cannot be opened.")
    return ref_a
   
def write_ref(ref, chr_name, fasta):
    # write the reference to a fasta file
    line_len = 60
    for i in [1, 2]:
        file = open(fasta + str(i) + ".fa", "w")
        for j in range(len(ref[i-1])):
            file.write(">"+chr_name[j]+"\n")
            if len(ref[i-1][j]) < line_len:
                file.write(ref[i-1][j]+"\n")
            else:
                tmp_str = ref[i-1][j]
                a = 0
                while a + line_len < len(tmp_str):
                    file.write(tmp_str[a:a+line_len]+"\n")
                    a = a + line_len
                if len(tmp_str) > 0:
                    file.write(tmp_str[a:] + "\n")
        file.close()
