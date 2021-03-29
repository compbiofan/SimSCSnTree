# definition of copy number class: only one copy number defined. 

class CN:
    def __init__(self, CN_Ale, CN_Del, CN_chromosome, CN_p1, CN_p2, CN_amp_num, corres):
        self.CN_Ale = CN_Ale
        self.CN_Del = CN_Del
        self.CN_chromosome = CN_chromosome
        self.CN_p1 = CN_p1
        self.CN_p2 = CN_p2
        self.CN_amp_num = CN_amp_num
        self.corres = corres
    def get_CN_Ale(self):
        return self.CN_Ale
    def get_CN_Del(self):
        return self.CN_Del
    def get_CN_chromosome(self):
        return self.CN_chromosome
    def get_CN_position(self):
        #print self.CN_p1, self.CN_p2
        return self.CN_p1, self.CN_p2
    def get_CN_amp_num(self):
        return self.CN_amp_num
