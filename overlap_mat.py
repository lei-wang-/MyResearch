from numpy import *
import re

set_printoptions(threshold=nan)


class overlap_mat(object):

    """a class to extract the overlap matrix S"""

    def __init__(self, string, nbf):

        self.regEX = "1eov\s+(\d+)\s+\w+\s+\w+\s+(\d+)\s+\w+\s+\w+\s+(-?\d+\.\d+)(\s+)?"
        self.overlap = zeros([nbf,nbf],float64)
        self.nbf = nbf
        self.string = string

    def getter_overlap(self):

        return self.overlap


    def extract_coef(self):

        f = open(self.string,'r')
        line = f.readline()

        while line != '':
            entry = re.match(self.regEX,line)
            if entry != None:
                self.overlap[int(entry.group(1))-1, int(entry.group(2))-1] = entry.group(3)
                #self.overlap[int(entry.group(2))-1, int(entry.group(1))-1] = entry.group(3)
                #print entry.group(0), entry.group(1), entry.group(2), entry.group(3)
            line = f.readline()

        f.close()

#        for i in range(self.nbf):
#            for j in range(self.nbf):
#                if i == j:
#                    print i+1, j+1, self.overlap[i, j]
#
#P = overlap_mat("Desktop\working-space\uo2f4-overlap.out",297)
#P.extract_coef()
