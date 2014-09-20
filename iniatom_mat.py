from numpy import *
import re

set_printoptions(threshold=nan)


class init_atom_coef(object):

    """a class to extract the coefficients in the original atomic orbitals"""

    def __init__(self, string, nbf):

        self.regEX   = "(\s+(-?\d+\.\d+)E(\+|-)(\d+))+(\s+)?"
        self.regEX_0 = "\s+(-?\d+\.\d+)E(\+|-)(\d+)\s+(-?\d+\.\d+)E(\+|-)(\d+)\s+(-?\d+\.\d+)E(\+|-)(\d+)(\s+)?"
        self.regEX_1 = "\s+(-?\d+\.\d+)E(\+|-)(\d+)(\s+)?"
        self.regEX_2 = "\s+(-?\d+\.\d+)E(\+|-)(\d+)\s+(-?\d+\.\d+)E(\+|-)(\d+)(\s+)?"
        self.mat = zeros([2*nbf + 4, nbf],float64)
        self.matrix_1 = zeros([nbf, nbf],float64)
        self.matrix_2 = zeros([nbf, nbf],float64)
        self.eigenvalue = zeros([1,nbf],float64)
        self.N = nbf%3
        self.NN = nbf/3 if self.N == 0 else nbf/3 + 1
        self.lines = 0
        self.total_line = 0
        self.nbf = nbf
        self.string = string


    def str_to_num(self, integer, symbol, decimal):
        
        if symbol == "+":
            return double(integer)*10**double(decimal)
        else:
            return double(integer)*10**(-double(decimal))


    def getter_spin_1(self):

        return self.matrix_1


    def getter_spin_2(self):

        return self.matrix_2
        
    def getter_energy_1(self):
        
        #print self.eigenvalue_1.shape
        return self.eigenvalue_1 * 27.211396132
        
    def getter_energy_2(self):
        
        return self.eigenvalue_2 * 27.211396132
        
    def getter_occupancy_1(self):
        
        return self.occupy_1
        
    def getter_occupancy_2(self):
        
        return self.occupy_2

    def extract_coef(self):

        f = open(self.string,'r')
        line = f.readline()
        while line != '':
            entry = re.match(self.regEX,line)
            if entry != None:
                self.total_line += 1    
            line = f.readline()

        f.seek(0,0)
        line = f.readline()
        ##print lines, total_line, f.readline()

        while line != '':
            entry = re.match(self.regEX,line)
            if entry != None:
                self.lines += 1
                if self.lines == self.total_line:break
                i = (self.lines-1)/self.NN
                j = (self.lines-1)*3%(3*self.NN)
                if self.lines%self.NN != 0 or self.N == 0:
                    entry = re.match(self.regEX_0,line)
                    self.mat[i,j    ] = self.str_to_num(entry.group(1),entry.group(2),entry.group(3))
                    self.mat[i,j + 1] = self.str_to_num(entry.group(4),entry.group(5),entry.group(6))
                    self.mat[i,j + 2] = self.str_to_num(entry.group(7),entry.group(8),entry.group(9))
                else:
                    if self.N == 1:
                        entry = re.match(self.regEX_1,line)
                        self.mat[i,j] = self.str_to_num(entry.group(1),entry.group(2),entry.group(3))
                    else:
                        entry = re.match(self.regEX_2,line)
                        self.mat[i,j    ] = self.str_to_num(entry.group(1),entry.group(2),entry.group(3))
                        self.mat[i,j + 1] = self.str_to_num(entry.group(4),entry.group(5),entry.group(6))
            line = f.readline()

        f.close()
        self.matrix_1 = self.mat[2:(self.total_line-1)/(2*self.NN),:].transpose()
        self.matrix_2 = self.mat[(self.total_line-1)/(2*self.NN)+2:(self.total_line-1)/self.NN,:].transpose()
        self.eigenvalue_1 = self.mat[1,:]
        self.eigenvalue_2 = self.mat[(self.total_line-1)/(2*self.NN)+1,:]
        self.occupy_1 = self.mat[0,:]
        self.occupy_2 = self.mat[(self.total_line-1)/(2*self.NN),:]
        #print self.lines, self.matrix_1.shape, self.matrix_2.shape, self.eigenvalue_2, self.occupy_2
        #for i in range(self.nbf-5):
            #print dot(self.matrix_1[:,i],self.matrix_1[:,i+1].transpose())
        #print dot(self.matrix_1.transpose(), self.matrix_1)
        
##
##
#
#P = init_atom_coef("Desktop\working\uo2f4-2.txt", 297)
#P.extract_coef()
#print P.getter_energy_1()
#print P.getter_energy_2()