import overlap_mat as ol_mat
import iniatom_mat as ia_mat
from numpy import *


class Lowdin(object):
    
    """a class to transform coefficients expressed from
       atomic orbital basis to Lowdin orbital basis"""

    def __init__(self, str1, str2, nbf):

        self.nbf = nbf
        p = ia_mat.init_atom_coef(str1, nbf)
        p.extract_coef()
        self.Coef_spin_1 = p.getter_spin_1()
        self.Coef_spin_2 = p.getter_spin_2()
        self.energy_1 = p.getter_energy_1()
        self.energy_2 = p.getter_energy_2()
        self.occupy_1 = p.getter_occupancy_1()
        self.occupy_2 = p.getter_occupancy_2()

        pp = ol_mat.overlap_mat(str2, nbf)
        pp.extract_coef()
        self.S = pp.getter_overlap()
        #print self.S.shape
        
        

        S_value, S_vector = linalg.eigh(self.S)
        self.dimension = S_value.shape[0]
        self.Coef = zeros([2*self.dimension, 2*self.dimension],float64)
        self.energy = zeros([2*self.dimension, 1],float64)
        self.occupy = zeros([2*self.dimension, 1],float64)
        self.Coef[ :self.dimension, :self.dimension] = self.Coef_spin_1
        self.Coef[self.dimension: , self.dimension:] = self.Coef_spin_2
        self.energy[ :self.dimension, 0] = self.energy_1
        self.energy[self.dimension: , 0] = self.energy_2
        self.occupy[ :self.dimension, 0] = self.occupy_1
        self.occupy[self.dimension: , 0] = self.occupy_2
        #print S_value

        
        S_eigen_invsqrt = zeros([self.dimension, self.dimension],float64)
        for i in range(self.dimension):
            S_eigen_invsqrt[i,i] = pow(S_value[i],-0.5)
        #print S_eigen_invsqrt
        
        self.X_sym = dot(dot(S_vector, S_eigen_invsqrt), S_vector.transpose())
        self.X_can = dot(S_vector, S_eigen_invsqrt)

        self.Lowdin_1 = dot(dot(self.X_sym.transpose(), self.S), self.Coef_spin_1)
        self.Lowdin_2 = dot(dot(self.X_sym.transpose(), self.S), self.Coef_spin_2)
        self.Lowdin = zeros([2*self.dimension, 2*self.dimension],float64)
        self.Lowdin[ :self.dimension, :self.dimension] = self.Lowdin_1
        self.Lowdin[self.dimension: , self.dimension:] = self.Lowdin_2

        
        self.Canonical_1 = dot(dot(self.X_can.transpose(), self.S), self.Coef_spin_1)
        self.Canonical_2 = dot(dot(self.X_can.transpose(), self.S), self.Coef_spin_2)
        self.Canonical = zeros([2*self.dimension, 2*self.dimension],float64)
        self.Canonical[ :self.dimension, :self.dimension] = self.Canonical_1
        self.Canonical[self.dimension: , self.dimension:] = self.Canonical_2
        
        #print self.Lowdin.shape, self.Canonical.shape
        self.X = zeros([2*self.dimension, 2*self.dimension], float64)
        self.X[ :self.dimension, :self.dimension] = self.X_can
        self.X[self.dimension: , self.dimension:] = self.X_can

        self.lap = zeros([2*self.dimension, 2*self.dimension], float64)
        self.lap[ :self.dimension, :self.dimension] = self.S
        self.lap[self.dimension: , self.dimension:] = self.S

        
    def getter_X(self):
        
        return self.X
        
        
    def getter_overlap(self):
        
        return self.lap
        

    def getter_Lowdin(self):

        return self.Lowdin

    def getter_Canonical(self):

        return self.Canonical
        #return self.Coef
        
    def getter_spectrum(self):
        
        return self.energy
        
    def getter_occupancy(self):
        
        return self.occupy

    def check(self):

        print "--------------------------------------------"
        print "Check the identity C.transpose * S * C = I"
        print "--------------------------------------------"
        check = dot(dot(self.Coef_spin_1.transpose(), self.S), self.Coef_spin_1)
        for i in range(self.dimension):
            for j in range(self.dimension):
                if abs(check[i, j]) >= 0.00001:
                    print i, j, check[i, j]

        print "--------------------------------------------"
        print "Check the identity X.transpose * S * X = I"
        print "--------------------------------------------"
        check = dot(dot(self.X_sym.transpose(), self.S), self.X_sym)
        for i in range(self.dimension):
            for j in range(self.dimension):
                if abs(check[i, j]) >= 0.00001:
                    print i, j, check[i, j]
        
        print "---------------------------------------------"
        print "Check the orthonormality of Lowdin orbitals"
        print "---------------------------------------------"
        for i in range(self.dimension-1):
            print dot(self.Canonical_1[:,i], self.Canonical_1[:,i].transpose())
            print dot(self.Lowdin_1[:,i], self.Lowdin_1[:,i].transpose())

        for i in range(self.dimension-1):
            print dot(self.Canonical_1[:,i], self.Canonical_1[:,i+1].transpose())
            print dot(self.Lowdin_1[:,i], self.Lowdin_1[:,i+1].transpose())



#L = Lowdin("Desktop\working\uo2f4-2.txt", "Desktop\working\uo2f4-overlap.out", 297)
#L.check()
#L.transform()
