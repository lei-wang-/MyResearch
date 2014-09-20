import pickle
from numpy import *
import iniatom_mat as ia_mat
import overlap_mat as ol_mat

class matrix(object):
    
    """a class for mapping 49*2 gaussian f functions onto 7*2 atomic f functions"""

    def __init__(self, s1, s2,l):
        
        f = ia_mat.init_atom_coef(s1, l)
        f.extract_coef()
        Coef_spin_1 = f.getter_spin_1()
        Coef_spin_2 = f.getter_spin_2()
        #Coef_spin_2 = f.getter_spin_1()


        ff = ol_mat.overlap_mat(s2, l)
        ff.extract_coef()
        S = ff.getter_overlap()
        S_value, S_vector = linalg.eigh(S)
        self.dimension = S_value.shape[0]
        S_eigen_invsqrt = zeros([self.dimension, self.dimension],float64)
        for i in range(self.dimension):
            S_eigen_invsqrt[i,i] = pow(S_value[i],-0.5)
        X_can = dot(S_vector, S_eigen_invsqrt)
        self.Canonical_1 = dot(dot(X_can.transpose(), S), Coef_spin_1)
        self.Canonical_2 = dot(dot(X_can.transpose(), S), Coef_spin_2)
        
        self.X = zeros([2*self.dimension, 2*self.dimension], float64)
        self.X[ :self.dimension, :self.dimension] = X_can
        self.X[self.dimension: , self.dimension:] = X_can
        
        
        
    def f_orbital(self):
        
        "F orbitals order: -3, -2, -1, 0, 1, 2, 3"
        self.F_alpha = zeros([7, 159], float64)
        self.F_beta  = zeros([7, 159], float64)
        self.F_alpha[0, :] = self.Canonical_1[:, 44-1].transpose()
        self.F_alpha[1, :] = self.Canonical_1[:, 45-1].transpose()
        self.F_alpha[2, :] = self.Canonical_1[:, 46-1].transpose()
        self.F_alpha[3, :] = self.Canonical_1[:, 47-1].transpose()
        self.F_alpha[4, :] = self.Canonical_1[:, 48-1].transpose()
        self.F_alpha[5, :] = self.Canonical_1[:, 49-1].transpose()
        self.F_alpha[6, :] = self.Canonical_1[:, 50-1].transpose()
 
        #53-1,...,59-1 for U+0
        #44-1,...,50-1 for U+6
        
        self.F_beta[0, :] = self.Canonical_2[:, 44-1].transpose()
        self.F_beta[1, :] = self.Canonical_2[:, 45-1].transpose()
        self.F_beta[2, :] = self.Canonical_2[:, 46-1].transpose()
        self.F_beta[3, :] = self.Canonical_2[:, 47-1].transpose()
        self.F_beta[4, :] = self.Canonical_2[:, 48-1].transpose()
        self.F_beta[5, :] = self.Canonical_2[:, 49-1].transpose()
        self.F_beta[6, :] = self.Canonical_2[:, 50-1].transpose()       
               
                             
        

    
    def check(self):
        
        for i in arange(43,50,1):
            print dot(self.Canonical_1[110:159,i] , self.Canonical_1[110:159,i].transpose())
            
        for i in arange(43,50,1):
            print dot(self.Canonical_1[:,i], self.Canonical_1[:,i].transpose())
    
    
        
    def save(self):
        
        f1 = open('./F_alpha.pkl', 'w')
        f2 = open('./F_beta.pkl', 'w')
        f3 = open('./F_lap.pkl', 'w')
        pickle.dump(self.F_alpha, f1)
        pickle.dump(self.F_beta, f2)
        pickle.dump(self.X, f3)
        #pickle.dump(self.Canonical_1.transpose(), f1)
        #pickle.dump(self.Canonical_2.transpose(), f2)
        f1.close()
        f2.close()
        f3.close()
        
        
    def load(self):
    
        f1 = open('./F_alpha.pkl', 'r')
        f2 = open('./F_beta.pkl', 'r')
        t1 = pickle.load(f1)
        t2 = pickle.load(f2)
        print t1.shape, t2.shape



h = matrix("./U+6.txt", "./u-overlap.out", 159)
h.f_orbital()
h.save()
h.load()
#h.check()
