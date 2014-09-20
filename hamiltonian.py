from numpy import *
import lowdin
import iniatom_mat as ia_mat
import pickle
from pylab import *
import random

class Hamiltonian(object):
    
    """a class to project the original Hamiltonian onto the reduced basis set"""
    
    def __init__(self, str1, str2, nbf, nrbf):
        
        
        L = lowdin.Lowdin(str1, str2, nbf)
        self.nrbf = nrbf * 2
        self.spectrum = L.getter_spectrum()
        self.matrix = L.getter_Canonical()
        self.occupy = L.getter_occupancy()
        self.size = self.spectrum.shape[0]
        self.hamiltonian = zeros([self.size, self.size],float64)
        fill_diagonal(self.hamiltonian, [self.spectrum[i] for i in range(self.size)])
        self.electron_1 = 0
        self.electron_2 = 0
        for i in range(self.size/2):
            self.electron_1 += int(self.occupy[i])
            self.electron_2 += int(self.occupy[self.size/2 + i])
        print self.electron_1, self.electron_2, self.size/2, type(self.matrix), self.size
        self.X_molecule = L.getter_X()
        self.lap = L.getter_overlap()
        #print self.spectrum

        f1 = open('./F_alpha.pkl', 'r')
        f2 = open('./F_beta.pkl', 'r')
        f3 = open('./F_lap.pkl', 'r')
        self.f_alpha = pickle.load(f1)
        self.f_beta  = pickle.load(f2)
        self.X_atom = pickle.load(f3)
        f1.close()
        f2.close()
        f3.close()
        self.size_1 = self.X_atom.shape[0]
        self.overlap = zeros([self.size_1 , self.size], float64)
        self.overlap[:self.size_1/2 , :] = self.lap[:self.size_1/2 , :]
        self.overlap[self.size_1/2: , :] = self.lap[self.size/2 : self.size/2+self.size_1/2 , :]
        self.temp = dot(dot(self.X_atom.transpose(), self.overlap), self.X_molecule)
        
        
    def sanity_check(self):
        
        
        New_hamiltonian = dot(dot(self.matrix, self.hamiltonian), self.matrix.transpose())
        New_spectrum, New_composition = linalg.eigh(New_hamiltonian)
        print New_spectrum
        
        
    def Trans_hamiltonian(self):
        
        
        self.redbasis = zeros([self.nrbf, self.size], float64)
        self.basis = zeros([14, self.size_1],float64)
        self.basis[:7 , :self.size_1/2] = self.f_alpha
        self.basis[7: , self.size_1/2:] = self.f_beta

        self.temp_1 = dot(self.basis, self.temp)
        self.redbasis[self.nrbf-14: , :] = self.temp_1
        
        self.hybridization = []
        hybridize = dot(self.temp_1, self.matrix)
        for i in range(len(self.matrix)):
            self.hybridization.append(sum(hybridize[:, i] ** 2))
        index = argsort(array(self.hybridization))
        print "Index of 14 KS orbitals with highest 5f hybridization"
        print index[-14:]
        print "-----------------------------------------------------"
        
        self.dicts = {}
        for i in range(len(self.hybridization)):
            self.dicts[ravel(self.spectrum)[i] + random.random()/10**6] = self.hybridization[i]
        
        div = (self.nrbf-14)/4.0
        
        #index[-14:]; index[120:134]
        if div > self.electron_1:
            k = 0
            kk = 0
            for i in range(int(2*div)):
                if i+k in index[120:134]:
                    print i+k
                    k += 1
                    kk = k 
                    for j in range(7):
                        if i+kk+j in index[120:134]:
                            print i+k
                            k += 1    
                        else:
                            break            
                    self.redbasis[i, :] = self.matrix[:, i+k].transpose()
                else:
                    self.redbasis[i, :] = self.matrix[:, i+k].transpose()     
        else:
            k = 0
            kk = 0
            for i in range(int(2*div)):
                if self.electron_1-int(div)+i+k in index[120:134]:
                    print self.electron_1-int(div)+i+k
                    k += 1
                    kk = k
                    for j in range(7):
                        if self.electron_1-int(div)+i+kk+j in index[120:134]:
                            print self.electron_1-int(div)+i+k
                            k += 1
                        else:
                            break
                    self.redbasis[i, :] = self.matrix[:, self.electron_1-int(div)+i+k].transpose()
                else:
                    self.redbasis[i, :] = self.matrix[:, self.electron_1-int(div)+i+k].transpose()
        
        
        
        if div > self.electron_2:
            k = 0
            kk = 0
            for i in range(int(2*div)):
                if self.size/2+i+k in index[120:134]:
                    print self.size/2+i+k
                    k += 1
                    kk = k 
                    for j in range(7):
                        if self.size/2+i+kk+j in index[120:134]:
                            print self.size/2+i+k
                            k += 1    
                        else:
                            break            
                    self.redbasis[2*div+i, :] = self.matrix[:, self.size/2+i+k].transpose()
                else:
                    self.redbasis[2*div+i, :] = self.matrix[:, self.size/2+i+k].transpose()     
        else:
            k = 0
            kk = 0
            for i in range(int(2*div)):
                if self.size/2+self.electron_2-int(div)+i+k in index[120:134]:
                    print self.size/2+self.electron_2-int(div)+i+k
                    k += 1
                    kk = k
                    for j in range(7):
                        if self.size/2+self.electron_2-int(div)+i+kk+j in index[120:134]:
                            print self.size/2+self.electron_2-int(div)+i+k
                            k += 1
                        else:
                            break
                    self.redbasis[2*div+i, :] = self.matrix[:, self.size/2+self.electron_2-int(div)+i+k].transpose()
                else:
                    self.redbasis[2*div+i, :] = self.matrix[:, self.size/2+self.electron_2-int(div)+i+k].transpose()
                       

        self.basisup = zeros([self.nrbf/2, self.size], float64)
        self.basisdown = zeros([self.nrbf/2, self.size], float64)
        self.basisup[ :2*div , :] = self.redbasis[ :2*div , :]
        self.basisup[2*div:self.nrbf/2 , :] = self.redbasis[4*div:4*div+7 , :]
        self.basisdown[ :2*div , :] = self.redbasis[ 2*div:4*div , :]
        self.basisdown[2*div:self.nrbf/2 , :] = self.redbasis[4*div+7:self.nrbf , :]
        
        #self.orthbasisup = self.gram_schmidt_1(self.nrbf/2, 7, self.basisup)
        #self.orthbasisdown = self.gram_schmidt_1(self.nrbf/2, 7, self.basisdown)
        
        self.orthbasisup = self.gram_schmidt_2(self.nrbf/2, self.basisup)
        self.orthbasisdown = self.gram_schmidt_2(self.nrbf/2, self.basisdown)
        
        self.orthbasis = zeros([self.nrbf, self.size], float64)
        self.orthbasis[:self.nrbf/2 , :] = self.orthbasisup
        self.orthbasis[self.nrbf/2: , :] = self.orthbasisdown
        
        temp_hamiltonian = dot(dot(self.matrix, self.hamiltonian), self.matrix.transpose())
        self.New_hamiltonian = dot(dot(self.orthbasis, temp_hamiltonian), self.orthbasis.transpose())
        self.New_spectrum, self.New_composition = linalg.eigh(self.New_hamiltonian)
        print shape(self.New_hamiltonian)
        
        
    def getter_hamiltonian(self):
        
        
        return self.New_hamiltonian
        
    
    def getter_Uo(self):
        
        
        return dot(self.orthbasis, self.matrix)
        
        
    def getter_h(self):
        
        
        return self.hamiltonian
      
          
    def project(self, u, v):


        return dot(u, v.transpose())/dot(u, u.transpose()) * u
        
        
    def gram_schmidt_1(self, m, n, basis):
        
        
        orthbasis = zeros([basis.shape[0], basis.shape[1]], float64)
        orthbasis[-7:, :] = basis[-7:, :]
        for k in range(m - n):
            summ = mat(zeros((1, self.size),float64))
            for j in range(n + k):
                summ += self.project(orthbasis[m-1-j,:], basis[m-n-1-k,:])
            orthbasis[m-n-1-k,:] = basis[m-n-1-k,:] - summ
            summm = 0.0
            for i in range(self.size):
                summm += pow(orthbasis[m-n-1-k, i],2)
            orthbasis[m-n-1-k,:] = 1/pow(summm,0.5) * orthbasis[m-n-1-k,:]
            
        return orthbasis
        

    def gram_schmidt_2(self, k, basis):
        
        
        for i in arange(k-1, -1, -1):
            basis[i, :] = basis[i, :]/sqrt(dot(basis[i, :], basis[i, :]))
            for j in arange(i-1, -1, -1):
                basis[j, :] = basis[j, :] - self.project(basis[i, :], basis[j, :])
        
        return basis

        
    def save(self):
        
        
        f = open('./hamiltonian.pkl','w')
        pickle.dump(self.New_hamiltonian, f)
        f.close()
        f = open('./spectrum.pkl','w')
        pickle.dump(self.spectrum, f)
        f.close()
        f = open('./statesmatrix.pkl','w')
        pickle.dump(self.matrix, f)
        f.close()
        f = open('./Uo.pkl','w')
        pickle.dump(matrix(self.getter_Uo()), f)
        f.close()
        f = open('./Hdft.pkl','w')
        pickle.dump(matrix(self.getter_h()), f)
        f.close()
        f = open('./Nks.pkl','w')
        Nks = multiply(resize(self.occupy, self.getter_Uo().shape) , multiply(self.getter_Uo(), self.getter_Uo().conjugate())).sum()
        pickle.dump(Nks, f)
        print self.occupy.shape, self.getter_Uo().shape
        f.close()
        f = open('./Occ.pkl','w')
        pickle.dump(self.occupy.transpose(), f)
        f.close()
        #f = open('./self.full_r.pkl','w')
        #pickle.dump(self.full_r, f)
        #f.close()
        
        
    def load(self):
        
        
        f = open('./Uo.pkl','r')
        print pickle.load(f).shape
        f.close()
        f = open('./Hdft.pkl','r')
        print pickle.load(f).shape
        f.close()
        f = open('./Nks.pkl','r')
        print pickle.load(f)
        f.close()
        
        
    def check(self):
        
        
        print "--------------------------------------------"
        print "Check the orthonormality of the reduced basis set"
        print "--------------------------------------------"
        for i in range(self.nrbf):
            print dot(self.redbasis[i,:], self.redbasis[i,:].transpose())
            
        for i in range(self.nrbf-3):
            print dot(self.redbasis[i,:], self.redbasis[i+2,:].transpose())
            
        print "  ",dot(self.orth_redbasis[30,:], self.orth_redbasis[66,:].transpose())
        print "  ",dot(self.redbasis[30,:], self.redbasis[66,:].transpose())


    def rgb_to_hex(self, rgb):
        
        
        return '#' + '%02x%02x%02x' % rgb


    def plot(self):
        
        
        elist=[]  #one-particle spectrum
        ewid=[]
        eleft=[]
        ecolor=[]
        
        eelist=[]  #full spectrum
        eewid=[]
        eeleft=[]
        eecolor=[]
        self.full_r=[]

        self.spectrumm = sort(self.dicts.keys()) 
        
        
        for i in range(len(self.New_spectrum)):  #one-particle spectrum
            elist.append(self.New_spectrum[i])
            ewid.append(0.0032) #0.002.....0.0032
            
            if i == 0:
                eleft.append(0.032)
            elif abs(float(self.New_spectrum[i]) - float(self.New_spectrum[i-1])) <= 0.0050:  #0.038
                eleft.append(eleft[i-1] + 0.004) #0.003.....0.004
            else:
                eleft.append(0.032)
            
            r1 = sum(self.New_composition[-7:, i] ** 2) + sum(self.New_composition[-7+self.nrbf/2: self.nrbf/2, i] ** 2) 
            #r2 = sum(self.New_composition[:, i] ** 2)
            r2 = sum(ravel(self.New_composition[-7:, :]) ** 2) + sum(ravel(self.New_composition[-7+self.nrbf/2: self.nrbf/2, :]) ** 2)
            r = r1/r2 * 255 * 13
            rgb = self.rgb_to_hex((r,0,0))
            ecolor.append(rgb)
            

        for i in range(len(self.spectrumm)):  #full spectrum
            eelist.append(self.spectrumm[i])
            eewid.append(0.0032) #0.002.....0.0032
            
            if i == 0:
                eeleft.append(0.004)
            elif abs(self.spectrumm[i] - self.spectrumm[i-1] <= 0.0050):
                eeleft.append(eeleft[i-1] + 0.004) #0.003.....0.004
            else:
                eeleft.append(0.004) 
            
            r = self.dicts[self.spectrumm[i]]/sum(self.dicts.values()) * 255 * 13
            rgb = self.rgb_to_hex((r,0,0))
            self.full_r.append(rgb)
            eecolor.append(rgb)
            
            
        figure(1)
        title("Energy Spectrum")
        xlabel("")
        ylabel("Energy(eV)")
        ax1 = gca()
        ax1.set_xlim(0, 0.05)
        # (0.20, 9.50) for uo2f4-2; (-5.20, 1.50) for uo2f4-1;
        # (-10.8, -4.0) for uf5; (-1.0, 5.5) for uf5-1
        # (-13.0, 1.0) for uf6; (-8.0, 4.0) for uf6-1
        # (-2.0, 7.0) for uo2cl4-2; (-5.8, 3.0) for uo2cl4-1
        ax1.set_ylim(0.20, 9.50)
        tick_lcs = array([0.012, 0.04])
        tick_lbs = array(['(a)','(b)'])
        xticks(tick_lcs, tick_lbs)
        draw()
        #3.52 for uo2f4-2; -3.270 for uo2f4-1; -6.19 for uf5; 0.65 for uf5-1
        #-8.75 for uf6; -0.53 for uf6-1; 1.93 for uo2cl4-2; -3.47 for uo2cl4-1
        axhline(y=3.52, xmin=0.20, xmax=0.80, ls='--')
        barh(eelist, eewid, height=0.001, left=eeleft, color=eecolor, edgecolor=eecolor)   #full spectrum
        barh(elist, ewid, height=0.001, left=eleft, color=ecolor, edgecolor=ecolor)     #one-particle spectrum
        show()
        
        
        figure(2)
        ax2 = gca()
        ax2.set_xlim(0, self.nrbf/2)
        ax2.set_ylim(0, 0.08)
        p = self.hybridization / sum(self.hybridization)
        #print sum(p)
        plot(p)
        xlabel("Index of Kohn-Sham orbital state")
        ylabel("5f orbital content")
        axvline(x=self.electron_1 - 1, color='red')
        axvline(x=self.nrbf/2 + self.electron_2 - 1, color='red')
        show()
        
        
        figure(3)
        ax3 = gca()
        ax3.set_xlim(-20, 20)
        ax3.set_ylim(0, 0.08)
        p = self.dicts.values() / sum(self.dicts.values())
        bar(self.dicts.keys(), p, width=0.01, color='blue', edgecolor='blue')
        axvline(x=3.52, color='red')
        xlabel("Energy of Kohn-Sham orbital state(eV)")
        ylabel("5f orbital content")
        show()
        
        

##(297, 47) for uo2f4-1(-2); (274, 47) for uf5(-1)
##(313, 47) for uo2cl4-1(-2); (297, 47) for uf6(-1)
h = Hamiltonian("./uo2f4-2.txt", "./uo2f4-2-overlap.out", 297, 47)
h.Trans_hamiltonian()
#h.plot()
h.save()
h.load()
