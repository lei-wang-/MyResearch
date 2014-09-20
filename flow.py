from matplotlib import cm
from numpy import *
from pylab import *
import pickle
import copy
import random
import timeit

class Flowequation(object):

    """Class that reduces the size of Hamiltonian by
       flow equation using 4th order Runge-Kutta algorithm"""

    dt = 5e-20
    set_printoptions(threshold=nan)

    def __init__(self, endtime = 1.5e+01, n = 47, n00 = 14):


        f = open('./hamiltonian.pkl', 'r')
        self.H0 = pickle.load(f)
        f.close()
        f = open('./spectrum.pkl','r')
        self.spectrumm = pickle.load(f)
        s = self.spectrumm.shape[0]
        f.close()
        #f = open('./self.full_r.pkl','r')
        #self.full_r = pickle.load(f)
        #f.close()
        self.spect = zeros([s, s], float64)
        fill_diagonal(self.spect, self.spectrumm)
        self.spectrum, a = linalg.eigh(self.spect)
        
        self.M1 = zeros([n, n], float64)
        self.M2 = zeros([n, n], float64)
        self.M1 = self.H0[:n , :n]
        self.M2 = self.H0[n: , n:]

        self.alpha = 0.000
        self.endtime = endtime
        self.n00 = n00
        self.n = n
        self.f = 7
        self.nf = n - 7
        self.nn = n - n00
        self.P = zeros([n, n], float64)
        dia_elements = [1.0 for i in range(self.nf)] + [0 for i in range(self.f)]
        #dia_elements = [1.0 for i in range(self.n)]
        fill_diagonal(self.P, dia_elements)        
        self.t = 0.0
        
        self.e = self.spectrum   #full 
        self.spectrum_size = self.e.shape[0]
        
        self.precision = 1.0e-9
        self.magnitude = []


    def split_M1(self):
        
        
        self.H = copy.deepcopy(self.M1)
        
        
    def split_M2(self):
        
        
        self.H = copy.deepcopy(self.M2)


    def derivative(self, H):


        Hd = zeros([self.n, self.n],float64)
        #H = linalg.inv(H)
        fill_diagonal(Hd, [(H[i, i]) for i in range(self.n)])
        
        
        #if (self.t < 1.0e-1):    
        #    fill_diagonal(Hd, [(H[i, i]) for i in range(self.n)])
        #    #if(self.t > 1.0e-1 - self.dt):
        #        #self.dt = 1.0e-20
        #else:
        #    fill_diagonal(Hd, [(H[i, i]) for i in range(self.n)])
        
        
        g = dot(Hd,H) - dot(H,Hd)
        G = dot(dot(self.P,g),self.P)
        rhs = dot(G,H) - dot(H,G)
        return rhs


    def step(self, dt, H):


        k1 = dt * self.derivative(H)
        k2 = dt * self.derivative(H + 0.5*k1)
        k3 = dt * self.derivative(H + 0.5*k2)
        k4 = dt * self.derivative(H + k3)
        H += 1.0/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)
        return H


    def algo(self):


        while(self.t < self.endtime):                                                                                              
            
            H_temp = copy.deepcopy(self.H)
            H_temp = self.step(self.dt, H_temp)

            self.H = self.step(0.5*self.dt, self.H)
            self.H = self.step(0.5*self.dt, self.H)
            
            error = max(abs(ravel(H_temp - self.H)))
            
            #magnitudes = mean(sum(abs(ravel(self.H[:self.nf, :self.nf]))) - sum(abs(self.H[:self.nf, :self.nf].diagonal())))
            #print "magnidue ", magnitudes
            h = zeros([self.nf, self.nf])
            fill_diagonal(h, [self.H[i, i] for i in range(self.nf)])
            magnitudes = mean(ravel(abs(self.H[:self.nf, :self.nf] - h)))
            self.magnitude.append(magnitudes)
            if magnitudes < 0.0003:
                break
            self.t += self.dt
            print "error ", error, "at time", self.t, "by time step ", self.dt
            print "magnitude ", magnitudes
            if(error < self.precision):
                self.dt = self.dt*pow(self.precision/error, 0.20)
            else:
                self.dt = self.dt*pow(self.precision/error, 0.25)
            
        self.t = 0.0
        self.dt = 5e-20
        
        
    def rearrange1(self):


        thresh = [0.0 for k in range(self.nf)]
        lists = []
        
        for p in range(self.nf):
            for q in range(self.n):
                if p != q and q >= self.nf:
                    thresh[p] += abs(self.H[q , p]) ** 2
                if p != q and q < self.nf:
                    thresh[p] += self.alpha * abs(self.H[q , p]) ** 2
                
        for p in range(self.nf):
            if (thresh[p] < 1e-04):
                thresh[p] = 0.0   
                
        #print thresh
                
        for p in range(self.nf):
            k = 0
            for q in range(self.nf):
                if thresh[p] < thresh[q]:
                    k = k + 1
            lists.append(k)
            
        #print lists
        
        for p in range(self.nf):
            k = 0
            for q in arange(p, self.nf, 1):
                if (lists[p] == lists[q]) and p != q:
                    k += 1
            lists[p] += k
                    
        #print lists
        
        index = array(lists)
        extra = self.nf - self.nn
        temp = 0
        for p in arange(self.nn, self.nf, 1):
            if index[p] >= extra:
                for q in arange(temp, self.nn, 1):
                    if index[q] < extra:
                        #print "p", p, "q", q, "change"
                        self.exchange(self.H[:,p], self.H[:,q])
                        self.exchange(self.H[p,:], self.H[q,:])
                        temp = q + 1
                        break

            
    def rearrange2(self, e):
        
        
        thresh = [0.0 for k in range(self.nf)]
        lists = []
        
        for p in range(self.nf):
            thresh[p] = abs(self.H[p, p] - e)
                
        for p in range(self.nf):
            if (thresh[p] < 1e-04):
                thresh[p] = 0.0   
                
        #print thresh
                
        for p in range(self.nf):
            k = 0
            for q in range(self.nf):
                if thresh[p] > thresh[q]:
                    k = k + 1
            lists.append(k)
            
        #print lists
        
        for p in range(self.nf):
            k = 0
            for q in arange(p, self.nf, 1):
                if (lists[p] == lists[q]) and p != q:
                    k += 1
            lists[p] += k
                    
        #print lists
        
        index = array(lists)
        extra = self.nf - self.nn
        temp = 0
        for p in arange(self.nn, self.nf, 1):
            if index[p] >= extra:
                for q in arange(temp, self.nn, 1):
                    if index[q] < extra:
                        #print "p", p, "q", q, "change"
                        self.exchange(self.H[:,p], self.H[:,q])
                        self.exchange(self.H[p,:], self.H[q,:])
                        temp = q + 1
                        break

            
    def exchange(self, a1, a2):
        
        
        for i in range(len(a1)):
            a1[i], a2[i] = a2[i], a1[i]
                                
            
    def reduce(self):
        
        
        diagonal_element = [self.H[k,k] for k in range(self.nf)]
        print "diagonal elements before truncation."
        print sort(array(diagonal_element))
        self.new_size = self.n        
        
        #up_bound = 0.0
        #low_bound = 2.20
        ##3.52 for uo2f4-2; -3.270 for uo2f4-1; -6.19 for uf5; 0.65 for uf5-1
        #fermi = 3.52
        #offset = 0
        #kk = 0
        #print "truncation window from: ", fermi-low_bound, "to:", fermi+up_bound
        self.rearrange2(e=3.52)
        self.H_save = copy.deepcopy(self.H)
        
        #for k in range(self.nf):           
        #    if (fermi-low_bound) > diagonal_element[k] or (fermi+up_bound) < diagonal_element[k]:
        #        for i in range(self.new_size):
        #            for j in range(self.new_size):
        #                kk = k - offset
        #                if i < kk and j >= kk and j < self.new_size - 1:
        #                    self.H[i, j] = self.H[i, j+1]
        #                elif j < kk and i >= kk and i < self.new_size - 1:
        #                    self.H[i, j] = self.H[i+1, j]
        #                elif i >= kk and i < self.new_size - 1 and j >= kk and j < self.new_size - 1:
        #                    self.H[i, j] = self.H[i+1, j+1]
        #                else:
        #                    self.H[i, j] = self.H[i, j]
        #        self.H = self.H[:self.new_size-1 , :self.new_size-1]
        #        self.new_size = self.H.shape[0]
        #        offset += 1
                
        #print self.H.shape , offset
        
        self.H = self.H[self.nn: , self.nn:]
        a, b = linalg.eigh(self.H)
        #print "diagonal elements after truncation."
        #print [self.H[i,i] for i in range(len(b))]

            
    def spec(self):
        
        
        e2, v2 = linalg.eigh(self.H)
        
        print "********eigenvalues**********"
        print e2
        print "*****************************"

        dicts = {}
        for i in range(len(e2)):
            dicts[str(e2[i]+random.random()/10**6)] = v2[:, i]
        
        #print dicts
        print "******************"
        return dicts

        
    def save(self, s):
        
        
        f = open('./Ham_M1.pkl','w')
        pickle.dump(self.M1, f)
        f.close()
        f = open('./Ham_M2.pkl','w')
        pickle.dump(self.M2, f)
        f.close()
        f = open('./Ham_flow_'+s+'.pkl','w')
	pickle.dump(self.H_save, f)
	f.close()
	f = open('./Ham_reduce_'+s+'.pkl','w')
	pickle.dump(self.H, f)
	f.close()
    
    def ssave(self):
        
        f1 = open('./Ham_reduce_M1.pkl','r')
        m1 = pickle.load(f1)
        f1.close()
        f2 = open('./Ham_reduce_M2.pkl','r')
        m2 = pickle.load(f2)
        f2.close()
        m = zeros([len(m1)+len(m2), len(m1)+len(m2)], float64)
        m[:len(m1), :len(m1)] = m1
        m[len(m1):, len(m1):] = m2
        f = open('./Ham_reduce.pkl','w')
        pickle.dump(m, f)
        f.close()

        
    def hex_to_rgb(self, value):
        
        
        value = value.lstrip('#')
        lv = len(value)
        return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3))
        
        
    def rgb_to_hex(self, rgb):
        
        
        return '#' + '%02x%02x%02x' % rgb
                    
        
    def plot(self, f1, f2):
        
        
        elist=[]  #flow
        ew=[]
        left_1=[]
        color_1 = []
        
        eelist=[]  #full
        eew=[]
        left_2=[]
        color_2 = []
        
        f1.update(f2)
        print type(f1)
        e = sort(f1.keys())
        print '++++++++++++++++++++'
        print e
        print '+++++++++++++++++++'

        for i in range(len(e)):
            elist.append(float(e[i]))   #flow
            ew.append(0.003)
            
            if i == 0:
                left_1.append(0.032)
            elif abs(float(e[i]) - float(e[i-1])) <= 0.0050:  #0.038
                left_1.append(left_1[i-1] + 0.004)
            else:
                left_1.append(0.032)
            
            r =  sum(f1[e[i]][-7:] ** 2)/sum(ravel(f1.values()[-7:]) ** 2) * 255 * 5
            #print r
            rgb = self.rgb_to_hex((r,0,0))
            color_1.append(rgb)
            #print sum(f1[e[i]][-7:] ** 2)/sum(f1[e[i]] ** 2)
            
        color_2 = self.full_r    
        for i in range(self.spectrum_size):
            eelist.append(self.e[i])      #full
            eew.append(0.003)
            
            if i == 0:
                left_2.append(0.004)
            elif abs(self.e[i] - self.e[i-1]) <= 0.0050:
                left_2.append(left_2[i-1] + 0.004)
            else:
                left_2.append(0.004)
                                 
           
        figure(1)
        title("Flow-equation Spectrum")
        ylabel("Energy(eV)")
        xlabel("")
        ax = gca()
        ax.set_xlim(0, 0.05)
        # (0.20, 9.50) for uo2f4-2; (-5.20, 1.50) for uo2f4-1;
        # (-10.8, -4.0) for uf5; (-4.8, 2.2) for uf5-1
        ax.set_ylim(0.20, 9.50)
        tick_lcs = array([0.013, 0.037])
        tick_lbs = array(['(a)','(b)'])
        xticks(tick_lcs, tick_lbs)
        draw()
        #3.52 for uo2f4-2; -3.270 for uo2f4-1; -6.19 for uf5; 0.65 for uf5-1
        axhline(y=3.52, xmin=0.20, xmax=0.80, ls='--')
        barh(elist, ew, height=0.001, left = left_1, color = color_1, edgecolor = color_1)        #flow 0.032
        barh(eelist, eew, height=0.001, left = left_2, color = color_2, edgecolor = color_2)      #full  0.008
        show()


    def plot2(self):
        
        
        figure(2)
        plot(arange(len(self.magnitude)), self.magnitude)
        title("Matrix Element Evolution")
        xlabel("stepsize")
        ylabel("Largest Non-diagonal Element of Non-f Block Matrix")
        show()        
        
        
        

f = Flowequation()
start = timeit.default_timer()
f.split_M1()
f.algo()
end = timeit.default_timer()
f.reduce()
f1 = f.spec()
f.save("M1")
#f.plot2()


f.split_M2()
f.algo()
f.reduce()
f2 = f.spec()
f.save("M2")
#f.plot(f1, f2)
f.ssave()
print "time consumed:", end-start
