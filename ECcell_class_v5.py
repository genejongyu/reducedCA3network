#############################################################################################
# ECcell_class.py
# by Gene Yu
# 01/20/2012
#
# This function creates the entorhinal cortex (EC) cell class. 
#############################################################################################

from neuron import h
import random
import pickle
import numpy as np
from scipy.optimize import fsolve
import warnings
class EC_cell:
	def __init__(self, cellIdent, loc, input_type_list, gridCell=False):
		
		self.ranGen=random.Random()	
		self.ranGen.seed(cellIdent)
		self.ID = cellIdent
		self.loc = loc
		#self.assign_param(input_type_list)
		self.gridCell = gridCell
		if gridCell == False:
			self.assign_param(input_type_list)
		else:
			self.rate_map = __import__('rate_map')
			self.coeffa = input_type_list[0]
			self.coeffd = input_type_list[1]

	def connect_pre(self, post, wt, dly):
		# Creates the NetCon object which will connect the EC cell to the 
		# synapse of the postsynaptic cell.
		#self.netcon = h.NetCon(self.stim, post, sec=self.stim)
		nc = h.NetCon(self.stim, post)
		nc.weight[0] = wt
		nc.delay = dly
		return nc
		
	def assign_param(self, input_type_list):
		if len(input_type_list) > 1:
			self.meanFreq = self.ranGen.gauss(input_type_list[2][0],input_type_list[2][1])
			if self.meanFreq < 0.0:
				self.meanFreq = 0.0
			
			self.initialDelay = self.ranGen.randint(input_type_list[4][0], input_type_list[4][1])
		if input_type_list[0] == 'NetStim':
			self.stim = h.NetStim()
			self.stim.interval = 1000./input_type_list[1]
			self.stim.number = 1.5*self.meanFreq/self.interval
			self.stim.start = 0
			self.stim.noise = input_type_list[3]
		
		elif input_type_list[0] == 'VecStim':
			if input_type_list[1] == 'Poisson':
				self.event_times = self.Poisson(self.meanFreq, input_type_list[3], self.initialDelay)
				#self.event_times = h.Vector()
				#for ii in range(len(self.spike_times)):
				#	self.event_times.append(self.spike_times[ii])
			
			if input_type_list[1] == 'Load':
				self.f = open(self.meanFreq)
				self.spike_times = pickle.load(f)
				self.f.close()
				
				self.event_times = h.Vector()
				for ii in range(len(spike_times)):
					self.event_times.append(spike_times[ii])
		
			if input_type_list[1] == 'Multi Poisson':
				self.Poisson_intervals = [ [] for i in range(self.meanFreq) ]
				for ii in range(self.meanFreq):
					self.Poisson(input_type_list[3][ii], self.initialDelay, input_type_list[5][ii])
					for jj in range(len(self.spike_times)):
						self.Poisson_intervals[ii].append(self.spike_times[jj])
				self.event_times = h.Vector()
				for ii in range(self.stim.meanFreq):
					for jj in range(len(self.Poisson_intervals[ii])):
						if ii < self.stim.meanFreq - 1:
							if self.Poisson_intervals[ii][jj] < input_type_list[5][ii+1]:
								self.event_times.append(self.Poisson_intervals[ii][jj])
							else:
								break
						else:
							if self.Poisson_intervals[ii][jj] >= input_type_list[5][ii]:
								self.event_times.append(self.Poisson_intervals[ii][jj])
			
			if input_type_list[1] == 'Input':
				self.event_times = h.Vector(self.stim.meanFreq)
									
			self.stim = h.VecStim()
			self.stim.play(self.event_times)
	
	def Poisson(self, frequency, tstop, delay, refractory=True):
		if frequency == 0:
			spike_times = [tstop+10]
		else:
			tres = 0.00025 # seconds
			tau = 0.035 # seconds
			#num_bins = int(tstop/tres)
			init_refract = -1e9
			spike_times = h.Vector()
			time = 0
			#for ii in range(num_bins):
			while time < tstop:
				roll = self.ranGen.uniform(0,1)
				refract = 1-np.exp(-(time-init_refract)/tau)
				if (roll <= frequency*tres*refract):
					spike_times.append(time*1000+delay)
					if refractory:
						init_refract = time+tres
				time += tres
			
			return spike_times
	
	def make_gridmap(self,loc,position,gridCoeffs):
		min_y = 0
		max_y = 2.4080397823282595
		y = max_y-loc[1]
		ysig = self.sig(y,min_y-0.32560712688527405,max_y,2.15,-1.7858005014930183,0.01)
		self.dist = gridCoeffs['slope']['gridSpacing']*(ysig-min_y) + gridCoeffs['intercept']['gridSpacing']
		self.area = 0.8*(gridCoeffs['slope']['fieldSizes']*(ysig-min_y) + gridCoeffs['intercept']['fieldSizes'])
		self.theta = 2*(self.dist - gridCoeffs['intercept']['orientationSpacing'])/gridCoeffs['slope']['orientationSpacing']
		self.theta += self.ranGen.gauss(0,15.49)
                #self.theta = self.ranGen.uniform(0,60)
		
                with warnings.catch_warnings():
                       warnings.simplefilter("ignore")
                       self.a, self.d = fsolve(self.gridParamEquations,(0.01,150),args=(self.area,self.dist))	
		
		self.a = np.abs(self.a)
		self.d = np.abs(self.d)
		#self.a = 4.2159442960765331e-05
                #self.d = 37.198982852044438
                
		self.x0 = self.ranGen.uniform(-self.d,self.d)
		self.y0 = self.ranGen.uniform(-self.d,self.d)
		c = np.array([[self.x0],[self.y0]])
		
		X = np.array(position['x'])
		Y = np.array(position['y'])
		rates = self.rate_map.output_rate(X,Y,self.a,c,self.d,self.theta)
		
		rates = np.array(rates) - np.min(rates)
		rates = rates/np.max(rates)*50
		rates[rates==0] += 0.1
		
		return rates
	
	def sig(self,x,A,K,B,M,v):
		numerator = K-A
		denominator = (1+np.exp(-B*(x-M)))
		return A + numerator/denominator**(1/v)
	
	def gridParamEquations(self,p,A,d):
		x, y = p
		x = np.abs(x)
		y = np.abs(y)
		a = np.array([1,x,y,x**2,y**2,x**2*y,x**2*y**2,y**2,x*y**2,x*y])
		a1 = a*self.coeffa
		a2 = a*self.coeffd
		return (sum(a1)-A,sum(a2)-d)
	
	def hetero_poisson(self,rates,rate_time,tau,refractory=True):
		self.stim = h.VecStim()
		
		# Convert to seconds
		tstop = np.max(rate_time) # seconds
		tres = 0.00025 # seconds
		
		# Interpolate frequencies to be the same resolution as time
		frequency = rates
		#num_bins = int(np.ceil(tstop/tres))
		#time = tres*np.arange(num_bins)
		#self.frequency = np.interp(time,rate_time,frequency)
		#interp = interp1d(rate_time,frequency,'linear')
		#self.frequency = interp(time)
		#del time
		# Generate spikes
		self.evec = h.Vector()
		init_refract = -1e9
		t = 0
		#ind = 0
		#for ii in range(num_bins):
                rate_res = rate_time[1]-rate_time[0]
		while t < tstop:
                        ind = int(t/rate_res)
			roll = self.ranGen.uniform(0,1)
			refract = 1-np.exp(-(t-init_refract)/tau)
			#refract = 1
			if (roll <= frequency[ind]*tres*refract):
				self.evec.append(t*1000)
				if refractory:
					init_refract = t+tres
			t += tres
			#ind += 1
			#if ind >= num_bins:
			#	ind = num_bins-1
		
		self.stim.play(self.evec)
        
        def hetero_poisson_old(self,rates,rate_time,tau,refractory=True):
                self.stim = h.VecStim()
                tstop = np.max(rate_time)
                frequency = rates
                mu = 1./frequency
                rmax = max(frequency)
                mumax = 1./rmax
                elapsed_time = 0
                self.thin_spikes = []
                while elapsed_time < tstop:
                        interval = -mumax*np.log(self.ranGen.uniform(0,1))
                        elapsed_time += interval
                        self.thin_spikes.append(elapsed_time)
                
                extra = self.thin_spikes.pop(-1)
                self.thin_spikes = np.array(self.thin_spikes)
                rn = [ self.ranGen.uniform(0,1) for ii in range(len(self.thin_spikes)) ]
                idx = np.searchsorted(rate_time,self.thin_spikes)-1
                spike_rate = frequency[idx]
                self.thin_spikes = self.thin_spikes[rn<spike_rate/rmax]
                self.evec = h.Vector()
                for spike in self.thin_spikes:
                        self.evec.append(1000*spike)
                
                self.stim.play(self.evec)
                
        def hetero_poisson2(self,rates,rate_time,tau,refractory=True):
                self.stim = h.VecStim()
                tstop = np.max(rate_time)
                rate_res = rate_time[1]-rate_time[0]
                frequency = np.copy(rates)
                mu = 1./frequency
                rmax = np.max(frequency)
                mumax = 1./rmax
                
                L = int(10*tau/rate_res)
                refract = 1-np.exp(-rate_time[:L]/tau)
                elapsed_time = 0
                pass_idx = []
                self.thin_spikes = []
                while elapsed_time < tstop:
                        interval = -mumax*np.log(self.ranGen.uniform(0,1))
                        elapsed_time += interval
                        self.thin_spikes.append(elapsed_time)
                
                self.thin_spikes = np.array(self.thin_spikes[:-1])
                rn = np.random.uniform(0,1,len(self.thin_spikes))
                idx = np.searchsorted(rate_time,self.thin_spikes)-1
                for ii in range(len(self.thin_spikes)):
                        if rn[ii] < frequency[idx[ii]]/rmax:
                                pass_idx.append(ii)
                                if (idx[ii]+L) > len(frequency):
                                        frequency[idx[ii]:] *= refract[:len(frequency)-idx[ii]]
                                else:
                                        frequency[idx[ii]:idx[ii]+L] *= refract
                
                self.thin_spikes = self.thin_spikes[np.array(pass_idx)]
                self.evec = h.Vector()
                for spike in self.thin_spikes:
                        self.evec.append(1000*spike)
                
                self.stim.play(self.evec)
