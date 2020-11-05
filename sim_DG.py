########################################################################
# Main --> DG simulation code
# 
# 2015.09.01
# Version 4 
# This version incorporates the new cell class structure Cell
# Also has many flags to specify connectivities
# Ability to save out somatic or segment voltages
#
# Old News:
# Version 2 includes basket cells
# 
########################################################################

import sys
sys.dont_write_bytecode = True

from neuron import h #This has everything we need for a NEURON-based simulation.
pc = h.ParallelContext()

import time # So we can time code execution.
import random # We want a pseudo-random number generator so that the simulations are repeatable.
import subprocess # For OS-level functions, including reads/writes to disk.
import cPickle # To easily allow us to write data structures to disk & read them back in again.
import numpy as np

for path in sys.path:
        if 'scratch' in path:
                scratchdir = path

START = time.time()

import ECcell_class_v5 as ECcell_class
import cell
from optParam_reduced import * 			# Load simulation parameters
from dirName import * 					# Load directory and connectivity information
h.xopen("ca3b.hoc")
# End module import

def euclidean(loc1,loc2):
	return ((loc1[0]-loc2[0])**2 + (loc1[1]-loc2[1])**2)**0.5

ca3type = 'ca3pyramidalcell_Bakerb'

############################
# Processing the NMDA flag #
############################
nmda_flag = int(sys.argv[15])

#############################
# Processing the Assoc flag #
#############################
if len(sys.argv) > 16:
	assoc_Factor = float(sys.argv[18])
	weight_dict[ca3type]['rad'] *= assoc_Factor
        weight_dict['granulecell']['LEC'] *= MEC_GC
	weight_dict['granulecell']['MEC'] *= MEC_GC
	weight_dict['granulecell']['BC'] *= BC_GC
	weight_dict['basketcell']['GC'] *= GC_BC
	weight_dict['basketcell']['MEC'] *= EC_BC
	weight_dict['basketcell']['LEC'] *= EC_BC

#############################
# Processing grid cell flag #
#############################
if grid_flag:
	# File containing the virtual mouse position data
	# generated from random_walk.py
	#with open('rat_position_2017.pickle') as f:
        with open('rat_position_80cm.pickle') as f:
        #with open('rat_position_brown2.pickle') as f:
		position = cPickle.load(f)
        
        position['x'] = np.array(position['x'])
        position['y'] = np.array(position['y'])
        position['time'] = np.array(position['time'])
        '''
        pos_tstop = tstop/1000
        tres = 0.1
        num_bins = int(np.ceil(pos_tstop/tres))
        t_interp = tres*np.arange(num_bins)
        pos_interp = {}
        pos_interp['time'] = t_interp
        pos_interp['x'] = np.interp(t_interp,position['time'],position['x'])
        pos_interp['y'] = np.interp(t_interp,position['time'],position['y'])
        '''
        with open('gridMapCoefficients.pickle') as f:
                gridMapCoeffs = cPickle.load(f)
	
	# Contains grid cell parameters generated from
	# gridCellParam.py
	#with open('gridParam2.pickle') as f:
	with open('gridCellTopography.pickle') as f:
		gridMapTopography = cPickle.load(f)

#############################
# Processing model flag #
#############################
if model_flag:
	model = 'Single'
else:
	model = 'Multi'

h.load_file("stdrun.hoc")

# Initialize the random number generator.
ranGen = random.Random()
ranGen.seed(299792458)

MECCells = []
LECCells = []
postCells = []
basketCells = []
CA3pyramidalCells = []
CA3intCells = []

net_cons = []

# Write information about the simulation to the output file
if pc.id() == 0:
        sim_files = subprocess.check_output(['ls']).split('\n')[:-1]
        sim_files = [fname for fname in sim_files if not fname.startswith('Con')]
        filetypes = ['mod','py','pickle','hoc','swc']
	_=subprocess.call(['mkdir',outputDir])
        for fname in sim_files:
                for ftype in filetypes:
                        if fname.endswith(ftype):
                                _=subprocess.call(['cp',fname,outputDir])
        
        print "%f, %f, %f"%(weight_dict['basketcell']['GC'],weight_dict['basketcell']['MEC'],weight_dict['granulecell']['BC'])
        print "Number of MEC cells: " + str(N_MEC) + "."
	print "Number of LEC cells: " + str(N_LEC) + "."
	print "Number of Granule cells: " + str(N_GC) + "."
	print "Number of Basket cells: " + str(N_BC) + "."
	print "Simulation time: " + str(tstop) + "."
	print time.time() - START
	
	with open(outputDir + '/system_params.txt','w') as f:
		f.write('System Connectivity Configuration\n')
		f.write('ECtoGC ' + str(ECtoGC)+'\n')
		f.write('ECtoCA3 ' + str(ECtoCA3)+'\n')
		f.write('ECtoBC ' + str(ECtoBC)+'\n')
		f.write('GCtoBC ' + str(GCtoBC)+'\n')
		f.write('GCtoCA3 ' + str(GCtoCA3)+'\n')
		f.write('CA3toCA3 ' + str(CA3toCA3)+'\n')
		f.write('CA3inhibition ' + str(CA3inhibition)+'\n')
		f.write('NMDA ' + str(nmda_flag)+'\n')
		f.write('Slice ' + str(slice_flag)+'\n')
		f.write('Save Voltages ' + str(volt_flag)+'\n')
		f.write('Grid Cells ' + str(grid_flag)+'\n')
		f.write('Single Compartment ' + str(model_flag))

#############################
# Load in connectivity data #
#############################
filedir = '/staging/tb/geneyu/Con500_v3'
filedir = scratchdir
#with open("/scratch/simulation/newCon%.3d.pickle" % int(pc.id())) as f:
with open(filedir+'/Con%.3d.pickle' % int(pc.id())) as f:
#with open(filedir+"/reducedCon%i.pickle" % int(pc.id())) as f:
	Conv = cPickle.load(f)

##############################################
# Load in granule cell and CA3 cell location #
# and synapse data                           #
##############################################
with open(filedir+"/data_cellLocations.pickle") as f:
	locs = cPickle.load(f)

with open(filedir+"/data_IDranges.pickle") as f:
        IDranges = cPickle.load(f)

with open(filedir+"/data_MFpathlength.pickle") as f:
        MFpathlength = cPickle.load(f)

with open(filedir+"/data_CAregions.pickle") as f:
        CA3regions = cPickle.load(f)

with open(filedir+"/data_numSynCA3.pickle") as f:
        numSynCA3 = cPickle.load(f)

if pc.id() == 0:
        print 'Files loaded'

IDs = np.array(Conv.keys())
idx = np.all([IDs>=IDranges['DG']['GC'][0],IDs<=IDranges['DG']['GC'][1]],axis=0)
IDs_DG_GC = IDs[idx]

idx = np.all([IDs>=IDranges['DG']['BC'][0],IDs<=IDranges['DG']['BC'][1]],axis=0)
IDs_DG_BC = IDs[idx]

idx = np.all([IDs>=IDranges['CA3']['PC'][0],IDs<=IDranges['CA3']['PC'][1]],axis=0)
IDs_CA3_PC = IDs[idx]

idx = np.all([IDs>=IDranges['CA3']['BC'][0],IDs<=IDranges['CA3']['BC'][1]],axis=0)
IDs_CA3_BC = IDs[idx]

###############################
# Create the Entorhinal Cells #
###############################
MECCells = []
LECCells = []
# Need to convert tstop into seconds
#input_type_list = ['VecStim','Poisson', (11.27,1.13), tstop/1000, (0,0)]
for i in range(int(pc.id()), N_MEC, int(pc.nhost())):
	if grid_flag:
                input_type_list = [gridMapCoeffs['area coeffs'],gridMapCoeffs['dist coeffs']]
		MECCells.append(ECcell_class.EC_cell(i, locs['EC']['MEC'][i], input_type_list, grid_flag))
		if MECCells[-1].gridCell:
			rates = MECCells[-1].make_gridmap(MECCells[-1].loc,position,gridMapTopography)
			MECCells[-1].hetero_poisson2(rates,position['time'],0.035)
	else:
                input_type_list = ['VecStim','Poisson', (5.0,0.01), tstop/1000, (0,0)]
		MECCells.append(ECcell_class.EC_cell(i, locs['EC']['MEC'][i], input_type_list))
        
        #if (i % 1000) == 0:
        #        with open(outputDir+'/%i.txt'%i,'w') as f:
        #                pass
        
	pc.set_gid2node(i,pc.id()) # This associates this cell with this node
	pc.cell(i,MECCells[-1].connect_pre(None,0,0)) 

if pc.id() == 0:
        print 'MEC made'

for i in range(int(pc.id()), N_LEC, int(pc.nhost())):
	input_type_list = ['VecStim','Poisson', (5.0,0.01), tstop/1000, (0,0)]
	LECCells.append(ECcell_class.EC_cell(i+N_MEC, locs['EC']['LEC'][i+N_MEC], input_type_list))
        #if (i % 1000) == 0:
        #        with open(outputDir+'/%i.txt'%(i+N_MEC),'w') as f:
        #                pass
        
	pc.set_gid2node(i+N_MEC,pc.id()) # This associates this cell with this node
	pc.cell(i+N_MEC,LECCells[-1].connect_pre(None,0,0)) 

if pc.id() == 0:
        print 'LEC made'

##############################
# Create the Granule Cells   #
##############################
# synapse count structure:
# 0: MEC
# 1: LEC
# 2: BC
if GC_flag == 1:
	if GC_spike_flag == 0:
		if GC_rand_flag == 0:
			for ID in IDs_DG_GC:
				fileName = "output0_updated.swc"
				location = locs['DG']['GC'][ID]
				connect_flag = 0
				if slice_flag == 0:
					postCells.append(cell.Cell(ID, location, synvars, 'granulecell', fileName,model))
					connect_flag = 1
				else:
					if Lower <= location[1] < Upper:
						postCells.append(cell.Cell(ID, location, synvars, 'granulecell', fileName,model))
						connect_flag = 1
				if connect_flag:
					pc.set_gid2node(ID, pc.id())
					pc.cell(ID,postCells[-1].connect_pre(None,0,0))
		
		# If GC_rand_flag is set to one, then treat the granule cells as independent,
		# random Poisson spike generators
		elif GC_rand_flag == 1:
			for ID in IDs_DG_GC:
				if locs['DG']['GC'][ID][0] < 0:
					input_type_list = ['VecStim','Poisson', (0.8,0.1), tstop/1000, (0,0)]
				else:
					input_type_list = ['VecStim','Poisson', (0.8,0.1), tstop/1000, (0,0)]
				postCells.append(ECcell_class.EC_cell(ID, locs['DG']['GC'][ID], input_type_list))
				pc.set_gid2node(ID,pc.id()) # This associates this cell with this node
				pc.cell(ID,postCells[-1].connect_pre(None,0,0))
	
	# if GC_spike_flag is set to one, then load in the GC spike
	# times and use them as spike generators instead of creating cells
	# to reduce computation time
	elif GC_spike_flag == 1:
		GC_evecs = {}
		with open(scratchdir+'/spikeTimes', 'rb') as f:
			tmp = cPickle.load(f)
		
		for ID in IDs_DG_GC:
			if ID in tmp:
				GC_evecs[ID] = h.Vector(tmp[ID])
			else:
				GC_evecs[ID] = h.Vector([])
				
		
		del tmp
                
		GC_stim = {}
		for ID in GC_evecs:
			GC_stim[ID] = h.VecStim()
			GC_stim[ID].play(GC_evecs[ID])
			pc.set_gid2node(ID, pc.id())
			
			# Doing the whole connect_pre thing
			nc = h.NetCon(GC_stim[ID],None)
			nc.weight[0] = 0
			nc.delay = 0
			
			pc.cell(ID,nc)
			net_cons.append(nc)

if pc.id() == 0:
        print 'GC made'

##############################
# Create the Basket Cells   #
##############################
# synapse count structure:
# 0: MEC
# 1: LEC
# 2: GC
if BC_flag == 1:
	for ID in IDs_DG_BC:
		location = locs['DG']['BC'][ID]
		connect_flag = 0
		if slice_flag == 0:
			basketCells.append(cell.Cell(ID, location, synvars, 'basketcell_DG', '',model))
			connect_flag = 1
		else:
			if Lower <= location[1] < Upper:
				basketCells.append(cell.Cell(ID, location, synvars, 'basketcell_DG', '',model))
				connect_flag = 1
		if connect_flag:
			pc.set_gid2node(ID, pc.id())
			pc.cell(ID,basketCells[-1].connect_pre(None,0,0))

if pc.id() == 0:
        print 'DG BC made'

##################################
# Create the CA3 pyramidal cells #
##################################
# synapse count structure
# 0: GC lucidum
# 1: radiatum
# 2: oriens
# 3: GC basal
# 4: MEC
# 5: LEC
if CA3_flag == 1:
	for ID in IDs_CA3_PC:
		fileName = "CA3output0_updated.swc"
		location = locs['CA3']['PC'][ID]
		connect_flag = 0
		if slice_flag == 0:
			CA3pyramidalCells.append(cell.Cell(ID, location, synvars, ca3type, fileName,model))
			connect_flag = 1
		else:
			if Lower <= location[1] < Upper:
				CA3pyramidalCells.append(cell.Cell(ID, location, synvars, ca3type, fileName,model))
				connect_flag = 1
		if connect_flag:
                        CA3pyramidalCells[-1].region = CA3regions['CA3']['PC'][ID]
			pc.set_gid2node(ID,pc.id()) # This associates this cell
			pc.cell(ID, CA3pyramidalCells[-1].connect_pre(None,0,0))

if pc.id() == 0:
        print 'Number of CA3 PC: %i'%len(CA3pyramidalCells)
        print 'CA3 PC made'

###############################
# Create the CA3 interneurons #
###############################
if CA3inhibition == 1:
	for ID in IDs_CA3_BC:
		location = locs['CA3']['BC'][ID]
		connect_flag = 0
		if slice_flag == 0:
			CA3intCells.append(cell.Cell(ID, location, synvars, 'basketcell_CA3', '',model))
			connect_flag = 1
		else:
			if Lower <= location[1] < Upper:
				CA3intCells.append(cell.Cell(ID, location, synvars, 'basketcell_CA3', '',model))
				connect_flag = 1
		if connect_flag:
                        CA3intCells[-1].region = CA3regions['CA3']['BC'][ID]
			pc.set_gid2node(ID,pc.id())
			pc.cell(ID, CA3intCells[-1].connect_pre(None,0,0))

if pc.id() == 0:
        print 'CA3 BC made'

######################################
# Now, we need to connect the cells. #
# Start with EC-Granule projection   #
######################################
ratio = 0.0775389348453*0.872
tau1 = 20.3877436154
tau2 = 26.6830234133
tau3 = 158.729359569
wtau2 = 0.963468100127
wtau3 = 1-wtau2

t0 = np.linspace(0,100,10000)
left = -np.log(tau1)-t0/tau1
right = np.log(wtau2/tau2*np.exp(-t0/tau2) + wtau3/tau3*np.exp(-t0/tau3))
tp = t0[np.argmin(np.abs(left-right))]
factor = -np.exp(-tp/tau1) + wtau2*np.exp(-tp/tau2) + wtau3*np.exp(-tp/tau3)
factor = 1/factor

if ECtoGC == 1:
	for iCell in postCells:
		# This is the set of synapses in the middle 1/3, where MEC cells project to
                Syn = Conv[iCell.ID][Conv[iCell.ID] < 66000]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			# Need to choose a segment in the right layer
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['middleThird'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['middleThird'][choice])
			iCell.synGroups['AMPA']['middleThird'][choice].tau2 = tau_dict['granulecell']['MEC']
                        iCell.synGroups['AMPA']['middleThird'][choice].tau1 = 0.5
                        iCell.synGroups['AMPA']['middleThird'][choice].e = 0
			delay = baseSynDelay_MEC_GC+1.7+euclidean(iCell.location,(0,locs['DG']['MEC'][inCellNum]))/MECvel
                        nc.weight[0] = weight_dict['granulecell']['MEC']
			nc.delay = delay
			net_cons.append(nc)
                        
                        if nmda_flag:
                                nc = pc.gid_connect(inCellNum, iCell.synGroups['NMDA']['middleThird'][choice])
                                iCell.synGroups['NMDA']['middleThird'][choice].tau1 = tau1
                                iCell.synGroups['NMDA']['middleThird'][choice].tau2 = tau2
                                iCell.synGroups['NMDA']['middleThird'][choice].tau3 = tau3
                                iCell.synGroups['NMDA']['middleThird'][choice].wtau2 = wtau2
                                iCell.synGroups['NMDA']['middleThird'][choice].factor = factor
                                iCell.synGroups['NMDA']['middleThird'][choice].e = 0
                                nc.weight[0] = weight_dict['granulecell']['MEC']*ratio
                                nc.delay = delay
                                net_cons.append(nc)
                        
		
		# This is the set of synapses in the outer 1/3, where LEC cells project to
                idx = np.all([Conv[iCell.ID]>=66000,Conv[iCell.ID]<112000],axis=0)
                Syn = Conv[iCell.ID][idx]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['outerThird'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['outerThird'][choice])
			#iCell.synGroups['AMPA']['outerThird'][choice].tau1 = tau_dict['granulecell']['LEC'][0]
			iCell.synGroups['AMPA']['outerThird'][choice].tau2 = tau_dict['granulecell']['LEC']
                        iCell.synGroups['AMPA']['outerThird'][choice].tau1 = 0.5
                        iCell.synGroups['AMPA']['outerThird'][choice].e = 0
			delay = baseSynDelay_LEC_GC+1.7+euclidean(iCell.location,(0,locs['DG']['LEC'][inCellNum]))/LECvel
                        nc.weight[0] = weight_dict['granulecell']['LEC']
			nc.delay = delay
			net_cons.append(nc)
                        
                        if nmda_flag:
                                nc = pc.gid_connect(inCellNum, iCell.synGroups['NMDA']['outerThird'][choice])
                                iCell.synGroups['NMDA']['outerThird'][choice].tau1 = tau1
                                iCell.synGroups['NMDA']['outerThird'][choice].tau2 = tau2
                                iCell.synGroups['NMDA']['outerThird'][choice].tau3 = tau3
                                iCell.synGroups['NMDA']['outerThird'][choice].wtau2 = tau2
                                iCell.synGroups['NMDA']['outerThird'][choice].factor = factor
                                iCell.synGroups['NMDA']['outerThird'][choice].e = 0
                                nc.weight[0] = weight_dict['granulecell']['LEC']*ratio
                                nc.delay = delay
                                net_cons.append(nc)

if pc.id() == 0:
        print 'EC-GC connected'

####################################
# Feedforward Inhibition: EC to BC #
####################################
if ECtoBC == 1:
	for iCell in basketCells:
		# This is the set of synapses in the middle 1/3, where MEC cells project to
                Syn = Conv[iCell.ID][Conv[iCell.ID]<66000]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
                        nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][0])
			iCell.synGroups['AMPA']['soma'][0].tau2 = tau_dict['basketcell']['MEC']
                        nc.weight[0] = weight_dict['basketcell']['MEC']
                
                idx = np.all([Conv[iCell.ID]>=66000,Conv[iCell.ID]<112000],axis=0)
                Syn = Conv[iCell.ID][idx]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][0])
			iCell.synGroups['AMPA']['soma'][0].tau2 = tau_dict['basketcell']['LEC']
                        nc.weight[0] = weight_dict['basketcell']['LEC']
			nc.delay = baseSynDelay_LEC_GC+1.7+euclidean(iCell.location,(0,locs['DG']['LEC'][inCellNum]))/LECvel
			net_cons.append(nc)

if pc.id() == 0:
        print 'EC-DGBC connected'

#################################
# Feedback Inhibition: GC to BC #
#################################
if GCtoBC == 1:
	for iCell in basketCells:
                idx = np.all([Conv[iCell.ID]>=IDranges['DG']['GC'][0],Conv[iCell.ID]<=IDranges['DG']['GC'][1]],axis=0)
                Syn = Conv[iCell.ID][idx]
                for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][1])
			iCell.synGroups['AMPA']['soma'][1].tau2 = tau_dict['basketcell']['GC']
                        nc.weight[0] = weight_dict['basketcell']['GC']
			nc.delay = baseSynDelay_MF+euclidean(iCell.location,locs['DG']['GC'][inCellNum])/MFvel
			net_cons.append(nc)

if pc.id() == 0:
        print 'GC-DGBC connected'

#################################
# Feedback Inhibition: BC to GC #
#################################
if BCtoGC == 1:
	for iCell in postCells:
                idx = np.all([Conv[iCell.ID]>=IDranges['DG']['BC'][0],Conv[iCell.ID]<=IDranges['DG']['BC'][1]],axis=0)
                Syn = Conv[iCell.ID][idx]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
                        choice = iCell.ranGen.randint(0,len(iCell.synGroups['GABA']['soma'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['GABA']['soma'][choice])
			nc.weight[0] = weight_dict['granulecell']['BC']
			iCell.synGroups['GABA']['soma'][choice].tau2 = tau_dict['granulecell']['BC']
			iCell.synGroups['GABA']['soma'][choice].e = -75
			nc.delay = baseSynDelay_BC_GC+euclidean(iCell.location,locs['DG']['BC'][inCellNum])/BCvel
			net_cons.append(nc)

if pc.id() == 0:
        print 'DGBC-GC connected'

##########################################################
# Let's connect some mossy fibers to some CA3 pyramidals #
##########################################################
for iCell in CA3pyramidalCells:
	if GCtoCA3 == 1:
                idx = np.all([Conv[iCell.ID]>=IDranges['DG']['GC'][0],Conv[iCell.ID]<=IDranges['DG']['GC'][1]],axis=0)
                Syn = Conv[iCell.ID][idx]
                if iCell.region == 'c':
                        Syn = [ID for ID in Syn if locs['DG']['GC'][ID][0] >= 0]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			distance = MFpathlength[iCell.ID][inCellNum]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['lucidum'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['lucidum'][choice])
			iCell.synGroups['AMPA']['lucidum'][choice].tau2 = tau_dict[ca3type]['GC']
			nc.weight[0] = weight_dict[ca3type]['GC']
			nc.delay = baseSynDelay_MF + distance/MFvel
			net_cons.append(nc)
		
		##################################
		# Some basal dendrite action too #
		##################################
                Syn = Conv[iCell.ID][idx]
                if iCell.region == 'c':
                        Syn = [ID for ID in Syn if locs['DG']['GC'][ID][0] < 0]
                else:
                        Syn = []
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			distance = MFpathlength[iCell.ID][inCellNum]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['oriensProximal'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['oriensProximal'][choice])
			iCell.synGroups['AMPA']['oriensProximal'][choice].tau2 = tau_dict[ca3type]['GCbasal']
			nc.weight[0] = weight_dict[ca3type]['GCbasal']
			nc.delay = baseSynDelay_MF + distance/MFvel
			net_cons.append(nc)

if pc.id() == 0:
        print "Mossy-CA3 connected"
        print "Number net_cons: %i"%len(net_cons)

for iCell in CA3pyramidalCells:	
	############################################
	# Connecting EC to CA3 pyramidal cells now #
	############################################
	if ECtoCA3 == 1:
                Syn = Conv[iCell.ID][Conv[iCell.ID]<66000]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['lacunosumMEC'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['lacunosumMEC'][choice])
			iCell.synGroups['AMPA']['lacunosumMEC'][choice].tau2 = tau_dict[ca3type]['MEC']
			nc.weight[0] = weight_dict[ca3type]['MEC']
			nc.delay = baseSynDelay_MEC_GC+1.7+2+euclidean(iCell.location,(0,locs['DG']['MEC'][inCellNum]))/MECvel
			net_cons.append(nc)
		
                idx = np.all([Conv[iCell.ID]>=66000,Conv[iCell.ID]<112000],axis=0)
                Syn = Conv[iCell.ID][idx]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['lacunosumLEC'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['lacunosumLEC'][choice])
			iCell.synGroups['AMPA']['lacunosumLEC'][choice].tau2 = tau_dict[ca3type]['LEC']
			nc.weight[0] = weight_dict[ca3type]['LEC']
			nc.delay = baseSynDelay_LEC_GC+1.7+2+euclidean(iCell.location,(0,locs['DG']['LEC'][inCellNum]))/LECvel
			net_cons.append(nc)

if pc.id() == 0:
        print "Number net_cons: %i" % len(net_cons)
        print "EC-CA3 connected"

for iCell in CA3pyramidalCells:
	###############################################
	# Connecting auto-associative radiatum spines #
	###############################################
	if CA3toCA3 == 1:
                idx = np.all([Conv[iCell.ID]>=IDranges['CA3']['PC'][0],Conv[iCell.ID]<=IDranges['CA3']['PC'][1]],axis=0)
                Syns = Conv[iCell.ID][idx]
                Syn = Syns[:numSynCA3['rad'][iCell.ID]]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['radiatum'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['radiatum'][choice])
			iCell.synGroups['AMPA']['radiatum'][choice].tau2 = tau_dict[ca3type]['rad']
			nc.weight[0] = weight_dict[ca3type]['rad']
			nc.delay = baseSynDelay_MF+euclidean(locs['CA3']['PC'][iCell.ID],locs['CA3']['PC'][inCellNum])/Assvel
			net_cons.append(nc)
	
		#############################################
		# Connecting auto-associative oriens spines #
		#############################################
                Syn = Syns[numSynCA3['rad'][iCell.ID]:]
		for iSyn in range(len(Syn)):
			inCellNum = Syn[iSyn]
			choice = iCell.ranGen.randint(0,len(iCell.synGroups['AMPA']['oriensDistal'])-1)
			nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['oriensDistal'][choice])
			iCell.synGroups['AMPA']['oriensDistal'][choice].tau2 = tau_dict[ca3type]['ori']
			nc.weight[0] = weight_dict[ca3type]['ori']
			nc.delay = baseSynDelay_MF+euclidean(locs['CA3']['PC'][iCell.ID],locs['CA3']['PC'][inCellNum])/Assvel
			net_cons.append(nc)

if pc.id() == 0:
        print 'Excitatory-CA3 connected'
	
######################################################
# Connecting CA3 interneurons to CA3 pyramidal cells #
######################################################
# 0.075 probability to be in the lucidum/pyramidal layer
# 0.6 probability to be in the radiatum
# 0.075 probability to be in the lacunosem-moleculare
# 0.25 probability to be in the oriens

if CA3inhibition == 1:
	for iCell in CA3pyramidalCells:
                idx = np.all([Conv[iCell.ID]>=IDranges['CA3']['BC'][0],Conv[iCell.ID]<=IDranges['CA3']['BC'][1]],axis=0)
                Syn = Conv[iCell.ID][idx]
		for inCellNum in Syn:
                        choice = iCell.ranGen.randint(0,len(iCell.synGroups['GABA']['soma'])-1)
                        nc = pc.gid_connect(inCellNum, iCell.synGroups['GABA']['soma'][choice])
                        iCell.synGroups['GABA']['soma'][choice].tau2 = tau_dict[ca3type]['BC']
                        iCell.synGroups['GABA']['soma'][choice].e = -75
                        nc.weight[0] = weight_dict[ca3type]['BC']
                        nc.delay = baseSynDelay_MF+euclidean(locs['CA3']['PC'][iCell.ID],locs['CA3']['BC'][inCellNum])/Assvel
                        net_cons.append(nc)

if pc.id() == 0:
        print 'CA3 BC-PC connected'

################################################
# Now connecting input to the CA3 interneurons #
################################################
if CA3inhibition  == 1:
	for iCell in CA3intCells:
		if GCtoCA3 == 1:
                        idx = np.all([Conv[iCell.ID]>=IDranges['DG']['GC'][0],Conv[iCell.ID]<=IDranges['DG']['GC'][1]],axis=0)
                        Syn = Conv[iCell.ID][idx]
                        if iCell.region == 'c':
                                Syn = [ID for ID in Syn if locs['DG']['GC'][ID][0] >= 0]
                        
			for iSyn in range(len(Syn)):
				inCellNum = Syn[iSyn]
				distance = MFpathlength[iCell.ID][inCellNum]
				nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][0])
				iCell.synGroups['AMPA']['soma'][0].tau2 = tau_dict['ca3interneuron']['GC']
				nc.weight[0] = weight_dict['ca3interneuron']['GC']
				nc.delay = baseSynDelay_MF + distance/MFvel
				net_cons.append(nc)
                        
                        Syn = Conv[iCell.ID][idx]
                        if iCell.region == 'c':
                                Syn = [ID for ID in Syn if locs['DG']['GC'][ID][0] < 0]
                        else:
                                Syn = []
                        
                        for iSyn in range(len(Syn)):
                                inCellNum = Syn[iSyn]
                                distance = MFpathlength[iCell.ID][inCellNum]
                                nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][1])
                                iCell.synGroups['AMPA']['soma'][1].tau2 = tau_dict['ca3interneuron']['GCbasal']
                                nc.weight[0] = weight_dict['ca3interneuron']['GCbasal']
                                nc.delay = baseSynDelay_MF + distance/MFvel
                                net_cons.append(nc)
                
		if CA3toCA3 == 1:
                        idx = np.all([Conv[iCell.ID]>=IDranges['CA3']['PC'][0],Conv[iCell.ID]<=IDranges['CA3']['PC'][1]],axis=0)
                        Syn = Conv[iCell.ID][idx]
			for iSyn in range(len(Syn)):
				inCellNum = Syn[iSyn]
				nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][2])
				iCell.synGroups['AMPA']['soma'][2].tau2 = tau_dict['ca3interneuron']['rad']
				nc.weight[0] = weight_dict['ca3interneuron']['rad']
				nc.delay = baseSynDelay_MF+euclidean(locs['CA3']['BC'][iCell.ID],locs['CA3']['PC'][inCellNum])/Assvel
				net_cons.append(nc)
                
		if ECtoCA3 == 1:
                        Syn = Conv[iCell.ID][Conv[iCell.ID]<66000]
			for iSyn in range(len(Syn)):
				inCellNum = Syn[iSyn]
				nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][3])
				iCell.synGroups['AMPA']['soma'][3].tau2 = tau_dict['ca3interneuron']['MEC']
				nc.weight[0] = weight_dict['ca3interneuron']['MEC']
				nc.delay = baseSynDelay_LEC_GC+1.7+2+euclidean(iCell.location,(0,locs['DG']['MEC'][inCellNum]))/MECvel
				net_cons.append(nc)
                        
                        idx = np.all([Conv[iCell.ID]>=66000,Conv[iCell.ID]<112000],axis=0)
                        Syn = Conv[iCell.ID][idx]
			for iSyn in range(len(Syn)):
				inCellNum = Syn[iSyn]
				nc = pc.gid_connect(inCellNum, iCell.synGroups['AMPA']['soma'][3])
				iCell.synGroups['AMPA']['soma'][3].tau2 = tau_dict['ca3interneuron']['LEC']
				nc.weight[0] = weight_dict['ca3interneuron']['LEC']
				nc.delay = baseSynDelay_LEC_GC+1.7+2+euclidean(iCell.location,(0,locs['DG']['LEC'][inCellNum]))/LECvel
				net_cons.append(nc)		

if pc.id() == 0:
        print 'CA3 BC input connected'

# Instrument cells to record spikes
tvec = h.Vector()
idvec = h.Vector()
'''
for i in MECCells:
	nc = i.connect_pre(None,0,0)
	nc.record(tvec, idvec, i.ID)
for i in LECCells:
	nc = i.connect_pre(None,0,0)
	nc.record(tvec, idvec, i.ID)
'''
for i in postCells:
	nc = i.connect_pre(None,0,0)
	nc.record(tvec, idvec, i.ID)
for i in basketCells:
	nc = i.connect_pre(None,0,0)
	nc.record(tvec, idvec, i.ID)
for i in CA3pyramidalCells:
	nc = i.connect_pre(None,0,0)
	nc.record(tvec, idvec, i.ID)
for i in CA3intCells:
	nc = i.connect_pre(None,0,0)
	nc.record(tvec, idvec, i.ID)

##############################################
# Instrument cells to record somatic voltage #
##############################################
if volt_flag:
	# Force a default flag so it won't error
	flags = ['soma','seg']
	if volt_flag not in flags:
		volt_flag = 'soma'
	
	if volt_flag == 'soma':
		# Recording only at the soma
		voltages = {}
		for cell in postCells:
			voltages[cell.ID] = h.Vector()
			voltages[cell.ID].record(cell.c.soma[0](0.5)._ref_v,0.25)
                
		for cell in basketCells:
			voltages[cell.ID] = h.Vector()
			voltages[cell.ID].record(cell.soma(0.5)._ref_v,0.25)	
                
		for cell in CA3pyramidalCells:
			voltages[cell.ID] = h.Vector()
			voltages[cell.ID].record(cell.soma(0.5)._ref_v,0.25)
	
	elif volt_flag == 'seg':
		# Recording at every segment
		vflag = True
		iflag = False
		timestep = 0.25
		output_data = {}
		if pc.id() == 0:
			for ii in range(len(postCells)):
				if (ii%10) == 0:
					cell = postCells[ii]
					output_data[cell.ID] = {}
					data = cell.setdata(timestep)
					output_data[cell.ID]['v'] = data['v']
					output_data[cell.ID]['loc'] = data['loc']
			for ii in range(len(CA3pyramidalCells)):
				if (ii%10) == 0:
					cell = CA3pyramidalCells[ii]
					output_data[cell.ID] = {}
					data = cell.setdata(timestep)
					output_data[cell.ID]['v'] = data['v']
					output_data[cell.ID]['loc'] = data['loc']

######################
# Run the simulation #
######################
#if pc.id() == 0:
#voltage = h.Vector()
#voltage.record(CA3pyramidalCells[0].soma(0.5)._ref_v)
voltages = {}
for cell in postCells:
        voltages[cell.ID] = h.Vector()
        voltages[cell.ID].record(cell.soma(0.5)._ref_v,0.25)

pc.set_maxstep(10.0)

if pc.id() == 0:
        print "Max Step Set"

h.stdinit()
#h.cvode.queue_mode(1,0)
h.celsius = 35.0
tfin = 0
if pc.id() == 0:
	print "Initialization Complete, Simulation Starting..."

pc.barrier()
tStart = time.time()
for state in range(num_states):
	tvec.resize(0)
	idvec.resize(0)
	tfin += save_interval
	h.continuerun(tfin)
	
	spikes = {}
	for ii in range(len(tvec)):
		try:
			spikes[idvec[ii]].append(tvec[ii])
		except KeyError:
			spikes[idvec[ii]] = [tvec[ii]]
	
	for ID in spikes:
		spikes[ID] = np.array(spikes[ID])
	
	#time.sleep(int(pc.id())*0.1) # File writing issues; delays seem to work
	with open(scratchdir+"/spikeTimes"+str(int(pc.id())),'wb') as f:
		cPickle.dump(spikes,f,protocol=-1)
	
        #np.save(scratchdir+'/voltage%i.npy'%int(pc.id()),np.array(voltage))
        
	pc.barrier()
	
	if pc.id() == 0:
		files = subprocess.check_output(['ls',scratchdir]).split('\n')[:-1]
		files = [ fname for fname in files if fname.startswith('spikeTimes') ]
		combined_spikes = {}
		for fname in files:
			with open(scratchdir+'/'+fname,'rb') as f:
				tmp = cPickle.load(f)
		
			combined_spikes = dict(combined_spikes.items()+tmp.items())
			del tmp
	
		with open(outputDir+'/spikeTimes-'+str(state),'wb') as f:
			cPickle.dump(combined_spikes,f,protocol=-1)
		
		del combined_spikes
        
	pc.barrier()

for ID in voltages:
        voltages[ID] = np.array(voltages[ID])

with open(outputDir+'/voltages%i.pickle'%int(pc.id()),'wb') as f:
        cPickle.dump(voltages,f,protocol=-1)

#voltages = {}
#for ii in range(int(pc.nhost())):
#       voltages[ii] = np.load(scratchdir+'/voltage%i.npy'%ii)
#with open(outputDir+'/voltages-%i.pickle'%int(state),'wb') as f:
#        cPickle.dump(voltages,f,protocol=-1)

tEnd = time.time()
if pc.id() == 0:
        print "Wall-time (simulation only) = " + str(tEnd-tStart) + ".\n"
        #voltage = np.array(voltage)
        #np.save(outputDir+'/voltage.npy',voltage)

pc.barrier()

# Combine spike time file
if pc.id() == 0:
        spikeTimes = {}
        for state in range(num_states):
                with open(outputDir+'/spikeTimes-'+str(state)) as f:
                        tmp = cPickle.load(f)
                
                for ID in tmp:
                        try:
                                spikeTimes[ID] = np.concatenate([spikeTimes[ID],tmp[ID]],axis=0)
                        except KeyError:
                                spikeTimes[ID] = tmp[ID]
                
                del tmp
                
                # Delete intermediate file
                _=subprocess.call(['rm',outputDir+'/spikeTimes-'+str(state)])
        
        with open(outputDir+'/spikeTimes','wb') as f:
                cPickle.dump(spikeTimes,f,protocol=-1)

pc.barrier()

pc.done()
# End file
