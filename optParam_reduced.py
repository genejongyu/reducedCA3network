# This file contains all relevant parameters for the EPSP testing code

N_MEC = 66000
N_LEC = 46000
N_GC = 120000
N_BC = 894
N_CA3 = 25000
N_CA3int = 10000

tstop = 5000
save_interval = 1000
num_states = int(tstop/save_interval)

########################
# Define Synaptic Type #
########################
synvars = {}
norm_syntype = ['E2','E2_Prob','STDPE2_Clo','E2_STP_Prob','E2-NMDA2']
plas_syntype = ['STDPE2','STDPE2_STP','STDPE2_Prob']
synvars['type'] = 'E2'

#############################################
# Define physical dimensions of hippocampus #
#############################################
L = 10.0

######################################################
# Define propagation velocities of action potentials #
######################################################
# Differential Conduction Velocities in Perforant Path Fibres in Guinea Pig
# Tielen, Lopes da Silva, Mollevanger. 1981
LECvel = 0.32 # (mm/ms)
MECvel = 0.32 # (mm/ms)

# High Threshold, Proximal Initiation and Slow Conduction Velocity
# of Action potentials in Dentate Granule Neuron Mossy Fibers
MFvel = 0.27 # (mm/ms)

# The hippocampal lamella hypothesis revisited
Assvel = 0.39 # (mm/ms)

BCvel = 0.3 # Estimate but all velocities seem to be around 0.3 anyways

###################################
# Transmission delays at synapses #
###################################
baseSynDelay_MEC_GC = 0.7       # (ms)
baseSynDelay_MEC_BC = 0.7       # (ms)
baseSynDelay_LEC_GC = 0.7       # (ms)
baseSynDelay_LEC_BC = 0.7       # (ms)
baseSynDelay_GC_BC = 0.7        # (ms)
baseSynDelay_MF = 0.7           # (ms)
baseSynDelay_BC_GC = 0.7        # (ms)

###################
# Synapse Weights #
###################
weight_dict = {}
weight_dict['granulecell'] = {}
weight_dict['basketcell'] = {}
weight_dict['ca3pyramidalcell_Bakerb'] = {}
weight_dict['ca3pyramidalcell_Bakerc'] = {}
weight_dict['ca3pyramidalcell_Bakerd'] = {}
weight_dict['ca3interneuron'] = {}

weight_dict['granulecell']['MEC'] = 4.335266e-04
weight_dict['granulecell']['LEC'] = 2.262146e-04
weight_dict['granulecell']['BC'] = 1.074458e-03

weight_dict['basketcell']['MEC'] = 1.438995e-04
weight_dict['basketcell']['LEC'] = 1.438995e-04
weight_dict['basketcell']['GC'] = 3.298706e-04

weight_dict['ca3pyramidalcell_Bakerb']['MEC'] = 2.993835e-04
weight_dict['ca3pyramidalcell_Bakerb']['LEC'] = 3.115784e-04
weight_dict['ca3pyramidalcell_Bakerb']['GC'] = 1.049365e-03
weight_dict['ca3pyramidalcell_Bakerb']['GCbasal'] = 1.171313e-03
weight_dict['ca3pyramidalcell_Bakerb']['rad'] = 3.603577e-04
weight_dict['ca3pyramidalcell_Bakerb']['ori'] = 3.603577e-04
weight_dict['ca3pyramidalcell_Bakerb']['BC'] = 1.782459e-03

weight_dict['ca3pyramidalcell_Bakerc']['MEC'] = 3.115784e-04
weight_dict['ca3pyramidalcell_Bakerc']['LEC'] = 3.176758e-04
weight_dict['ca3pyramidalcell_Bakerc']['GC'] = 1.024976e-03
weight_dict['ca3pyramidalcell_Bakerc']['GCbasal'] = 1.146924e-03
weight_dict['ca3pyramidalcell_Bakerc']['rad'] = 3.847473e-04
weight_dict['ca3pyramidalcell_Bakerc']['ori'] = 3.786499e-04
weight_dict['ca3pyramidalcell_Bakerc']['BC'] = 1.855700e-03

weight_dict['ca3pyramidalcell_Bakerd']['MEC'] = 2.932861e-04
weight_dict['ca3pyramidalcell_Bakerd']['LEC'] = 3.054810e-04
weight_dict['ca3pyramidalcell_Bakerd']['GC'] = 1.390820e-03
weight_dict['ca3pyramidalcell_Bakerd']['GCbasal'] = 1.561548e-03
weight_dict['ca3pyramidalcell_Bakerd']['rad'] = 3.298706e-04
weight_dict['ca3pyramidalcell_Bakerd']['ori'] = 3.359680e-04
weight_dict['ca3pyramidalcell_Bakerd']['BC'] = 1.660390e-03

weight_dict['ca3interneuron']['LEC'] = 1.438995e-04
weight_dict['ca3interneuron']['MEC'] = 1.438995e-04
weight_dict['ca3interneuron']['rad'] = 4.152344e-04
weight_dict['ca3interneuron']['ori'] = 4.152344e-04
weight_dict['ca3interneuron']['GC'] = 2.323120e-04
weight_dict['ca3interneuron']['GCbasal'] = 1.835327e-04

###################
# AMPA parameters #
###################
tau_dict = {}
tau_dict['granulecell'] = {}
tau_dict['basketcell'] = {}
tau_dict['ca3pyramidalcell_Bakerb'] = {}
tau_dict['ca3pyramidalcell_Bakerc'] = {}
tau_dict['ca3pyramidalcell_Bakerd'] = {}
tau_dict['ca3interneuron'] = {}

tau_dict['granulecell']['MEC'] = 1.777344
tau_dict['granulecell']['LEC'] = 4.109375
tau_dict['granulecell']['BC'] = 6.441406

tau_dict['basketcell']['MEC'] = 16.423438
tau_dict['basketcell']['LEC'] = 16.423438
tau_dict['basketcell']['GC'] = 1.661328

tau_dict['ca3pyramidalcell_Bakerb']['MEC'] = 12.660156
tau_dict['ca3pyramidalcell_Bakerb']['LEC'] = 12.660156
tau_dict['ca3pyramidalcell_Bakerb']['GC'] = 144.031250
tau_dict['ca3pyramidalcell_Bakerb']['GCbasal'] = 144.031250
tau_dict['ca3pyramidalcell_Bakerb']['rad'] = 8.384766
tau_dict['ca3pyramidalcell_Bakerb']['ori'] = 8.773438
tau_dict['ca3pyramidalcell_Bakerb']['BC'] = 8.384766

tau_dict['ca3pyramidalcell_Bakerc']['MEC'] = 11.882812
tau_dict['ca3pyramidalcell_Bakerc']['LEC'] = 11.882812
tau_dict['ca3pyramidalcell_Bakerc']['GC'] = 128.484375
tau_dict['ca3pyramidalcell_Bakerc']['GCbasal'] = 131.593750
tau_dict['ca3pyramidalcell_Bakerc']['rad'] = 7.607422
tau_dict['ca3pyramidalcell_Bakerc']['ori'] = 7.996094
tau_dict['ca3pyramidalcell_Bakerc']['BC'] = 7.413086

tau_dict['ca3pyramidalcell_Bakerd']['MEC'] = 16.546875
tau_dict['ca3pyramidalcell_Bakerd']['LEC'] = 16.546875
tau_dict['ca3pyramidalcell_Bakerd']['GC'] = 122.265625
tau_dict['ca3pyramidalcell_Bakerd']['GCbasal'] = 122.265625
tau_dict['ca3pyramidalcell_Bakerd']['rad'] = 11.494141
tau_dict['ca3pyramidalcell_Bakerd']['ori'] = 11.882812
tau_dict['ca3pyramidalcell_Bakerd']['BC'] = 13.437500

tau_dict['ca3interneuron']['MEC'] = 16.423438
tau_dict['ca3interneuron']['LEC'] = 16.423438
tau_dict['ca3interneuron']['rad'] = 1.661328
tau_dict['ca3interneuron']['ori'] = 1.661328
tau_dict['ca3interneuron']['GC'] = 13.335938
tau_dict['ca3interneuron']['GCbasal'] = 11.020312

# End of file
