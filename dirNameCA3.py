##################################
# Creating output directory name #
##################################
import sys
import datetime

##################################
# Import connectivity parameters #
##################################
ECtoGC = int(sys.argv[5])
ECtoCA3 = int(sys.argv[6])
ECtoBC = int(sys.argv[7])
GCtoBC = int(sys.argv[8])
BCtoGC = int(sys.argv[9])
GCtoCA3 = int(sys.argv[10])
CA3toCA3 = int(sys.argv[11])
CA3inhibition = int(sys.argv[12])

GC_spike_flag = int(sys.argv[13])
GC_rand_flag = int(sys.argv[14])

if GC_spike_flag == 1:
        ECtoGC = 0

# Setting flags so that certain cell types aren't instantiated if
# they are disconnected from the simulation
if (ECtoGC==1)|(GCtoCA3==1):
	GC_flag = 1
else:
	GC_flag = 0

if (ECtoBC==1)|(GCtoBC==1):
	BC_flag = 1
else:
	BC_flag = 0

if (ECtoCA3==1)|(GCtoCA3==1)|(CA3toCA3==1):
	CA3_flag = 1
else:
	CA3_flag = 0

nmda_flag = int(sys.argv[15])
slice_flag = int(sys.argv[16])

assoc_Factor = 1
volt_flag = 0
if len(sys.argv) > 16:
	GCspike_file = sys.argv[17]
	assoc_Factor = float(sys.argv[18])
	volt_flag = int(sys.argv[19])
	grid_flag = int(sys.argv[20])
	model_flag = int(sys.argv[21])
	MEC_GC = float(sys.argv[22])
	GC_BC = float(sys.argv[23])
	EC_BC = float(sys.argv[24])
	BC_GC = float(sys.argv[25])
        EC_CA3_BC = float(sys.argv[26])
        CA3_BC = float(sys.argv[27])
        GC_CA3_BC = float(sys.argv[28])
        BC_CA3 = float(sys.argv[29])
        EC_CA3 = float(sys.argv[30])
        GC_CA3_PC = float(sys.argv[31])
        CA3_CA3 = float(sys.argv[32])

##################################
# Creating output directory name #
##################################
outputDir_prefix = "/auto/rcf-proj2/tb/geneyu/Big_Model/EC-DG-CA3_2019/{" 
outputDir_suffix = str(datetime.date.today())+"}"
'''
# Flags that prevent the instantiation of CA3 pyramidal cells not in
# the specified CA3 region
if CA3c_flag == 1:
	outputDir_suffix = 'CA3c_' + outputDir_suffix
	CA3_reject = 	range(N_MEC+N_LEC+N_GC+numBasketCells,
				N_MEC+N_LEC+N_GC+numBasketCells+N_CA3a+N_CA3b)

if CA3b_flag == 1:
	outputDir_suffix = 'CA3b_' + outputDir_suffix
	CA3_reject = 	range(N_MEC+N_LEC+N_GC+numBasketCells
				,N_MEC+N_LEC+N_GC+numBasketCells+N_CA3a)
				+range(N_MEC+N_LEC+N_GC+numBasketCells+N_CA3a+N_CA3b,
				N_MEC+N_LEC+N_GC+numBasketCells+N_CA3a+N_CA3b+N_CA3c)

if CA3a_flag == 1:
	outputDir_suffix = 'CA3a_' + outputDir_suffix
	CA3_reject = 	range(N_MEC+N_LEC+N_GC+numBasketCells+N_CA3a,
				N_MEC+N_LEC+N_GC+numBasketCells+N_CA3a+N_CA3b+N_CA3c)
'''
if model_flag == 1:
	outputDir_suffix = 'Single_' + outputDir_suffix

if grid_flag == 1:
	outputDir_suffix = 'gridCell_' + outputDir_suffix

if assoc_Factor > 1:
	outputDir_suffix = 'Assoc-'+str(assoc_Factor)+'_'+outputDir_suffix

if GC_spike_flag == 1:
	outputDir_suffix = GCspike_file + '_' + outputDir_suffix

if CA3inhibition == 1:
	outputDir_suffix = 'i_' + outputDir_suffix

if CA3toCA3 == 1:
    outputDir_suffix = 'CA3-CA3_' + outputDir_suffix

if GCtoCA3 == 1:
    outputDir_suffix = 'GC-CA3_' + outputDir_suffix

if ECtoCA3 == 1:
    outputDir_suffix = 'EC-CA3_' + outputDir_suffix

if ECtoGC == 1:
    outputDir_suffix = 'EC-GC_' + outputDir_suffix

if ECtoBC == 1:
	outputDir_suffix = 'FFI_' + outputDir_suffix

if GCtoBC == 1:
	outputDir_suffix = 'FBI_' + outputDir_suffix

if BCtoGC == 0:
	outputDir_suffix = 'OpenLoop_' + outputDir_suffix

if nmda_flag == 1:
	outputDir_suffix = 'NMDA_' + outputDir_suffix

if slice_flag == 1:
	outputDir_suffix = 'slice_' + outputDir_suffix

outputDir_suffix = outputDir_suffix + 'EC-CA3-' + str(EC_CA3) + '_'
outputDir_suffix = outputDir_suffix + 'GC-CA3-' + str(GC_CA3_PC) + '_'
outputDir_suffix = outputDir_suffix + 'EC-BC-' + str(EC_CA3_BC) + '_'
outputDir_suffix = outputDir_suffix + 'CA3-BC-' + str(CA3_BC) + '_'
outputDir_suffix = outputDir_suffix + 'GC-BC-' + str(GC_CA3_BC) + '_'
outputDir_suffix = outputDir_suffix + 'BC-CA3-' + str(BC_CA3) + '_'
outputDir_suffix = outputDir_suffix + 'CA3-CA3-' + str(CA3_CA3)

outputDir = outputDir_prefix + outputDir_suffix
