import os
# the directory structure should be prepared before the run
# the codes should be placed in ./Codes with write permission
#os.putenv ("PYTHONPATH","/opt/nicadd/contrib/piot/lib/python2.4/site-packages/")
#os.environ['PYTHONPATH'] = '/opt/nicadd/contrib/piot/lib/python2.4/site-packages/'

# you have to do a 
#    setenv PYTHONPATH /opt/nicadd/contrib/piot/lib/python2.4/site-packages/
#    export PYTHONPATH=/opt/nicadd/contrib/piot/lib/python2.4/site-packages/
# before running this scripts

import sys


from S2E_tools import *

#configuration
runtype='new'; #'restart','new'
fig_dir='./Figures/';
report_file=fig_dir+'report.html';

inputfile=sys.argv[1]

#inputfile='inputCRP.txt_ASTRA2BOTH'
#inputfile='inputCRP_testUnifiedPlotter.txt'
#inputfile='inputCRP_testSequence_2.txt'
#inputfile='inputCRP_ELEGANT2FMT1.txt'
#inputfile='inputCRP_astra_downsample_fmt1.txt'
#inputfile='inputCRP.txt_plotTest'
#inputfile='input.txt_CSRtrack_downSample_match_run'
#inputfile='input.txt_FULL_ImpactZ_CSRtrack_beamline'
#inputfile='input.txt_LmatchingTest'
#inputfile='input.txt_dechirp_ForElegant'
#inputfile='input.txt_SliceAnalysisTest'
#inputfile='input_test_CSR.txt'
#inputfile='input_ASTRA2ELEGANT.txt'
#inputfile='input_ImpactBC1_Linearize_Chirped.txt'
#inputfile='input_CSRtrack_Linearize_Chirped_VariableSigma.txt'
#inputfile='input_Impact2SDDS.txt'
#inputfile='input_ImpactTesting.txt'
#inputfile='input_ImpactBC1_Chirped_NotLinear.txt'
#inputfile='input_ImpactConvergenceTesting.txt'
#inputfile='input_CSRtrack_ConvergenceCompare.txt'
#inputfile='input_reCenter_test.txt'
#inputfile='input_CSRtrack_Curved_Chirped.txt'
#inputfile='input_ImpactBC1_320testing.txt'
#inputfile='input_ImpactBC1_410testing.txt'
#inputfile='input_ImpactBC1_410convergence.txt'
#inputfile='input_fourier_test.txt'
#inputfile='input_CSRtrack.txt'
#inputfile='input_CSRtrack_Linearize.txt'
#inputfile='input_CSRtrack_Linearize_Chirped.txt'
#inputfile='input_CSRtrack_Linearize_Chirped_PtoP.txt'
#inputfile='input_CSRtrack_1nC.txt'


#inputfile='input_CSRtrack_250pC.txt'


#inputfile='input_CSRtrack_020pC.txt'


#inputfile='input.txt_FULL_IMPACTZ_rematchingDCM1'

#inputfile='input.txt_FMT3testing'
#inputfile='inputCRP.txt_IMPACTZtest'
#inputfile='inputCRP.txt_ASTRA2FMT1'
#inputfile='inputCRP_justCSR.txt'
#inputfile='inputCRP_astra2fmt1.txt'
#inputfile='inputCRP.txt_ASTRA2FMT1'
outputfile='burblah.txt'

#print "starting conversion"

# main code
if runtype=='new':
    html=open(report_file,'w');
    html.write('<HTML><BODY BGCOLOR="white">\n')
    start,stop,indfig=0,1000,0
if runtype=='restart':
    html=open(report_file,'a+');
    start,stop,indfig=FindRestartParameters(html)

actions=ReadTask(inputfile);
WriteTask(outputfile,actions)
indfig=RunTask(actions,fig_dir,html,start=start,stop=stop,indfig=indfig)
html.write('</BODY></HTML>\n'); html.close()

