import os
# the directory structure should be prepared before the run
# the codes should be placed in ./Codes with write permission
#os.putenv ("PYTHONPATH","/opt/nicadd/contrib/piot/lib/python2.4/site-packages/")
#os.environ['PYTHONPATH'] = '/opt/nicadd/contrib/piot/lib/python2.4/site-packages/'

# you have to do a 
#    setenv PYTHONPATH /opt/nicadd/contrib/piot/lib/python2.4/site-packages/
# before running this scripts

from S2E_tools import *

#configuration
runtype='new'; #'restart','new'
fig_dir='./Figures/';
report_file=fig_dir+'report.html';
inputfile='inputAPLCTRslits.txt'
outputfile='inputAPLCTRslits.txt'

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

