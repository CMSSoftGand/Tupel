import os
import sys
import optparse
from Tupel.Tupel.storeTools import getEOSlslist
from subprocess import Popen, PIPE

"""
steer the script
"""
def main():

    cmsswBase=os.environ['CMSSW_BASE']
#    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',      dest='inDir',       help='input directory with files',               default=None,   type='string')
    parser.add_option('-o', '--outDir',     dest='outDir',      help='output directory with files',              default=None,   type='string')
    parser.add_option('-q', '--queue',      dest='queue',       help='batch queue',                              default='2nd',  type='string')
    (opt, args) = parser.parse_args()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir
    pwd=os.getcwd()
    os.system("cp $X509_USER_PROXY %s" % pwd)
    if os.path.exists("tmp_combine"):
        os.system("rm -rf tmp_combine/")
    os.system("mkdir tmp_combine/")
    os.chdir("tmp_combine")
    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    for dset in dset_list:
        dsetname=dset.split('/')[-1]
        print 'dsetname<>>>>>>>>>>>>>>>> %s'% dsetname
        with open('job_%s.sh'% dsetname, 'w') as fout:
          fout.write("#!/bin/sh\n")
          fout.write("pwd=$PWD\n")
          fout.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
          fout.write("cd %s%s\n"%(pwd,"/../../"))
          fout.write("eval `scram runtime -sh`\n")
          fout.write("pwd\n")
          fout.write("export X509_USER_PROXY=%s\n"% pwd)
          fout.write("python %s/scripts/checkProductionIntegrity.py -i %s%s -o %s --nocheck 0 \n"% (pwd,opt.inDir,dsetname,opt.outDir))
        os.system("chmod 755 job_%s.sh"% dsetname)
   
   ###### sends bjobs ######
        os.system("qsub -q localgrid@cream02 -o script.stdout -e script.stderr job_%s.sh"% dsetname)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

