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

<<<<<<< HEAD

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir
    pwd=os.getcwd()
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
          fout.write("cd /user/mgul/Higgs_tottbar/anlyzer808/CMSSW_8_0_11/src\n")
          fout.write("eval `scram runtime -sh`\n")
          fout.write("pwd\n")
          fout.write("export X509_USER_PROXY=/user/mgul/Higgs_tottbar/anlyzer808/CMSSW_8_0_11/src/Tupel/Tupel\n")
          fout.write("python %s/scripts/testing_checkProductionIntegrity.py -i /store/user/mgul/tuple_8011_new1/1c09f31/%s -o /store/user/mgul/test1 --nocheck 0 \n"% (pwd,dsetname))
        os.system("chmod 755 job_%s.sh"% dsetname)
   
   ###### sends bjobs ######
        os.system("qsub -q localgrid@cream02 -o script.stdout -e script.stderr job_%s.sh"% dsetname)
=======
#    Popen([eos_cmd, ' -b fuse mount', 'eos'],stdout=PIPE).communicate()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir
#    Popen([eos_cmd, 'mkdir', '/eos/cms/'+opt.outDir],stdout=PIPE).communicate()
#    Popen(['mkdir', '/pnfs/iihe/cms'+opt.outDir],stdout=PIPE).communicate()

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    for dset in dset_list:
        dsetname=dset.split('/')[-1]

        pub_list=getEOSlslist(directory=dset,prepend='')
        for pubDir in pub_list:

            if not 'crab' in pubDir:
                print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
                continue
            pub=pubDir.split('/crab_')[-1]

#            localMerge='python /user/mgul/Higgs_tottbar/anlyzer808/CMSSW_8_0_11/src/Tupel/Tupel/scripts/checkProductionIntegrity.py -i %s -o %s --nocheck --only %s'%(opt.inDir,opt.outDir,pub)
            localMerge='python /user/mgul/Higgs_tottbar/anlyzer808/CMSSW_8_0_11/src/Tupel/Tupel/scripts/checkProductionIntegrity.py -i %s -o %s --nocheck --only %s'%(opt.inDir,opt.outDir,pub)
#            cmd='bsub -q %s %s/src/Tupel/Tupel/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localMerge)
            cmd='qsub -q localgrid@cream02 -o script.stdout -e script.stderr /user/mgul/Higgs_tottbar/anlyzer808/CMSSW_8_0_11/src/Tupel/Tupel/scripts/wrapLocalAnalysisRun.sh %s %s'% (cmsswBase,localMerge)
            os.system(cmd)

#    Popen([eos_cmd, ' -b fuse umount', 'eos'],stdout=PIPE).communicate()
>>>>>>> 22e0d1f5b734148c99f2f8f885198e7750e4829b

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

