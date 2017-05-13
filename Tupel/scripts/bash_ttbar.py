#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import optparse
import ROOT
import pickle
import time
from datetime import datetime
from Tupel.Tupel.storeTools import getEOSlslist
def main():

    #configuration                                                                                                               
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--in',  dest='input', help='input directory', default=None,  type='string')
    (opt, args) = parser.parse_args()

    queue = "1nh"
    from subprocess import Popen, PIPE
    path = os.getcwd()
    os.system("rm -r root_files")
    os.system("mkdir root_files")
    os.system("rm -r tmp")
    os.system("mkdir tmp")
    os.chdir("tmp/")
    bigjob = open('bigjobSub.txt', 'w')
    directory=opt.input
    dir_list=getEOSlslist(directory=opt.input,prepend='')
    for dir in dir_list:
      dsetname=dir.split('/')[-1]
      os.system("mkdir "+path+"/root_files/"+dsetname)
      os.system("mkdir "+dsetname)
      print "current time is: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S')
      print 'looking into: '+directory+'/'+dsetname+'...'
      prepend='dcap://maite.iihe.ac.be/pnfs/iihe/cms'
      print 'this is jjpath: %s '% path
      print 'do not worry about folder creation:'
      data = Popen(['ls', '/pnfs/iihe/cms/'+dir],stdout=PIPE)
      out,err = data.communicate()
      full_list = []
      for line in out.split('\n'):
        if len(line.split()) == 0 or line == "failed" or line == "log": continue
        full_list.append(prepend + directory + '/' + dsetname + '/' + line)
      input_list=[]
      input_list=full_list
      lent=len(input_list)
#use lent for testing purpose
#      lent=3
      x=1
      for ifile in xrange(0,lent):
        inF=input_list[ifile]
        outF="output_"+dsetname+"_"+str(x)
#        os.system("mkdir -p tmp/"+dsetname+"/"+str(x))
        os.system("mkdir -p "+dsetname+"/"+str(x))
#        os.chdir("tmp/"+dsetname+"/"+str(x))
        os.chdir(dsetname+"/"+str(x))
        with open(dsetname+'job.sh', 'w') as fout:
          fout.write("#!/bin/sh\n")
          fout.write("cd "+str(path)+"\n")
          fout.write("source $VO_CMS_SW_DIR/cmsset_default.sh \n")
          fout.write("eval `scram runtime -sh` \n")
          fout.write("root -l -b <<EOF \n")
#         fout.write("gSystem->Load("'"libFWCoreFWLite.so"'") \n")
#         fout.write("FWLiteEnabler::enable() \n")
          fout.write(".x analyzer/simpleReader.cc+ (" '"'+inF+'"' "," '"'+path+"/root_files"+"/"+dsetname+"/"+outF +'"' " ) \n")
#          fout.write(".x gen_plots/simpleReader.cc+ (" '"'+inF+'"' "," '"'+path+"/root_files"+"/"+dsetname+"/"+outF +'"' " ) \n")
          bigjob.write("qsub -q localgrid@cream02 "+dsetname+"/"+str(x)+"/"+dsetname+"job.sh \n")
          fout.write(".q; \n")
          fout.write("EOF\n")
#          fout.write("cp " +outF+".root "+str(path)+"/root_files/"+dsetname+"\n")
#          fout.write("rm "+str(path)+"/"+outF+".root \n")
          fout.write("cd - \n")
#         fout.write("echo 'STOP---------------'\n")
#         fout.write("echo\n")
        os.system("chmod 755 "+dsetname+"job.sh")
   
   ###### sends bjobs ######
#        os.system("qsub -q highmem@cream02 -o script.stdout -e script.stderr "+dsetname+"job.sh")
        #os.system("qsub -q localgrid@cream02 -o script.stdout -e script.stderr "+dsetname+"job.sh")
        #print "job nr " + str(x) + " submitted"
        os.chdir("../..")
        print
        print "your jobs:"
        print "output file : ", outF
        print
        print 'END'
        print
        x = x+1
        print " value of x", x
if __name__ == "__main__":
    sys.exit(main())
