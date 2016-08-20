import os
import ROOT
import sys
import optparse
from Tupel.Tupel.storeTools import getEOSlslist
from subprocess import Popen, PIPE

"""
steer the script
"""
def main():

#    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',               default=None,   type='string')
    parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory with files',              default=None,   type='string')
    parser.add_option('-c', '--cleanup',     dest='cleanup',     help='removes original crab directory',          default =False, action='store_true')
    parser.add_option(      '--nocheck',     dest='nocheck',     help='do not prompt user',                       default=False,  action='store_true')
    parser.add_option(      '--only',        dest='only',        help='only this tag',                            default=None,   type='string')
    (opt, args) = parser.parse_args()

#    Popen([eos_cmd, ' -b fuse mount', 'eos'],stdout=PIPE).communicate()

    #prepare output directory
    if opt.outDir is None: opt.outDir=opt.inDir
#    Popen(['mkdir', '/pnfs/iihe/cms'+opt.outDir],stdout=PIPE).communicate()

    dset_list=getEOSlslist(directory=opt.inDir,prepend='')
    os.system("mkdir temp_merged_files/")
    pub_list=getEOSlslist(directory=opt.inDir,prepend='')
    for pubDir in pub_list:
        dsetname=pubDir.split('/')[-1]

        if not 'crab' in pubDir:
            print 'Ambiguity found @ <publication-name> for <primary-dataset>=%s , bailing out'%dsetname
            continue
            #select if it doesn't match the required selection
        pub=pubDir.split('/crab_')[-1]
        if opt.only:
            if pub!=opt.only: 
                continue

            #check if it's an extension
        pubExt=None
        try:
            extSplit=pub.split('_ext')
            pubExt='ext'+extSplit[1]
            pub=extSplit[0]
            print 'Extension will be postfixed with ',pubExt
        except:
            print 'Core sample (no extension)'
              
        time_list=getEOSlslist(directory=pubDir,prepend='')
        if len(time_list)!=1:
            print 'Ambiguity found @ <time-stamp> for <primary-dataset>=%s , bailing out'%dsetname
            continue
        time_stamp=time_list[0].split('/')[-1]

        out_list=[]
        count_list=getEOSlslist(directory=time_list[0],prepend='')
        for count in count_list: out_list += getEOSlslist(directory=count,prepend='')
        file_list=[x for x in out_list if '.root' in x]

        newDir='%s/%s' % (opt.outDir,pub)        
        print '<primary-dataset>=%s <publication-name>=crab_%s <time-stamp>=%s has %d files' % (dsetname,pub,time_stamp,len(file_list) )
        if not opt.nocheck:
            choice = raw_input('Will move to %s current output directory. [y/n] ?' % newDir ).lower()
            if not 'y' in choice : continue
            
#           Popen([eos_cmd, 'mkdir', '/eos/cms/'+newDir],stdout=PIPE).communicate()
        Popen(['srmmkdir', 'srm://maite.iihe.ac.be:8443/pnfs/iihe/cms'+newDir],stdout=PIPE).communicate()
    
        moveIndividualFiles=True
        if len(file_list)>0:
               #subgroupMerge = int( raw_input('This set has %d files. Merge into groups? (enter 0 if no merging)' % len(file_list)) )
            subgroupMerge=100 if 'Data' in dsetname else 10 
            if 'TTJets' or 'WJets' or 'DYJets' in pub : subgroupMerge=50
#                if '/store/cmst3/group/hintt' in opt.inDir: 
#                    subgroupMerge=10 if '/data' in dsetname else 3

            if subgroupMerge>0:
                moveIndividualFiles=False

                splitFunc = lambda A, n=subgroupMerge: [A[i:i+n] for i in range(0, len(A), n)]
                split_file_lists = splitFunc( file_list )
                
                for ilist in xrange(0,len(split_file_lists)):
                    if pubExt:
                        mergedFileName='temp_merged_files/MergednTuple_%s_%d_%s.root '%(pub,ilist,pubExt)
                        mergedFileNameonly='MergednTuple_%s_%d_%s.root '%(pub,ilist,pubExt)
                    else:
                        mergedFileName='temp_merged_files/MergednTuple_%s_%d.root '%(pub,ilist)
                        mergedFileNameonly='MergednTuple_%s_%d.root '%(pub,ilist)
                    toAdd='%s ' % mergedFileName
                    for f in split_file_lists[ilist]:                            
                        toAdd += '/pnfs/iihe/cms%s '%f 

                    finalOutput='/pnfs/iihe/cms%s/%s'%(newDir,mergedFileName.replace('temp_merged_files/',''))
                    fIn=ROOT.TFile.Open(finalOutput)
                    try:
                        if fIn or not fIn.IsZombie():
                            print '%s already in pnfs, skipping'%finalOutput
                            fIn.Close()
                            continue
                    except:
                        pass
                        
                    os.system('hadd -f -k %s'%toAdd)
#                        os.system('cp %s /pnfs/iihe/cms%s/'%(mergedFileName,newDir))
                    os.system('srmcp file:///%s srm://maite.iihe.ac.be:8443/pnfs/iihe/cms%s/%s'%(mergedFileName,newDir,mergedFileNameonly))
                    os.system('rm %s'%mergedFileName)
                        #os.system('xrdcp  -f %s root://eoscms//eos/cms/%s/MergedMiniEvents_%d.root' %(mergedFileName,newDir,ilist))

                #if still needed copy individual files
            if moveIndividualFiles:
                 for f in file_list : 
                        #os.system('xrdcp  -f %s eos/cms/%s/' % (f, newDir) )
                    newF=f
                    if pubExt:
                        newF=f.replace('.root','_%s.root'%pubExt)

#                        os.system('cp %s /pnfs/iihe/cms%s/%s' % (f, newDir,newF) )
                    os.system('srmcp %s file:///%s srm://maite.iihe.ac.be:8443/pnfs/iihe/cms%s/%s' % (f, newDir,newF) )
                    print 'this is file: %s %s %s'% (f,newDir,newF)

        if not opt.nocheck and opt.cleanup : 
            choice = raw_input('Will remove output directory. [y/n] ?').lower()
            if 'y' in choice: 
                Popen(['rm', '-r /pnfs/iihe/cms'+pubDir],stdout=PIPE).communicate()

        print 'Crab outputs may now be found in %s' % newDir

    #Popen([eos_cmd, ' -b fuse umount', 'eos'],stdout=PIPE).communicate()
    print '-'*50
    print 'All done. In case errors were found check that the crab output structure is '
    print '<outLFNDirBase>/<primary-dataset>/<publication-name>/<time-stamp>/<counter>/<file-name>'
    print '-'*50
        


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

