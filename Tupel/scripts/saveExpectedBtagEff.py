import ROOT
from array import array
from Tupel.Tupel.storeTools import *
import optparse
import sys

"""
"""
def saveExpectedBtagEff(opt):

    csvWP='patJetPfAk04BDiscCSVv2_>%f'%opt.csv if opt.HiForest else 'patJetPfAk04BDiscCSVv2_>%f'%opt.csv
    inputDir=opt.input
    
    #open a file
    input_list=getEOSlslist(directory=inputDir)
    data=ROOT.TChain('ak4PFJetAnalyzer/t') if opt.HiForest else ROOT.TChain('tupel/MuonTree') 
    for i in xrange(0,min(50,len(input_list))):
    #for i in xrange(0,len(input_list)):
        data.Add(input_list[i])

    print 'Projecting tagging efficiency from ',data.GetEntries(),' events'
    print 'Working point is',csvWP


    #prepare histograms
    effgrs=[]
    ptBins = [0,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,800,1000]
    preTagH=ROOT.TH1F('preTagH',';p_{T} [GeV];',len(ptBins)-1,array('d',ptBins))
    preTagH.Sumw2()
    tagH=preTagH.Clone('tagH')


    #count number of tagged jets
    flavConds=[('b',"abs(patJetPfAk04PartonFlavour_)==5"),
              ('c',"abs(patJetPfAk04PartonFlavour_)==4"),
              ('udsg','abs(patJetPfAk04PartonFlavour_)!=5 && abs(patJetPfAk04PartonFlavour_)!=4'),
              ('pu','abs(patJetPfAk04PartonFlavour_)!=5 && abs(patJetPfAk04PartonFlavour_)!=4')]
#              ('pu','abs(patJetPfAk04PartonFlavour_)!=5 && abs(patJetPfAk04PartonFlavour_)!=4 && abs(j_g)<0')] #see the j_g
    ptVar='patJetPfAk04PtJERSmear'
    if opt.HiForest:
        flavConds=[('b',"abs(refparton_flavorForB)==5"),
                   ('c',"abs(refparton_flavorForB)==4"),
                   ('udsg',"abs(refparton_flavorForB)!=5 && abs(refparton_flavorForB)!=4")
                   ]
        ptVar='jtpt'
    for flav,cond in flavConds:
        print '\t computing for',flav,cond
        preTagH.Reset('ICE')
        tagH.Reset('ICE')
        data.Draw("%s >> preTagH"%ptVar,cond,'goff')
        data.Draw('%s >> tagH'%ptVar,cond + ' && ' + csvWP,'goff')
        effgrs.append(ROOT.TGraphAsymmErrors())
        effgrs[-1].SetName(flav)
        effgrs[-1].Divide(tagH,preTagH)

    #save efficiency graphs
    fOut=ROOT.TFile.Open(opt.output,'RECREATE')
    for gr in effgrs:
        gr.Write()
    fOut.Close()


"""
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(      '--HiForest',  dest='HiForest', help='use HiForest trees',          default=False, action='store_true')
    parser.add_option('-i', '--in',        dest='input',    help='input directory with files',  default='/store/user/mgul/tupel_8011/f862e06/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_MC13TeV_TTJets_powheg/160812_165322/0000', type='string')
    parser.add_option('-o', '--out',       dest='output',   help='output file',                 default='data/expTageff.root',                                       type='string')
    parser.add_option(      '--csv',       dest='csv',      help='csv cut',                     default=0.800,                                                       type=float)
    (opt, args) = parser.parse_args()

    saveExpectedBtagEff(opt)

"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
