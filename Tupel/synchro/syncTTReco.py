#!/usr/bin/env python

"""Compares ROOT files with information on tt reconstruction."""


import argparse
from collections import defaultdict
from functools import total_ordering
import os
from uuid import uuid4

import ROOT


@total_ordering
class EventID(object):
    """A class to describe event ID."""
    
    def __init__(self, runNumber=None, lumiSectionNumber=None, eventNumber=None):
        """Constructor with a complete initialization."""
        
        if runNumber is not None and lumiSectionNumber is not None and eventNumber is not None:
            self.run = int(runNumber)
            self.lumiSection = int(lumiSectionNumber)
            self.event = int(eventNumber)
        else:
            self.run = None
            self.lumiSection = None
            self.event = None
    
    
    def is_valid(self):
        """Check if the ID has been fully initialized."""
        
        return (self.run is not None and self.lumiSection is not None and self.event is not None)
    
    
    def __eq__(self, other):
        return self.run == other.run and self.lumiSection == other.lumiSection and \
            self.event == other.event
    
    
    def __lt__(self, other):
        """Comparison operator to define ordering."""
        
        if not self.is_valid() or not other.is_valid():
            raise RuntimeError('Invalid event ID')
        
        if self.run != other.run:
            return self.run < other.run
        elif self.lumiSection != other.lumiSection:
            return self.lumiSection < other.lumiSection
        else:
            return self.event < other.event
    
    
    def __str__(self):
        """String representation used in printing."""
        
        return '{}:{}:{}'.format(self.run, self.lumiSection, self.event)
    
    
    def __hash__(self):
        """Hash method needed to store objects in a set."""
        
        if self.is_valid():
            return hash((self.run, self.lumiSection, self.event))
        else:
            return 0



class Properties:
    """Aggregates properties of reconstructed events."""
    
    def __init__(self):
        self.massTT = None



if __name__ == '__main__':
    
    # Parse arguments
    argParser = argparse.ArgumentParser(
        epilog=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    argParser.add_argument(
        'input1', metavar='input1.root',
        help='Input ROOT file'
    )
    argParser.add_argument(
        'input2', metavar='input2.root',
        help='Input ROOT file'
    )
    argParser.add_argument(
        '-o', '--fig-dir', metavar='fig/', default='fig/', dest='figDir',
        help='Directory to store output'
    )
    args = argParser.parse_args()
    
    if not args.figDir.endswith('/'):
        args.figDir += '/'
    
    
    # Set up ROOT style
    ROOT.gStyle.SetHistMinimumZero(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetStripDecimals(False)
    ROOT.TGaxis.SetMaxDigits(3)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0.03, 'y')
    
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleFontSize(0.04)
    ROOT.gStyle.SetTitleFont(42, 'XYZ')
    ROOT.gStyle.SetTitleXOffset(1.0)
    ROOT.gStyle.SetTitleYOffset(1.1)
    ROOT.gStyle.SetTitleSize(0.045, 'XYZ')
    ROOT.gStyle.SetLabelFont(42, 'XYZ')
    ROOT.gStyle.SetLabelOffset(0.007, 'XYZ')
    ROOT.gStyle.SetLabelSize(0.04, 'XYZ')
    ROOT.gStyle.SetNdivisions(508, 'XYZ')
    
    ROOT.gROOT.SetBatch(True)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    
    
    # First compare lists of selected events and those in which the tt
    # reconstruction has succeeded
    events = [set(), set()]
    eventsReco = [set(), set()]
    
    for inputFileName, eventsCurFile, eventsRecoCurFile in [
        (args.input1, events[0], eventsReco[0]),
        (args.input2, events[1], eventsReco[1])
    ]:
        inputFile = ROOT.TFile(inputFileName)
        tree = inputFile.Get('sync')
        
        for entry in tree:
            eventID = EventID(entry.Run, entry.LumiSection, entry.Event)
            eventsCurFile.add(eventID)
            
            if entry.RecoSuccess:
                eventsRecoCurFile.add(eventID)
        
        inputFile.Close()
    
    
    print('Selected events:')
    
    print('  In the first file but not in the second:')
    diff = events[0] - events[1]
    
    if diff:
        for event in diff:
            print(' ' * 4 + str(event))
    else:
        print('    (none)')
    
    print('\n  In the second file but not in the first:')
    diff = events[1] - events[0]
    
    if diff:
        for event in diff:
            print(' ' * 4 + str(event))
    else:
        print('    (none)')
    
    print('\n  There are {} events in the overlap'.format(len(events[0] & events[1])))
    
    
    print('\n\nEvents with successful tt reconstruction:')
    
    print('  In the first file but not in the second:')
    diff = eventsReco[0] - eventsReco[1]
    
    if diff:
        for event in diff:
            print(' ' * 4 + str(event))
    else:
        print('    (none)')
    
    print('\n  In the second file but not in the first:')
    diff = eventsReco[1] - eventsReco[0]
    
    if diff:
        for event in diff:
            print(' ' * 4 + str(event))
    else:
        print('    (none)')
    
    eventsOverlap = eventsReco[0] & eventsReco[1]
    print('\n  There are {} events in the overlap'.format(len(eventsOverlap)))
    
    
    # Read properties of overlapping events that have been reconstructed
    # successfully
    properties = defaultdict(list)
    
    for inputFileName in [args.input1, args.input2]:
        
        inputFile = ROOT.TFile(inputFileName)
        tree = inputFile.Get('sync')
        
        for entry in tree:
            eventID = EventID(entry.Run, entry.LumiSection, entry.Event)
            
            if eventID not in eventsOverlap:
                continue
            
            p = Properties()
            p.massTT = float(entry.MassTT)
            
            properties[eventID].append(p)
        
        inputFile.Close()
    
    
    # Check for differences in m(tt)
    print('\n\nReconstructed events with differences in m(tt):')
    mttDiffFound = False
    
    for eventID, p in properties.items():
        if abs(p[0].massTT / p[1].massTT - 1) > 1e-3:
            print('  {:30}  ({:.3f} vs {:.3f})'.format(str(eventID), p[0].massTT, p[1].massTT))
    
    
    # Prepare for drawing
    canvas = ROOT.TCanvas('canvas', '', 1200, 900)
    canvas.SetTicks()
    
    if not os.path.exists(args.figDir):
        os.makedirs(args.figDir)
    
    
    # Draw correlations between reconstructed properties
    for varName, memberName, nBins, plotRange in [
        ('MassTT', 'massTT', 90, [300., 1100.])
    ]:
        
        hist = ROOT.TH2D(
            str(uuid4()), ';{varName}, input #1;{varName}, input #2'.format(varName=varName),
            nBins, plotRange[0], plotRange[1], nBins, plotRange[0], plotRange[1]
        )
        
        for p in properties.values():
            hist.Fill(getattr(p[0], memberName), getattr(p[1], memberName))
        
        hist.Draw('colz')
        canvas.Print(args.figDir + varName + '.pdf')
