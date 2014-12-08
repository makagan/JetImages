#!/usr/bin/env python
import time
import sys
import math
import sys
import os

#####
samples = {
    'h' : ('WHbb_up_to_2j', os.getenv('HomeDir', '.') + '/data/WHbb_up_to_2j.dat'),
    'g' : ('Wbb_up_to_2j' , os.getenv('HomeDir', '.') + '/data/Wbb_up_to_2j.dat' ),
    }

def runArgs(args):
    args = " ".join(args)
    print args
    result = subprocess.call(args, shell=True)
    if result != 0:
        raise ValueError('Returned %d for %s' % (result, args))
    return result

def modifyCards(configName, nevents, seed, jetminpt, jetminm):
    import re
    r_run = [ (r'(.*)[0-9]* *= nevents(.*)', '  %d = nevents' % nevents              ),
              (r'(.*)[0-9]* *= iseed(.*)'  , '  %d = iseed' % seed                   ),
              (r'(.*)[0-9]* *= mmbb(.*)'   , '  %d = mmbb' % max(jetminm - 20, 0)   ),
              (r'(.*)[0-9]* *= ptllmin(.*)', '  %d = ptllmin' % max(jetminpt - 30,0) ),
              ]

    run_card = open('%s/Cards/run_card.dat' % samples[configName][0]).read()
    for pattern, replacement in r_run:
        
        print "="*5
        print pattern
        print replacement
        sys.stdout.flush()
        
        m = re.search(pattern, run_card)
        if m:
            run_card = run_card.replace(m.group(0), replacement + '    ! local modification')
        else:
            raise ValueError('No patter matching %s' % pattern)

    open('%s/Cards/run_card.dat' % samples[configName][0], 'w').write(run_card)

    print "EAS: start run_card"
    print run_card
    print "EAS: end run_card"
    sys.stdout.flush()

def runMadgraph(configName, nevents, seed, jetminpt, jetmaxpt, jetminm, jetmaxm):
    print "runMadgraph %s:" % configName
    result = 0

    #result += runArgs(args)

    args = [os.getenv('MG5') + '/bin/mg5', '-f', samples[configName][1]]
    result += runArgs(args)

    modifyCards(configName, nevents, seed, jetminpt, jetminm)
    
    args = ['%s/bin/generate_events' % samples[configName][0], '-f']
    result += runArgs(args)

    args = ['gunzip', '%s/Events/run_01/unweighted_events.lhe.gz' % samples[configName][0]]
    result += runArgs(args)

    print "@@ Start LHE File"
    print open('%s/Events/run_01/unweighted_events.lhe' % samples[configName][0]).read()
    print "@@ End LHE File"

    return result

def runPythiaJets(seed, metype
                  , savemaxnjets, nevents
                  , cellv, trimjet
                  , jetminpt, jetmaxpt
                  , jetminm, jetmaxm
                  , rad, mu, smear
                  , wtagger, input_file
                 ):

    print "runPythiaJets:"

    args = ["pythiaJets.exe"]
    args += ["--seed=%s" % str(seed)                 , "--matching"                           ,
             "--metype=%s" % str(metype)             , "--LHEFile=%s" % input_file            ,
             "--savemaxnjets=%s" % str(savemaxnjets) , "--nevents=%s" % str(nevents)          ,
             "--cellv=%s" % str(cellv)               , "--trimjet=%s" % str(trimjet)          ,
             "--jetminpt=%s" % str(jetminpt)         , "--jetmaxpt=%s" % str(jetmaxpt)        ,
             "--jetminm=%s" % str(jetminm)           , "--jetmaxm=%s" % str(jetmaxm)          ,
             "--rad=%s" % str(rad)                   , "--mu=%s" % str(mu)                    ,
             ]
    
    if smear:
        args.append("--smear")
    if wtagger:
        args.append("--wtagger")

    result = runArgs(args)
    
    return result

def getSeed(optSeed):
    if optSeed > -1: return optSeed
    
    from time import time
    import os
    timeSeed = int(time() * pow(10,9))
    
    return abs(((timeSeed*181)*((os.getpid()-83)*359))%104729)

if __name__ == '__main__':
    import sys
    import argparse
    import subprocess

    print "Called with:"
    print " ".join(sys.argv)
    print "HomeDir: ", os.getenv("HomeDir", "-1")
    print ""
    
    
    parser = argparse.ArgumentParser(description='Run MadGraph + pythiaJets.')
    parser.add_argument('--seed'      , default=-1       , type=int                                         , help='The random number generator seed (-1 finds a random seed)')
    parser.add_argument('--savemaxnjets', default=-1     , type=int                                         , help='Max number of jets to save.  Default, -1, saves all jets (which satisfy other selectons)')
    parser.add_argument('--nevents'   , default=100      , type=int                                         , help='The generator number of events')
    parser.add_argument('--metype'    , default='h'      , type=str   , choices = ['h', 'g']                , help='Matrix element type')
    parser.add_argument('--lepcharge' , default=1        , type=int   , choices = [-1, 1, 0]                , help='')
    parser.add_argument('--meminpt'   , default=0        , type=int                                         , help='')
    parser.add_argument('--memaxpt'   , default=20000    , type=int                                         , help='')
    parser.add_argument('--jetminpt'  , default=0        , type=int                                         , help='Jet PTMin (GeV)')
    parser.add_argument('--jetmaxpt'  , default=20000    , type=int                                         , help='Jet PTMax (GeV)')
    parser.add_argument('--jetminm'   , default=0        , type=int                                         , help='Jet Min Mass (GeV)')
    parser.add_argument('--jetmaxm'   , default=20000    , type=int                                         , help='Jet Max Mass (GeV)')
    parser.add_argument('--cellv'     , default=1        , type=int   , choices = [0, 1]                    , help='Control how many and which cells to output.  0 shows no cells, 1 (default) shows only those that survive trimming and 2 shows all cells in a jet.')
    parser.add_argument('--trimjet'   , default=1        , type=int                                         , help='To run trimming or not.  Any value other than 1 (default) disables trimming.')
    parser.add_argument('--rad'       , default=0.4      , type=float                                       , help='Jet Finder Radius')
    parser.add_argument('--mu'        , default=0.0      , type=float                                       , help='Average Int./Cross')
    parser.add_argument('--smear'     , default=False    , action='store_true'                              , help='Turn on radial energy smearing')
    parser.add_argument('--wtagger'   , default=False    , action='store_true'                              , help='Turn on MVA Tagger output')
    
    args = parser.parse_args()

    args.seed = getSeed(args.seed)
    
    sys.stdout.flush()
    print ""
    print "-+"*50
    print "runMadgraph"
    print runMadgraph(args.metype, args.nevents, args.seed
                      , args.jetminpt, args.jetmaxpt
                      , args.jetminm, args.jetmaxm
                      )

    sys.stdout.flush()
    print ""
    print "-+"*50
    print "runPythiaJets"
    print runPythiaJets(seed=args.seed, metype=args.metype
                       , savemaxnjets=args.savemaxnjets, nevents=args.nevents
                       , cellv=args.cellv, trimjet=args.trimjet
                       , jetminpt=args.jetminpt, jetmaxpt=args.jetmaxpt
                       , jetminm=args.jetminm, jetmaxm=args.jetmaxm
                       , rad=args.rad, mu=args.mu, smear=args.smear
                       , wtagger=args.wtagger, input_file='%s/Events/run_01/unweighted_events.lhe' % samples[args.metype][0]
                       )
    sys.stdout.flush()
