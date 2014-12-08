#!/usr/bin/env python
import time
import sys
import math
import sys

#####
def writeBasicInfo(nevents=100, seed=-1):


    text = """################################################## 
# First add MaxKT cuts 
################################################## 
#cd /Herwig/Cuts 
#newdef JetKtCut:MaxKT 7000.0*GeV 
#newdef WBosonKtCut:MaxKT 7000.0*GeV 
#newdef ZBosonKtCut:MaxKT 7000.0*GeV 
#newdef TopKtCut:MaxKT 7000.0*GeV 
#


################################################## 
# Technical parameters for this run 
################################################## 
cd /Herwig/Generators 
set LHCGenerator:NumberOfEvents %d 
set LHCGenerator:RandomNumberGenerator:Seed %d
set LHCGenerator:DebugLevel 1 
set LHCGenerator:PrintEvent 10 
set LHCGenerator:MaxErrors 10000 
set LHCGenerator:EventHandler:StatLevel Full 
#set LHCGenerator:EventHandler:CascadeHandler NULL 
#Turning this off kills UE as well
set /Herwig/Shower/ShowerHandler:MPIHandler NULL 
#set LHCGenerator:EventHandler:DecayHandler NULL 
#set LHCGenerator:EventHandler:HadronizationHandler NULL 
set /Herwig/Analysis/Basics:CheckQuark 0 
set LHCGenerator:EventHandler:LuminosityFunction:Energy 8000.0 
#set /Herwig/Shower/Evolver:IntrinsicPtGaussian 3.85*GeV 
set /Herwig/Shower/Evolver:IntrinsicPtGaussian 2.2*GeV 

""" % (nevents, seed)
    
    return text


####
def writePDFInfo():
    text = """################################################## 
# PDF parameters 
################################################## 
cd /Herwig/Partons 
create ThePEG::LHAPDF LHAPDF ThePEGLHAPDF.so
set LHAPDF:PType Nucleon 
set LHAPDF:PDFName cteq6ll.LHpdf 
set LHAPDF:PDFNumber 10042 
set LHAPDF:RemnantHandler HadronRemnants 
set /Herwig/Particles/p+:PDF LHAPDF 
set /Herwig/Particles/pbar-:PDF LHAPDF 

"""
    return text


###
def writeUEtune():
    text = """##################################################
# Underlying Event Parameters
##################################################

"""

    return text


#####
def writeMatrixElementInfo(metype):
    text = """##################################################
# Matrix Elements
##################################################
cd /Herwig/MatrixElements/ 
"""

    if metype=='w':
        text +="""create Herwig::MEPP2VV MEPP2VV HwMEHadron.so
insert SimpleQCD:MatrixElements[0] MEPP2VV
set MEPP2VV:Process WW \n"""
    elif metype=="v":
        text +="""insert SimpleQCD:MatrixElements[0] MEWJet \n"""
    elif metype=="t":
        text +="""insert SimpleQCD:MatrixElements[0] MEHeavyQuark \n"""
    elif metype=="q":
        text +="""insert SimpleQCD:MatrixElements[0] MEQCD2to2 \n"""
    elif metype=="b":
        text +="""insert SimpleQCD:MatrixElements[0] MEHeavyQuark 
set MEHeavyQuark:QuarkType Bottom \n"""

    text += """\n"""

    return text


#####
def writeKinematicCutInfo(metype,meminpt,memaxpt):
    text = """##################################################
# Kinematic Cuts
##################################################"""

    if metype=='w' or metype=='v':
         text += """
set /Herwig/Cuts/WBosonKtCut:MinKT %d*GeV
#set /Herwig/Cuts/WBosonKtCut:MaxKT %d*GeV
set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV
""" % (meminpt, memaxpt)
    elif metype=="t":
        text += """
set /Herwig/Cuts/TopKtCut:MinKT    %d*GeV
#set /Herwig/Cuts/TopKtCut:MaxKT    %d*GeV
set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV
""" % (meminpt, memaxpt)
    elif metype=="q":
        text += """
set /Herwig/Cuts/JetKtCut:MinKT %d*GeV 
#set /Herwig/Cuts/JetKtCut:MaxKT %d*GeV
""" % (meminpt, memaxpt)

     #file.write("set /Herwig/Cuts/ZBosonKtCut:MinKT "+str(meminpt)+"*GeV \n")
    #file.write("set /Herwig/Cuts/ZBosonKtCut:MaxKT "+str(memaxpt)+"*GeV \n")
        

    text += """\n"""

    return text


#####
def writeParticleDecays(metype,lepcharge):

    Wp_lep = "Off"
    Wp_had = "Off"

    Wm_lep = "Off"
    Wm_had = "Off"

#    Z_lep = "Off"
#    Z_had = "Off"

    t_lep = "Off"
    t_had = "Off"

    tbar_lep = "Off"
    tbar_had = "Off"

    if lepcharge == 0:
        Wp_lep="On"
        Wm_lep="On"
        t_lep="On"
        tbar_lep="On"
    elif lepcharge == 1:
        Wp_lep="On"
        Wm_had="On"
        t_lep="On"
        tbar_had="On"
    elif lepcharge == -1:
        Wp_had="On"
        Wm_lep="On"
        t_had="On"
        tbar_lep="On"
    else:
        Wp_had="On"
        Wm_had="On"
        t_had="On"
        tbar_had="On"


    text="""################################################## 
# Particle Decays 
################################################## 
set /Herwig/Particles/W+/W+->nu_e,e+;:OnOff Off 
set /Herwig/Particles/W+/W+->nu_mu,mu+;:OnOff %s
set /Herwig/Particles/W+/W+->nu_tau,tau+;:OnOff Off 
set /Herwig/Particles/W+/W+->u,dbar;:OnOff %s 
set /Herwig/Particles/W+/W+->sbar,u;:OnOff %s 
set /Herwig/Particles/W+/W+->bbar,c;:OnOff %s 
set /Herwig/Particles/W+/W+->c,dbar;:OnOff %s 
set /Herwig/Particles/W+/W+->c,sbar;:OnOff %s 
""" % (Wp_lep, Wp_had, Wp_had, Wp_had, Wp_had, Wp_had)


    text += """set /Herwig/Particles/W-/W-->nu_ebar,e-;:OnOff Off 
set /Herwig/Particles/W-/W-->nu_mubar,mu-;:OnOff %s 
set /Herwig/Particles/W-/W-->nu_taubar,tau-;:OnOff Off 
set /Herwig/Particles/W-/W-->ubar,d;:OnOff %s 
set /Herwig/Particles/W-/W-->s,ubar;:OnOff %s
set /Herwig/Particles/W-/W-->b,cbar;:OnOff %s 
set /Herwig/Particles/W-/W-->cbar,d;:OnOff %s
set /Herwig/Particles/W-/W-->cbar,s;:OnOff %s 
""" % (Wm_lep, Wm_had, Wm_had, Wm_had, Wm_had, Wm_had)
  

#set /Herwig/Particles/Z0/Z0->u,ubar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->b,bbar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->c,cbar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->d,dbar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->e-,e+;:OnOff On 
#set /Herwig/Particles/Z0/Z0->mu-,mu+;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->nu_e,nu_ebar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->s,sbar;:OnOff Off 
#set /Herwig/Particles/Z0/Z0->tau-,tau+;:OnOff Off 



    text += """
set /Herwig/Particles/t/t->nu_e,e+,b;:OnOff Off 
set /Herwig/Particles/t/t->nu_mu,mu+,b;:OnOff %s 
set /Herwig/Particles/t/t->nu_tau,tau+,b;:OnOff Off 

set /Herwig/Particles/t/t->b,u,dbar;:OnOff %s
set /Herwig/Particles/t/t->b,c,sbar;:OnOff %s 
set /Herwig/Particles/t/t->b,sbar,u;:OnOff %s 
set /Herwig/Particles/t/t->b,c,dbar;:OnOff %s 
set /Herwig/Particles/t/t->b,bbar,c;:OnOff %s 
""" % (t_lep, t_had, t_had, t_had, t_had, t_had)

    text += """
set /Herwig/Particles/tbar/tbar->nu_ebar,e-,bbar;:OnOff Off 
set /Herwig/Particles/tbar/tbar->nu_mubar,mu-,bbar;:OnOff %s 
set /Herwig/Particles/tbar/tbar->nu_taubar,tau-,bbar;:OnOff Off 
set /Herwig/Particles/tbar/tbar->bbar,ubar,d;:OnOff %s 
set /Herwig/Particles/tbar/tbar->bbar,cbar,s;:OnOff %s 
set /Herwig/Particles/tbar/tbar->bbar,s,ubar;:OnOff %s 
set /Herwig/Particles/tbar/tbar->bbar,cbar,d;:OnOff %s 
set /Herwig/Particles/tbar/tbar->b,bbar,cbar;:OnOff %s 
""" % (tbar_lep, tbar_had, tbar_had, tbar_had, tbar_had, tbar_had)

    return text


#####
def writeFinalAnalysisHandlerInfo(runname, runJetAnalysis=False):
    if runJetAnalysis:
        ss= '''
create HerwigWrapper::Evgen Evgen Evgen.so
insert LHCGenerator:AnalysisHandlers 0 Evgen'''
    else:
        ss = '''
insert /Herwig/Generators/LHCGenerator:AnalysisHandlers[0] /Herwig/Analysis/HepMCFile
set /Herwig/Analysis/HepMCFile:PrintEvent 999999999
set /Herwig/Analysis/HepMCFile:Filename %s.hepmc
set /Herwig/Analysis/HepMCFile:Format GenEvent
set /Herwig/Analysis/HepMCFile:Units GeV_mm'''  % runname

    return """
#set /Herwig/Analysis/HepMCFile:PrintEvent 9999999999
##################################################
# Analysis Handler
##################################################
cd /Herwig/Generators
%s
saverun %s LHCGenerator

""" % (ss, runname)

#####
def makeRunname(metype):
    import socket
    import os
    runname = '_'.join((metype, 
                        socket.gethostname(), 
                        str(os.getpid()),
                        time.strftime('%b-%d-%Y-%H-%M-%S')))

    return runname

def createHerwigRunCard(outputfileName='run.in', runname='LHC-v', metype='v', lepcharge=1, meminpt=200, memaxpt=20000, seed=0, nevents=10):

    

    runcardText = writeBasicInfo(nevents, seed) \
        + writePDFInfo() \
        + writeUEtune() \
        + writeMatrixElementInfo(metype) \
        + writeKinematicCutInfo(metype,meminpt,memaxpt) \
        + writeParticleDecays(metype,lepcharge) \
        + writeFinalAnalysisHandlerInfo(runname)


    outfile = open(outputfileName, 'w')
    outfile.write(runcardText)
    outfile.close()

    return runname

def readRunCard(cardfileName = 'run.in'):
    print "readRunCard %s:" % cardfileName
    args = ["Herwig++", "read", cardfileName]
    args = " ".join(args)
    print args
    result = subprocess.call(args, shell=True)
    return result

def runHerwig(configName):
    print "runHerwig %s:" % configName
    args = ['Herwig++', 'run', configName]
    args = " ".join(args)
    print args
    result = subprocess.call(args, shell=True)
    return result

def runHepMCJets(seed, metype, lepcharge
                 , savemaxnjets, nevents
                 , cellv, trimjet
                 , jetminpt, jetmaxpt
                 , jetminm, jetmaxm
                 , rad, mu, smear
                 , wtagger, input_file
                 ):

    print "runHepMCJets:"

    args = ["hepmcJets.exe"]
    args += ["--seed=%s" % str(seed)                 , "--source=Herwig"                   ,
             "--metype=%s" % str(metype)             , "--lepcharge=%s" % str(lepcharge)   ,
             "--savemaxnjets=%s" % str(savemaxnjets) , "--nevents=%s" % str(nevents)       ,
             "--cellv=%s" % str(cellv)               , "--trimjet=%s" % str(trimjet)       ,
             "--jetminpt=%s" % str(jetminpt)         , "--jetmaxpt=%s" % str(jetmaxpt)     ,
             "--jetminm=%s" % str(jetminm)           , "--jetmaxm=%s" % str(jetmaxm)       ,
             "--rad=%s" % str(rad)                   , "--mu=%s" % str(mu)                 ,
             ]
    
    if smear:
        args.append("--smear")
    if wtagger:
        args.append("--wtagger")

    args.append(input_file)

    args = " ".join(args)
    print args
    result = subprocess.call(args, shell=True)
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
    
    parser = argparse.ArgumentParser(description='Create a Herwig++ run card.')
    parser.add_argument('--seed'      , default=-1       , type=int                                         , help='The random number generator seed (-1 finds a random seed)')
    parser.add_argument('--savemaxnjets', default=-1     , type=int                                         , help='Max number of jets to save.  Default, -1, saves all jets (which satisfy other selectons)')
    parser.add_argument('--nevents'   , default=100      , type=int                                         , help='The generator number of events')
    parser.add_argument('--metype'    , default='v'      , type=str   , choices = ['w', 'v', 't', 'q', 'b'] , help='Matrix element type')
#    parser.add_argument('--cardName'  , default='run.in' , type=str                                         , help='Output file name for Herwig++ run card name')
    parser.add_argument('--lepcharge' , default=1        , type=int   , choices = [-1, 1, 0]                , help='Charge of lepton desired in final state, only for ME ={v, t, W}.  1=l^+, -1=l^-, 0= All lepton flavors, any other value = no lepton (all hadronic)')
    parser.add_argument('--meminpt'   , default=200      , type=int                                         , help='Herwig++ MinKT (GeV) for v, t, and q')
    parser.add_argument('--memaxpt'   , default=20000    , type=int                                         , help='Herwig++ MaxKT (GeV) for v, t, and q')
    parser.add_argument('--jetminpt'  , default=200      , type=int                                         , help='Jet PTMin (GeV)')
    parser.add_argument('--jetmaxpt'  , default=20000    , type=int                                         , help='Jet PTMax (GeV)')
    parser.add_argument('--jetminm'   , default=200      , type=int                                         , help='Jet Min Mass (GeV)')
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
    print "makeRunname"
    runname = makeRunname(metype=args.metype)

    sys.stdout.flush()
    print ""
    print "-+"*50
    print "createHerwigRunCard"
    print createHerwigRunCard(outputfileName = runname + '.in', runname = runname
                        , metype = args.metype, lepcharge = args.lepcharge
                        , meminpt = args.meminpt, memaxpt = args.memaxpt
                        , seed = args.seed, nevents = args.nevents
                        )
    
    print subprocess.call(['which', 'Herwig++'])

    sys.stdout.flush()
    print ""
    print "-+"*50
    print "readRunCard"
    print readRunCard(runname + '.in')

    sys.stdout.flush()
    print ""
    print "-+"*50
    print "runHerwig"
    print runHerwig(runname + '.run')

    sys.stdout.flush()
    print ""
    print "-+"*50
    print "runHepMCJets"
    print runHepMCJets(seed=args.seed, metype=args.metype, lepcharge=args.lepcharge
                       , savemaxnjets=args.savemaxnjets, nevents=args.nevents
                       , cellv=args.cellv, trimjet=args.trimjet
                       , jetminpt=args.jetminpt, jetmaxpt=args.jetmaxpt
                       , jetminm=args.jetminm, jetmaxm=args.jetmaxm
                       , rad=args.rad, mu=args.mu, smear=args.smear
                       , wtagger=args.wtagger, input_file=runname + '.hepmc'
                      
                       )
    sys.stdout.flush()
