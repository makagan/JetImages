#
# $Id: pyrootmagic.py 340230 2013-10-15 00:37:57Z jcogan $
#
'''Make adding branches work in PyROOT for vectors and scalars
This code is experts only'''

__author__ = "Josh Cogan (jcogan@cern.ch), Emanuel Strauss (estrauss@slac.stanford.edu)"
__version__ = '$Revision: 340230 $'
__date__    = '$Date: 2013-10-15 02:37:57 +0200 (Tue, 15 Oct 2013) $'

import os
try:
    import simpleLogging
    scribe = simpleLogging.getChildLogger('main', __name__)
except ImportError:
    import logging
    scribe = logging.getLogger(__name__)

typeMap =   { 'O': ('bool'          , False )
            , 'I': ('int'           , -999  )
            , 'i': ('unsigned int'  , 999   )
            , 'S': ('short'         , -999  )
            , 's': ('unsigned short', 999   )
            , 'D': ('double'        , -999  )
            , 'F': ('float'         , -999  )
            }



### ITS VERY IMPORTANT THAT THESE MATCH ROOT DOCUMENTATION
typeRegexs = [
                (r'(unsigned\s+int|UInt_t)'     , 'i'),
                (r'(int|Int_t)'                 , 'I'),
                (r'(bool|Bool_t)'               , 'O'),
                (r'(short|Short_t)'             , 'S'),
                (r'(unsigned\s+short|UShort_t)' , 's'),
                (r'(float|Float_t)'             , 'F'),
                (r'(double|Double_t)'           , 'D'),
            ]

_namesTaken = set()
_nameConversions = [(
                     ('smart_el', 'smart_el'),
                     ('el_GSF'  , 'smart_el'),
                     ('el_gsf'  , 'smart_el'),
                     ('el'      , 'smart_el'),
                    ),
                    ( # this is a group of mutually exlusive redefs
                     ('smart_jet'               , 'smart_jet'),
                     ('jet_AntiKt4TopoEMJets'   , 'smart_jet'),
                     ('jet_akt4topoem'          , 'smart_jet'),
                     ('jet_AntiKt4TopoEM'       , 'smart_jet'),
                    ),
                    (
                     ('smart_mu'            , 'smart_mu'),
                     ('mu_staco'            , 'smart_mu'),
                     ('mu'                  , 'smart_mu'),
                    )
                   ]

## Generate the output tree code
def makeStructString(varNames, sname = "MyStruct"):
    '''Produce the ROOT/C code defining all the variables to be filled'''
    defs  = []
    inits = []
    for var, varDesc in varNames:
        cType = typeMap[varDesc.split(':')[0]][0]

        if ':' not in varDesc:
            line = '%s %s' % (cType, var)
        else:
            continue

        defs.append(line)

    s = "struct  %s{\n\t" % sname
    s += ';\n\t'.join(defs)
    s += ';\n};\n'

    if len(inits) > 0:
        s += ';\n'.join(inits) + ';'

    return  s

def typeLookUp(inTree, branch):
    import re

    ll = inTree.GetLeaf(branch)
    if ll is None:
        return None

    try:
        sdesc = ll.GetTypeName()
    except ReferenceError:
        scribe.error('RefError')
        return None

    #need to make this nVector
    nVectors = len(re.findall(r'vector\s*<', sdesc))
    for regStr, retVal in typeRegexs:
        if re.search(regStr, sdesc) is not None:
            return retVal + (':' * nVectors)

    scribe.error('%s %s failed' % branch, sdesc)
    return None

def inputTreeVars(inTree):
    try:
        branches = [xx.GetName() for xx in inTree.GetListOfBranches()]
    except TypeError as ee:
        scribe.critical("inTree, %s, has no braches" % str(inTree))
        raise ee

    branches = filter(lambda br: inTree.GetBranchStatus(br), branches)
    thisTreesNames = set(branches)

    #branches = doConversions(branches, thisTreesNames)
    newNameMap = {}
    for exclusiveGroup in _nameConversions:
        for prefix, unifiedName in exclusiveGroup:
            matchingBranches = [br for br in thisTreesNames if br.startswith(prefix)]
            if len(matchingBranches) < 10:
                #no branches match el_GSF, its okay to check el_vanilla
                #if just one branch does match el_GSF stop trying to match
                #ignore the few straggling matches, look for a block
                continue

            for mbr in matchingBranches:
                oldType = typeLookUp(inTree, mbr)
                if oldType is None:
                    scribe.warning('Couldnt ascertain type of %s.  Skipping it!' % mbr)
                    continue

                newNameMap[mbr] = (unifiedName + mbr[len(prefix):], oldType)
            
            #since a branch matched el_GSF and got mapped to smart_el
            #don't allow any el_vanilla branches to also get mapped to smart_el
            break
        
    remappedVars = set(newNameMap.keys())
    for br in thisTreesNames:
        if br in remappedVars:
            continue
        brType = typeLookUp(inTree, br)
        if brType is None:
            scribe.warning('Couldnt ascertain type of %s.  Skipping it!' % br)
        else:
            newNameMap[br] = (br, brType)

    return newNameMap

def bindInputToStruct(inTree, structObj, inTreeVars, tempDir='./'):
    import ROOT as rt
## this might not be necessary!!
## grep GenerateDictionary main_analysis.py
    loaderText ='''// File loader.C
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<unsigned int> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<float> >+;
#endif
'''

    fLoc = os.path.join(tempDir, 'prmLoader.C')
    open(fLoc, 'w').write(loaderText)
    rt.gROOT.ProcessLine('.L %s+'% fLoc)

    branches = []
    for oldName in inTreeVars:
        var, varDesc = inTreeVars[oldName]
        inBranch = inTree.GetBranch(oldName)

        if ':' in varDesc:
            inTree.SetBranchAddress(oldName, getattr(structObj, var))
            pass
        else:
            inTree.SetBranchAddress(oldName, rt.AddressOf(structObj, var))
            pass

def getAndLoadStruct(varNames, inTree=None, force_recompile=False, sname="MyStruct"):
    ''' Produce a python object (struct) correctly interfaced with ROOT'''

    #force_recompile = True
    import subprocess
    import ROOT as rt
    import shutil
    import os.path
    import tempfile


    inTreeVars = {}
    if inTree is not None:
        inTreeVars = inputTreeVars(inTree)

    allVarNames = varNames + inTreeVars.values()

    if force_recompile:
        tempDirName = tempfile.mkdtemp('_pyrootmagic')
    else:
        ##appears ROOT can't load files name ./blah.so only blah.so
        tempDirName = os.getcwd()

    fileBase = os.path.join(tempDirName, 'structCode')

    s = '#include <vector>\n%s\n' % makeStructString(allVarNames, sname)
    f = open(fileBase+'.h', "w")
    f.write(s)
    f.close()

    # ARA uses cintex, which uses reflex, which clashes with the
    # vector dictionaries that are made with LoadMacro
    # Gaaah!
    #if True:
    gccXMLPath='--gccxmlpath='+os.path.expanduser('~joshgc/public/cern_style_gccxml/i686-slc5-gcc43-opt/bin')

    scribe.debug("Forcing recompile?  %s" % force_recompile)
    scribe.debug("Shared object exists? %s" % os.path.exists(fileBase+'.so'))
    if not os.path.exists(fileBase+'.so') or force_recompile:
        scribe.warning('Need to rerun struct compilation on %s' % fileBase)
        subprocess.call("genreflex %s.h -o %s.cc -I $ROOTSYS/include --debug=1 %s" 
                        % (fileBase, fileBase, gccXMLPath), shell=True)
        subprocess.call(("g++ %s.h %s.cc -o %s.so `root-config --cflags --libs`"
                        " -lReflex -lCintex -shared")
                        % (fileBase, fileBase, fileBase), shell=True)
    else:
        scribe.warning('Skipping struct compilation')


    import PyCintex
    PyCintex.Cintex.Enable()
    rt.gSystem.Load(fileBase+'.so')

    structObj = getattr(rt, sname)()
    structObj.__dict__['_branches'] = list(allVarNames)

    # TBranch changes the vectors in the struct, so we use the hack below
    # Double Gaaah!
    for var, varDesc in allVarNames:
        if ':' in varDesc:
            cls = typeMap[varDesc.split(':')[0]][0]
            for ii in range(varDesc.count(':')):
                cls = rt.std.vector(cls)

            setattr(structObj, var, cls())

    bindInputToStruct(inTree, structObj, inTreeVars, tempDirName)

    if force_recompile:
        try:
            shutil.rmtree(tempDirName)
        except OSError:
            scribe.error('Unable to delete temp directory %s' % tempDirName)
            scribe.error('This will be bad later on.  Tell Josh.')


    return structObj, allVarNames

def clearStruct(sObj):
    '''Set all the struct variables to their defaults.
        These are defined in the type map.
    '''

    #do i have to handle the vector vector stuff here too?
    for var, varType in sObj._branches:
        if ':' in varType:
            getattr(sObj, var).clear()
            pass
        else:
            #this is probably slow, always go to 999?
            setattr(sObj, var, typeMap[varType.split(':')[0]][1])


def mapVariablesToBranches(vstruct, tree_out, listOfVars):
    import ROOT as rt
    for var, varDesc in listOfVars:
        cType = varDesc.split(':')[0]
        if ':' not in varDesc:
            tree_out.Branch(var, rt.AddressOf(vstruct, var), '%s/%s' % (var, cType))
            pass
        else:
            tree_out.Branch(var, "vector<%s>" % typeMap[cType][0], getattr(vstruct, var))
            pass
    pass #end of mapVariablesToBranches

listOfVars =[
            # Indices of Jets suriving overlap removal
              ('jet_index'         , 'I:')
            # Indices for H->Z(ll)ph reconstruction
            # One per reconstruction
            , ('l1_index'           , 'I:')
            , ('l2_index'           , 'I:')
            , ('ph_index'           , 'I:')
            
            # Event level variables for numbers of objects, indices for best candidates
            # One per event
            , ('n_ee'               , 'I')
            , ('n_mm'               , 'I')
            , ('n_jet'              , 'I')
            , ('n_tau'              , 'I')
            , ('best_ee_index'      , 'I')
            , ('best_mm_index'      , 'I')
            
            , ('truth_H_index'      , 'I')
            , ('truth_Z_index'      , 'I')
            , ('truth_ph_index'     , 'I')
            , ('truth_l1_index'     , 'I')
            , ('truth_l2_index'     , 'I')
            
            # Branch for keep track of species of particle and best candidates
            # One per reconstruction
            , ('mode'               , 'I:')

            #these appear to be never used in main_analysis.
            #Looks like is_best_Hee is preferred
            #, ('is_best_ee'         , 'I:')
            #, ('is_best_mm'         , 'I:')
            
            # Z(ll) and H(ll)ph kinematics
            # One per reconstruction
            , ('Z_pt'               , 'D:')
            , ('Z_eta'              , 'D:')
            , ('Z_phi'              , 'D:')
            , ('Z_E'                , 'D:')
            , ('Z_m'                , 'D:')
            , ('H_pt'               , 'D:')
            , ('H_eta'              , 'D:')
            , ('H_phi'              , 'D:')
            , ('H_E'                , 'D:')
            , ('H_m'                , 'D:')
            , ('H_cosCapTheta'      , 'D:')
            , ('H_cosTheta'         , 'D:')
            , ('H_cosBkgPhi'        , 'D:')
            , ('H_truth_cosCapTheta', 'D:')
            , ('H_truth_cosTheta'   , 'D:')
            , ('H_truth_cosBkgPhi'  , 'D:')
            , ('H_pT_t'             , 'D:')
            , ('H_category'         , 'I:')
            
            , ('brem_none_Z_m'      , 'D')
            , ('brem_none_H_m'      , 'D')
            , ('brem_all_Z_m'       , 'D')
            , ('brem_one_Z_m'       , 'D')
            , ('brem_one_H_m'       , 'D')
            , ('brem_clus_dr'       , 'D')
            , ('brem_clus_pt'       , 'D')
            , ('brem_clus_il'       , 'I')
            , ('brem_pre_mu_eta'    , 'D')
            , ('brem_pre_mu_pt'     , 'D')

            , ('brem_mass_llhood'   , 'D')
            , ('brem_ang_llhood'    , 'D')
            , ('brem_llhood'        , 'D')

            # DeltaR and DeltaPhi variables for H->Z(ll)ph reconstructions
            # One per reconstruction
            , ('l1l2_DR'            , 'D:')
            , ('l1ph_DR'            , 'D:')
            , ('l2ph_DR'            , 'D:')
            , ('l1Z_DR'             , 'D:')
            , ('l1H_DR'             , 'D:')
            , ('l2Z_DR'             , 'D:')
            , ('l2H_DR'             , 'D:')
            , ('phZ_DR'             , 'D:')
            , ('phH_DR'             , 'D:')
            , ('ZH_DR'              , 'D:')
            , ('l1l2_DPhi'          , 'D:')
            , ('l1ph_DPhi'          , 'D:')
            , ('l2ph_DPhi'          , 'D:')
            , ('l1Z_DPhi'           , 'D:')
            , ('l1H_DPhi'           , 'D:')
            , ('l2Z_DPhi'           , 'D:')
            , ('l2H_DPhi'           , 'D:')
            , ('phZ_DPhi'           , 'D:')
            , ('phH_DPhi'           , 'D:')
            , ('ZH_DPhi'            , 'D:')
            , ('ph_jet_DR'          , 'D:')
            , ('ph_tau_DR'          , 'D:')
            , ('ph_jet_index'       , 'I:')
            , ('ph_tau_index'       , 'I:')
            
            # Kinematic variables for objects so we don't have to look them up
            , ('l1_pt'              , 'D:')
            , ('l1_eta'             , 'D:')
            , ('l1_phi'             , 'D:')
            , ('l1_E'               , 'D:')
            , ('l1_charge'          , 'D:')
            , ('l2_pt'              , 'D:')
            , ('l2_eta'             , 'D:')
            , ('l2_phi'             , 'D:')
            , ('l2_E'               , 'D:')
            , ('l2_charge'          , 'D:')
            , ('H_ph_pt'            , 'D:')
            , ('H_ph_eta'           , 'D:')
            , ('H_ph_phi'           , 'D:')
            , ('H_ph_E'             , 'D:')
            
            , ('l1_etcone20'        , 'D:')
            , ('l1_etcone30'        , 'D:')
            , ('l1_etcone40'        , 'D:')
            , ('l1_nucone20'        , 'D:')
            , ('l1_nucone30'        , 'D:')
            , ('l1_nucone40'        , 'D:')
            , ('l1_ptcone20'        , 'D:')
            , ('l1_ptcone30'        , 'D:')
            , ('l1_ptcone40'        , 'D:')
            , ('l2_etcone20'        , 'D:')
            , ('l2_etcone30'        , 'D:')
            , ('l2_etcone40'        , 'D:')
            , ('l2_nucone20'        , 'D:')
            , ('l2_nucone30'        , 'D:')
            , ('l2_nucone40'        , 'D:')
            , ('l2_ptcone20'        , 'D:')
            , ('l2_ptcone30'        , 'D:')
            , ('l2_ptcone40'        , 'D:')
            , ('H_ph_etcone20'      , 'D:')
            , ('H_ph_etcone30'      , 'D:')
            , ('H_ph_etcone40'      , 'D:')
            , ('H_ph_ptcone20'      , 'D:')
            , ('H_ph_ptcone30'      , 'D:')
            , ('H_ph_ptcone40'      , 'D:')
            
            # Object tightness
            , ('l1_tightness'       , 'I:')
            , ('l2_tightness'       , 'I:')
            , ('ph_tightness'       , 'I:')
            
            , ('is_best_Hee'        , 'I:')
            , ('is_best_Hmm'        , 'I:')
            
            # Auxiliary variables
            # One per reconstruction
            , ('l1ph_m'             , 'D:') # Invariant mass of a lepton and the photon
            , ('l2ph_m'             , 'D:') # Invariant mass of a lepton and the photon
            
            # Truthmatching
            # One per reconstruction
            # Keep these as ints, we may want to contain more information than a bool at some point
            # H->Z(ee)ph
            , ('real_l1'            , 'I:') # =1 if the first  electron candidate in a recontsruction is really an electron
            , ('real_l2'            , 'I:') # =1 if the second electron candidate in a recontsruction is really an electron
            , ('real_ph'            , 'I:') # =1 if the photon candidate in a recontsruction is really a photon
            , ('l1_from_Z'          , 'I:') # =1 if the first  electron is real and comes from a real Z
            , ('l2_from_Z'          , 'I:') # =1 if the second electron is real and comes from a real Z
            , ('same_Z'             , 'I:') # =1 if both electrons come from the same real Z
            , ('Z_from_H'           , 'I:') # =1 if both electrons come from the same real Z and this Z comes from a real H
            , ('ph_from_H'          , 'I:') # =1 if the photon comes from a real Higgs
            , ('same_H'             , 'I:') # =1 if both electrons come from the same real Z, and that Z and the photon come from a real Higgs
            , ('ph_source'          , 'I:') # integer (enum-esque) indicating genre of photon
            , ('ph_source_weight'   , 'I:') # 0 if this photon is overlaps with sherpa ME calculation
            , ('ph_parent_index'    , 'I:') # mc_index of parent of photon
            , ('ph_parent'          , 'I:') # pdgId of parent of photon
            , ('ph_grandparent'     , 'I:') # pdgId of parent of parent of photon

            # Weights
            # Some weights are reconstruction dependent
            , ('pileup_pass'        , 'O' )
            , ('pileup_weight'      , 'D' )
            , ('trigger_weight'     , 'D' )
            , ('crosssection_weight', 'D' )
            , ('blind_weight'       , 'D:')
            , ('weight'             , 'D:') # Never used??
            , ('SF'                 , 'D:') # Multiplicative scale factor
            
            # Definitions of objects
            , ('is_cutflow_1'        , 'I:')
            , ('is_cutflow_2'        , 'I:')
            , ('H_definition_bitmap' , 'I:')
            , ('l1_definition_bitmap', 'I:')
            , ('l2_definition_bitmap', 'I:')
            , ('ph_definition_bitmap', 'I:')
            ]

for lepType in ('both', 'only_ph', 'only_el'):
    for clusMin in ('1000', '3500'):
        for vv in ('clus_pt', 'clus_dr', 'clus_il', 'one_Z_m', 'one_H_m'):
            listOfVars.append(('brem_%s_%s_%s' % (vv, lepType, clusMin), 'D'))

if __name__ == '__main__':
    import ROOT as rt

    inChain = rt.TChain("physics")
    inChain.Add("~arconde/dq2/HToZGamma/mc11_7TeV.128988.PowHegPythia_ggH125_Zgam_Zmumu.merge.NTUP_HSG2.e1063_s1372_s1370_r3043_r2993_p869_tid729809_00/*")
    #inChain.Add("output/skimmed/data11_7TeV_177986_physics_Egamma_r2603_p659_p761_p762.root")
    inChain.SetBranchStatus("*", 0)
    inChain.SetBranchStatus("mc_*", 1)
    #inChain.SetBranchStatus("el_*n", 1)
    #inChain.SetBranchStatus("el_*loose", 1)
    #inChain.SetBranchStatus("el_*eta", 1)

    ff = rt.TFile('output.root', 'recreate')
    tt = rt.TTree('arg', 'arg')

    vstruct, newListOfVars = getAndLoadStruct([], inChain)
    #vstruct, newListOfVars = getAndLoadStruct(listOfVars, inChain)

    mapVariablesToBranches(vstruct, tt, newListOfVars)
    for ii in range(5):
        clearStruct(vstruct)
        inChain.GetEntry(ii)
        print 'Entry %d ---------------------' % ii
        print 'direct', inChain.mc_n, inChain.mc_parents.size(),inChain.mc_parents[0].size()
        print 'struct', vstruct.mc_n, inChain.mc_parents.size(),inChain.mc_parents[0].size()
    print 'Calling destructors!'
