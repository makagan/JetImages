#!/usr/bin/env python

__author__ = 'Josh Cogan'

import os
import sys
import time
import InputCleaner as ic
#G{'ME': 'g', 'e': '62.3739', 'eta': '-0.519661', 'eventNumber': '0', 'm': '11.2668', 'nParts': '16', 'npv': '1', 'pdgIDHardParton': '21', 'phi': '2.17741', 'pt': '53.9043', 'radius': '0.4', 'rap': '-0.509537', 'seed': '10105', 'source': 'Pythia8', 'subjetRadius': '0.3', 'cells': [(0, 4.16334e-17, 8.32667e-17),(0, 0.1, 8.32667e-17),(0, 0.2, 8.32667e-17),(0, 0.3, 8.32667e-17),(0, 0.4, 8.32667e-17),(0, 0.5, 8.32667e-17),(0, 0.6, 8.32667e-17),(0, 0.7, 8.32667e-17),(0, 0.8, 8.32667e-17),(0, 4.16334e-17, 0.0997331),(0, 0.1, 0.0997331),(2.89078, 0.2, 0.0997331),(0, 0.3, 0.0997331),(0, 0.4, 0.0997331),(0, 0.5, 0.0997331),(0, 0.6, 0.0997331),(0, 0.7, 0.0997331),(0, 0.8, 0.0997331),(0, 4.16334e-17, 0.199466),(0, 0.1, 0.199466),(0, 0.2, 0.199466),(2.53236, 0.3, 0.199466),(1.06125, 0.4, 0.199466),(0, 0.5, 0.199466),(0, 0.6, 0.199466),(0, 0.7, 0.199466),(0, 0.8, 0.199466),(0, 4.16334e-17, 0.299199),(0, 0.1, 0.299199),(0.24269, 0.2, 0.299199),(18.4752, 0.3, 0.299199),(5.57112, 0.4, 0.299199),(0.480723, 0.5, 0.299199),(0, 0.6, 0.299199),(3.49651, 0.7, 0.299199),(2.14862, 0.8, 0.299199),(0, 4.16334e-17, 0.398932),(0, 0.1, 0.398932),(1.48241, 0.2, 0.398932),(0, 0.3, 0.398932),(0, 0.4, 0.398932),(0, 0.5, 0.398932),(0, 0.6, 0.398932),(0, 0.7, 0.398932),(0, 0.8, 0.398932),(0, 4.16334e-17, 0.498666),(0, 0.1, 0.498666),(0, 0.2, 0.498666),(0, 0.3, 0.498666),(1.24219, 0.4, 0.498666),(13.8829, 0.5, 0.498666),(0, 0.6, 0.498666),(1.23478, 0.7, 0.498666),(0, 0.8, 0.498666),(0, 4.16334e-17, 0.598399),(0, 0.1, 0.598399),(0, 0.2, 0.598399),(0, 0.3, 0.598399),(0, 0.4, 0.598399),(2.70018, 0.5, 0.598399),(4.10979, 0.6, 0.598399),(0, 0.7, 0.598399),(0, 0.8, 0.598399),(0, 4.16334e-17, 0.698132),(0, 0.1, 0.698132),(0, 0.2, 0.698132),(0, 0.3, 0.698132),(0, 0.4, 0.698132),(0.822497, 0.5, 0.698132),(0, 0.6, 0.698132),(0, 0.7, 0.698132),(0, 0.8, 0.698132),(0, 4.16334e-17, 0.797865),(0, 0.1, 0.797865),(0, 0.2, 0.797865),(0, 0.3, 0.797865),(0, 0.4, 0.797865),(0, 0.5, 0.797865),(0, 0.6, 0.797865),(0, 0.7, 0.797865),(0, 0.8, 0.797865),], 'n_cells': '81','jet_offset': (-0.95, 1.7952), 'subjets': [(26.9217, -0.641016, 2.07111),(21.8911, -0.427662, 2.32936),],'rel_subjets': [(26.9217, 0.308984, 0.275916),(21.8911, 0.522338, 0.534161),], 'pu_energies': [], 'pu_energy': 0}

_stringToInts = {'source': {'Pythia8':0, 'Herwig':1, 'HepMC':1},
                 'ME'    : {'q': 1, 'g':21, 'w':24, 'v':2, 'c':4, 'b': 5, 't':6, 'h': 25},
                }

_indexToNames = {
                #'parts'         : [('pdgId','i4'),('E' ,'f8'),('cosT','f8'),('phi','f8')],
                'cells'         : [('E' ,'f8'), ('eta', 'f8') , ('phi','f8')],
                'subjets'       : [('E' ,'f8'), ('eta', 'f8') , ('phi','f8')],
                'rel_subjets'   : [('E' ,'f8'), ('eta', 'f8') , ('phi','f8')],
                }

_integer_whitelist = set(('pdgIDHardParton', 'eventNumber', 'seed', 'nParts', 'npv'))

def handleCLI():
    from optparse import OptionParser

    usage = '%prog [options] [-o] outputFile input1 input2 ...'

    pp = OptionParser(usage=usage)
    pp.add_option('-f', '--format', dest='format', default=None,
                    help='Set the output format to either root or numpy.'
                         '  If not set, the ouput is a ROOT file if it ends in'
                         ' .root (case sensitive).')

    pp.add_option('-o', '--output', dest='output', default=None,
                    help='Set the output filename.  If not specified then'
                         ' the first filename after the options is used.')
    pp.add_option('-e', '--end', dest='end', type=int, default=-1,
                    help='Specify maximum number of rows to attempt to process'
                         ' -1 processes all rows available')

    pp.add_option('-s', '--search', dest='search', default=None,
                    help='Regex with possible grouping.  Each line in'
                    ' the input text must match this search before it gets'
                    ' converted.  If there are groups in this regex, then'
                    ' they may be used by the --repl option.  For'
                    ' example, some files may require -s \'r"^#G(.*)$"\''
                    ' -r \'r"\\1"\'.')

    pp.add_option('-r', '--repl', dest='repl', default=None,
                    help='Replacement regex with possible grouping. '
                    '  See --search')

    opts, args = pp.parse_args()

    if opts.output is None:
        try:
            opts.output = args.pop(0)
        except IndexError:
            print 'Please specify an output file!'
            pp.print_help()
            sys.exit(1)

    if len(args) == 0:
        print 'Please specify one or more input files!'
        pp.print_help()
        sys.exit(1)

    if opts.format is None:
        print opts.output
        if '.root' == os.path.splitext(opts.output)[1]:
            opts.format = 'root'
        else:
            opts.format = 'npy'

    return opts, args

def makeListOfVariables(record):
    lov = []
    #usedNames = set()
    #cannot trust InputCleaner to get 1.0 vs 1 right.
    #see makeListOfStructVariables
    print 'I WILL SCREW UP INT VS FLOAT'
    for kk, vv in record.iteritems():
        if isinstance(vv, list) or  isinstance(vv, tuple):
            try:
                nested = isinstance(vv[0], list) or isinstance(vv[0], tuple)
            except IndexError:
                nested = False

            if 'n_%s' % kk not in record.keys():
                lov.append(['n_%s' % kk, 'I'])

            if nested:
                #this should RIGHTLY break if kk not in _indexToNames!!
                for index, nameAndType in enumerate(_indexToNames[kk]):
                    name, tt = nameAndType
                    lov.append([kk + '_' + name, 'D:' if tt[0]=='f' else 'I:'])
            else:
                try:
                    tt = 'I' if isinstance(vv[0], int) else 'D'
                except IndexError: 
                    tt = 'D'
                lov.append([kk, tt+':'])

        elif isinstance(vv, str  ): lov.append([kk, 'I'])
        elif isinstance(vv, int  ): lov.append([kk, 'I'])
        elif isinstance(vv, float): lov.append([kk, 'D'])
        else:
            print "I don't know how to parse key value pair: %s, %s" % (kk, vv)
            sys.exit(1)

    return lov

def makeListOfStructVariables(record):
    jetType    = []
    #don't be so "smart" about type intuition
    easymap = {str: 'i4', int: 'f8', float: 'f8'}

    for kk, vv in record.iteritems():
        if 'n_' == kk[:2] or kk in _integer_whitelist:
            #don't be so "smart" about type intuition
            jetType.append((kk, int))
            continue

        try:
            jetType.append((kk, easymap[type(vv)]))
            continue
        except KeyError:
            pass

        #okay vv is list or tuple
        if kk in _indexToNames:
            ##this should return 10 for subjets 
            ##and the full length for cells
            nelems = max(10, len(vv)) 
            jetType.append((kk, _indexToNames[kk], nelems))
            if 'n_%s' % kk not in record.keys():
                jetType.append(('n_%s' % kk, 'i4'))
            continue

        if 'offset' in kk:
            jetType.append((kk, 'f8', 2))
        else: #puenergies
            jetType.append((kk, 'f8', 100))
            if 'n_%s' % kk not in record.keys():
                jetType.append(('n_%s' % kk, 'i4'))


    return jetType

def makeRootFile(fname, incoming, max_elems):
    first = [rr for rr in incoming.getGen(1)][0]
    listOfVars = makeListOfVariables(first)
    print 'Will produce a ROOT file with these fields ['
    for name, cType in sorted(listOfVars, key=lambda li: li[0]):
        print '    %s, %s' % (name.ljust(30), cType)
    print ']'

    import ROOT as rt
    fOut = rt.TFile(fname, 'recreate')
    tree_out = rt.TTree("physics", "physics")

    import pyrootmagic as prm
    vs, varsWithInTrees = prm.getAndLoadStruct(listOfVars, None, force_recompile=True)
    print 'pyroot struct written, compiled and loaded'
    ### Struct can't be garbage collected before tree
    ### http://root.cern.ch/phpBB3/viewtopic.php?f=14&t=8585
    tree_out._vs = vs
    prm.mapVariablesToBranches(vs, tree_out, listOfVars)

    printPoint = 0.1
    nRecord = incoming.size()

    print "BLARGH still need to handle n_subjets with subjets autogenerated??"
    print "BLARGH still need to handle n_subjets with subjets autogenerated??"
    print "BLARGH still need to handle n_subjets with subjets autogenerated??"
    for iRecord, record in enumerate(incoming.getGen(max_elems)):
        if iRecord % (nRecord * printPoint) == 0:
            print 'Working on the %d%%th event' % (100*iRecord/nRecord)
        prm.clearStruct(vs)
        for kk, vv in record.iteritems():
            if kk in _indexToNames: # its the pdgid/pt/eta/phi quadruplet
                if 'n_%s' % kk not in record.keys():
                    pass
                    #setattr(vs, 'n_%s' % kk, len(record[kk]))
                for vIndex, nameAndType in enumerate(_indexToNames[kk]):
                    name = nameAndType[0]
                    vec = getattr(vs, kk+'_'+name)
                    for partIndex in range(len(record[kk])):
                        vec.push_back(vv[partIndex][vIndex])

            elif kk in _stringToInts:
                setattr(vs, kk, _stringToInts[kk][vv])
            else:
                setattr(vs, kk, vv)

        tree_out.Fill()

    fOut.Write()
    fOut.Close()

def makeNumpyArray(inputcleaner, max_elems):
    import numpy as np

    first = [rr for rr in inputcleaner.getGen(1)][0]
    jetType = makeListOfStructVariables(first)
    #print jetType

    print "This 25000 is a hack for MemoryError"
    ra = np.zeros(min(inputcleaner.size(), max_elems if max_elems >= 0 else 25000), dtype=jetType)
    for iRec, rec in enumerate(inputcleaner.getGen(max_elems)):
        if iRec >= 25000:
            print "Limiting to 25000 jets per file"
            break
        for kk, vv in rec.iteritems():
            if kk in _indexToNames:
                namelist = [li[0] for li in _indexToNames[kk]]
                npvv = np.array(vv)
                for iname, name in enumerate(namelist):
                    ra[iRec][kk][name][:len(vv)] = npvv[:,iname]
                ra[iRec]['n_%s' % kk] = len(vv)
                continue

#                print npvv
#                print npvv.shape
#                print ra[iRec][kk]
#                print ra[iRec][kk]['E']
#                sys.exit(0)
#                for ielem, elem in enumerate(vv):
#                    ra[iRec][kk][ielem]
                #dType = _indexToNames[kk]
                #sublist = np.zeros(len(rec[kk]), dtype=dType)
                #for isubsublist, subsublist in enumerate(vv):
                #    for iElem, elem in enumerate(subsublist):
                #        sublist[isubsublist][_indexToNames[kk][iElem][0]] = elem
                #ra[iRec][kk] = sublist
                #continue

            if kk in _stringToInts:
                ra[iRec][kk] = _stringToInts[kk][vv]
                continue

            if isinstance(vv, tuple) or isinstance(vv, list):
                ra[iRec][kk][:len(vv)] = np.array(vv)
                continue

            ra[iRec][kk] = vv

    ra = ra[:iRec+1]
    return ra

def makeNumpyFile(outFile, inputcleaner, max_elems):
    import numpy as np

    ra = makeNumpyArray(inputcleaner, max_elems)
    fhandle = open(outFile, 'w')
    np.save(fhandle, ra)
    fhandle.close()

if __name__ == '__main__':
    opts, args = handleCLI()

    incoming = ic.InputCleanerGen(args, sys.stdin, search=opts.search,
                                                   repl  =opts.repl)

    if opts.format == 'root':
        makeRootFile (opts.output, incoming, opts.end)
    else:
        makeNumpyFile(opts.output, incoming, opts.end)

### Testing suite just in case
#    import numpy as np
#    ra = np.load(open(opts.output))
#    print
#    print ra.shape
#    print ra
#    print ra['npv'], ra['e']
#    print "ra['rel_subjets'][0]['E']"
#    print  ra['rel_subjets'][0]['E']
#    print "ra[0]['rel_subjets']['E']"
#    print  ra[0]['rel_subjets']['E']
