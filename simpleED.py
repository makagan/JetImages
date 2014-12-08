import sys
import math
import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.ticker    as tck 
import matplotlib.gridspec  as grdspc
import InputCleaner         as ic

PI = math.pi

def makeNCaloImages(pdfpages, listOfRecords):
    iRecord = 0
    key = 'cells'

    nRows = 4
    nCols = 4
    wRatio = 1
    combinedFig = plt.figure()
    combinedFig.suptitle('$<\\mu>=0')

    for record in listOfRecords:
        isRight  = ((iRecord+1)%nCols == 0)
        isLeft   = ((iRecord  )%nCols == 0)
        isTop    = ((iRecord/nRows) == 0)
        isBottom = ((iRecord/nRows) >= nRows - 1)

        #cells = sorted(record[key], key=lambda pp:
        #pts  = np.array([pp[0] for pp in record[key]])
        #etas = np.array([pp[1] for pp in record[key]])
        #phis = np.array([pp[2] for pp in record[key]])
        pts = np.array([pp[0] for pp in record[key]])

        n_on_side = math.sqrt(pts.shape[0])
        pts = pts.reshape(n_on_side, n_on_side)

        print pts
        print pts.shape
        ax = combinedFig.add_subplot(nRows, nCols, iRecord+1)

        print 'Im showing'
        ax.imshow(pts, interpolation='nearest')
        print 'Im showing'
        #ax.contourf(etas, phis, pts)

        if isLeft:
            ax.set_ylabel(r'$\phi$')
        else:
            ax.set_yticklabels([])

        if isBottom:
            ax.set_xlabel(r'$Cos \theta$')
        else:
            ax.set_xticklabels([])

        iRecord += 1
        if iRecord == nCols * nRows:
            break

    plt.savefig(pdfpages, format='pdf')

def makeNImages(pdfpages, listOfRecords, title, testFunction=None, doBoost=True):
    iRecord = 0
    key = ('bParts' if doBoost else 'parts')
    key = 'cells'

    nRows = 4
    nCols = 4
    wRatio = 1
    combinedFig = plt.figure()
    combinedFig.suptitle(('Unboosted ' if doBoost else 'Lab Frame ') + title)
    #gs_comb = grdspc.GridSpec(nRows,nCols,wspace=0.08,width_ratios=wRatio)
#    for ii, record in enumerate(listOfRecords):
#        print record['eventNumber']
#        print ['%.2f' % record[key][ii][2] for ii in range(record['nParts'])]
#        if ii >= 10: break
    for record in listOfRecords:

        isRight  = ((iRecord+1)%nCols == 0)
        isLeft   = ((iRecord  )%nCols == 0)
        isTop    = ((iRecord/nRows) == 0)
        isBottom = ((iRecord/nRows) >= nRows - 1)

        particles = sorted(record[key], key=lambda li:li[1], reverse=False)
        ids  = [pp[0] for pp in particles] 
        # if boosted, then these are really energies
        pts  = [pp[1] for pp in particles]
        # if boosted, then these are really cosThetas
        etas = [pp[2] for pp in particles]
        phis = [pp[3] for pp in particles]

#        ids  = [record[key][ii][0] for ii in range(record['nParts'])]
#        # if boosted, then these are really energies
#        pts  = [record[key][ii][1] for ii in range(record['nParts'])]
#        # if boosted, then these are really cosThetas
#        etas = [record[key][ii][2] for ii in range(record['nParts'])]
#        phis = [record[key][ii][3] for ii in range(record['nParts'])]

        ids   = np.array(ids )
        pts   = np.array(pts ) #/ (record['m'] if doBoost else record['pt'])
        etas  = np.array(etas) 
        phis  = np.array(phis) 

        pts = pts/np.cosh(etas)
        pts = np.minimum(pts, 10*np.ones(*pts.shape))
        ylims = (0, 6.29) #purposefuly a smidge larger than 2*pi
        xlims = (-2.5, 2.5)

        if not doBoost:
            #etas -= record['eta']
            #phis -= record['phi']
            #phis = np.where(phis > PI, phis-2*PI, phis)
            #phis = np.where(phis <-PI, phis+2*PI, phis)
            phis  = np.where(phis < 0, 2*PI+phis, phis)

            #etas *= 10
            #phis *= 10
            #ylims = (-10, 10) 
            #xlims = (-10, 10)
            pass

        aids = abs(ids)
        #charged pi/K/p
        hads = np.any([aids == 211, aids==321, aids==2212], axis=0)
        #charged e/mu
        leps = np.any([aids==11, aids==13], axis=0)

        cells = (aids == 0)

        markers = [
                   (cells    , 'b'      , 's',  0.1, 'Cells'),
                   (hads     , 'b'      , 'o', 0.05, 'Hadrons'),
                   #(hads     , 'b'      , 'o', 'Hadrons'),
                   #(ids == 22, '#FF9900', 's', 'Photons'),
                   #(leps     , '#A80000', 's', 'Leptons'),
                  ]

        ax = combinedFig.add_subplot(nRows, nCols, iRecord+1)
        lines  = []
        labels = []
        if False:#not doBoost:
            c1 =plt.Circle((record['eta'],record['phi']),record['radius'],
                                color='r',fill=False, zorder=100)
            c2 =plt.Circle((record['eta'],record['phi']+2*PI),record['radius'],
                                color='r',fill=False, zorder=100)
            c3 =plt.Circle((record['eta'],record['phi']-2*PI),record['radius'],
                                color='r',fill=False, zorder=10)
            plt.gcf().gca().add_artist(c1)
            plt.gcf().gca().add_artist(c2)
            plt.gcf().gca().add_artist(c3)

        for sliceObj, cc, mm, userLinSize, label in markers:

            if 0 == etas[sliceObj].shape[0]:
                continue


            ax.set_xlim(*xlims)
            ax.set_ylim(*ylims)
            ax.get_yaxis().set_major_locator(tck.MaxNLocator(nbins=5, integer=True))
            ax.get_xaxis().set_major_locator(tck.MaxNLocator(nbins=5, integer=True))

#            line = ax.scatter(etas[sliceObj], phis[sliceObj], c=pts[sliceObj],
#                            label=label, edgecolors='none', s=1, marker=mm, alpha=1.0)

            origTest = (xlims[0], ylims[0])
            dxTest   = [xx + userLinSize for xx in origTest] 
            dsp = ax.transData.transform([origTest, dxTest])
            dataSize = (dsp[0][0] - dsp[1][0])*(dsp[0][1]-dsp[1][1])
            print userLinSize, dsp, dataSize
            #dataSize *= 1000
            line = ax.scatter(etas[sliceObj], phis[sliceObj], c=pts[sliceObj],
                            label=label, edgecolors='none', s=dataSize, 
                            marker=mm, alpha=1.0, zorder=50)

            if isLeft:
                if doBoost:
                    ax.set_ylabel(r'$\phi$')
                else:
                    ax.set_ylabel(r'$\Delta \phi \times \mathit{10}$')
            else:
                #ax.tick_params(axis='y', length=0, width=0)
                ax.set_yticklabels([])

            if isBottom:
                if doBoost:
                    ax.set_xlabel(r'$Cos \theta$')
                else:
                    ax.set_xlabel(r'$\Delta \eta \times \mathit{10}$')
            else:
                #ax.tick_params(axis='x', length=0, width=0)
                ax.set_xticklabels([])

            lines.append(line)
            labels.append(label)

        #fig.legend(lines, labels)

        iRecord += 1
        if iRecord == nCols * nRows:
            break

    plt.savefig(pdfpages, format='pdf')

def makeGroups(fnames, rangeA=None, rangeB=None):
    retDi = {'q': {}, 'g':{}, 'w': {}}
    retDi = {'q': [], 'g':[], 'w': []}

    if rangeA is None:
        low, high = None, None
    elif rangeB is None:
        low, high = 0, rangeA
    else:
        low, high = rangeA, rangeB

#    for kk, di in retDi.iteritems():
#        for ii in range(low, high):
#            di[ii] = []

    for fname in fnames:
        records = np.load(open(fname))
        for rec in records:
            flavor = 'w' if rec['ME'] == 2 else ('g' if rec['pdgIDHardParton'] == 21 else 'q')
            nParts = rec['nParts']
            retDi[flavor].append(rec)
            #if high > nParts >= low:
            #    retDi[flavor][nParts].append(rec)

    return retDi

def powerSpectrum(records, lmax=40):
    import scipy.special as spe
    print '~'*50
    print records
    print '~'*50

    key = 'bParts'
    power = np.zeros(shape=(lmax, len(records)))
    for iRecord, rec in enumerate(records):
        es    = [rec[key][ii][1] for ii in range(rec['nParts'])]
        cosTs = [rec[key][ii][2] for ii in range(rec['nParts'])]
        phis  = [rec[key][ii][3] for ii in range(rec['nParts'])]
    
        thetas = np.arccos(cosTs)
        phis   = np.array(phis)
        es     = np.array(es)
        #es    /= np.sum(es*es)
        es    /= math.sqrt(np.sum(es*es))

        for ll in range(lmax):
            for mm in range(-ll, ll+1):
    #note scipy takes opposite unit names from Physics
    #phi (in scipy) is angle from z-axis
                flm = np.sum(                 spe.sph_harm(mm, ll, phis, thetas) * es)
                #flm = np.sum(np.sin(thetas) * spe.sph_harm(mm, ll, phis, thetas) * es)
                power[ll,iRecord] += np.abs(flm + np.conj(flm))



    power = np.sqrt(power)
    return power

class FourVector:
    def __init__(self, inputs, iType='pt'):#pt, eta, phi, m):
        if iType == 'pt':
            self.pt, self.eta, self.phi, self.m = inputs
            self.initializePx()
        elif iType == 'px':
            self.px, self.py, self.pz, self.E = inputs
            self.initializePt()

    def initializePx(self):
        from math import sqrt, sin, cos, cosh
        self.E  = self.pt * cosh(self.eta)
        self.px = self.pt * cos(self.phi)
        self.py = self.pt * sin(self.phi)
        self.pz = sqrt((E**2 - self.m**2 - self.pt**2))

    def initializePt(self):
        from math import sqrt, atan2, tan, log

        self.pt  = sqrt(self.px**2 + self.py**2)
        self.eta = -1 * log(tan(atan2(self.pt, self.pz)/2))
        self.phi = atan2(self.py, self.px)
        self.m   = sqrt(self.E**2 - self.pz**2 - self.pt**2)
        print self.pt, self.eta, self.phi, self.m

    def mass(self):
        return self.m

    def p(self):
        from math import sqrt
        return sqrt(self.pt**2 + self.pz**2)

    def pt(self):
        return self.pt

    def __add__(self, other):
        sumVec = [getattr(self, dd) + getattr(other, dd) for dd in ['px', 'py', 'pz', 'E']]
        print 'Sum pxpypzE is:', sumVec
        return self.__class__(sumVec, iType='px')

def pairCorr(allRecords):
    from itertools import combinations
    import ROOT as rt

    pdfs = {}
    names = [('m', 'f8'), ('pt', 'f8'), ('dr', 'f8'),
             ('phi', 'f8'), ('cosT', 'f8'), ('p', 'f8'),
             ('nParts','i4'), ('isQCD', 'i4')]

    records = allRecords[allRecords['nParts'] == 20]
    nPartss = records['nParts']
    nPairs  = np.sum((nPartss * (nPartss -1 ))/2)
    pairVals = np.zeros(nPairs, dtype=names)
    iTrk = 0
    for iRec, rec in enumerate(records):
#        if iRec >= 1000:
#            break
        aids = rec['parts']['pdgId']
        #charged = np.any([aids == 211, aids==321, aids==2212, 
        #                  aids==11, aids==13] , axis=0)
        #for trackPair in combinations(rec['parts'][charged], 2):
        for trackPair in combinations(rec['parts'], 2):
            vecLi = []
            for track in trackPair:
                mass = 0
                vec = rt.TLorentzVector()
                vec.SetPtEtaPhiM(track[1], track[2], track[3], mass)
                vecLi.append(vec)
    
            combVec = vecLi[0] + vecLi[1]
#            bVec    = rt.TLorentzVector(vecLi[0])
#            bVec.Boost( -combVec.BoostVector())

            pairVals[iTrk]['isQCD'] = rec['pdgIDHardParton'] != 24
            pairVals[iTrk]['nParts']  = rec['nParts']
            pairVals[iTrk]['m' ]   = combVec.M()
            pairVals[iTrk]['pt']   = combVec.Pt()/rec['pt']
            pairVals[iTrk]['p' ]   = combVec.P()/(rec['pt']*math.cosh(rec['eta']))
            pairVals[iTrk]['dr']   = vecLi[0].DeltaR(vecLi[1])
#            pairVals[iTrk]['phi']  = math.atan2(bVec.Py(),bVec.Px())
#            pairVals[iTrk]['cosT'] = bVec.Pz()/math.sqrt(bVec.Pt()**2 + bVec.Pz()**2)

            if np.any([np.isnan(pairVals[iTrk]['phi']), \
                np.isnan(pairVals[iTrk]['cosT'])]):
                print 'Found NAN... dropping whole pair in', iRec
                iTrk -= 1

    
            iTrk += 1
    pairVals = pairVals[pairVals['nParts'] != 0]
    return pairVals

def drawCorr(pairVals):
    import scipy.stats as sta

    nBins = 20
    edges = {'m': [], 'dr': [], 'pt': []}
    for kk in edges:
        edges[kk] = [sta.scoreatpercentile(pairVals[kk], ii * 100.0/nBins)
                    for ii in range(nBins + 1)]

    
    cut = pairVals['isGluon'] == 1
    H, xedges, yedges = np.histogram2d( pairVals['m'][cut], pairVals[cut]['dr'],
                                        bins=[edges['m'], edges['dr']]
                                      )
    preAspect = (yedges[-1] - yedges[0])/(xedges[-1] - xedges[0])
    aspect = 1.2/preAspect

    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    plt.imshow(H, extent=extent, aspect=aspect, interpolation='nearest', origin='lower')
    plt.xticks(xedges)
    plt.yticks(yedges)
    plt.show()
    
def scaledTicks(edges):

    nBins = len(edges) - 1
    rng = edges[-1] - edges[0]
    low = edges[0]
    locs  = [float(ii)*rng/nBins  + low for ii in range(nBins + 1)]
    names = ['%.2f' % xx for xx in edges]

    return locs, names

def corr_stuff():
#    records = np.load(open(sys.argv[1]))
    pp = PdfPages('pairCorr.pdf')
    wrec = np.load(open('/u/eb/joshgc/mynfs/wjets.npy'))
    grec = np.load(open('/u/eb/joshgc/mynfs/gjets.npy'))

    wPairVals = pairCorr(wrec)
    gPairVals = pairCorr(grec)
    
    #drawCorr(pairVals)

    import scipy.stats as sta

    nBins2 = 20
    nBins1 = 100
    names = {'m'   : 'Pair Mass (GeV)', 
             'dr'  : 'Pair $\\Delta R$',
             'pt'  : 'Pair Sum Fractional $p_{T}$',
             'p'   : 'Pair Sum Fractional $p$',
             'phi' : 'Decay Angle $\\phi$',
             'cosT': 'Decay Angle $Cos \\theta$',
            }
    edges2 = dict((kk, []) for kk in names.keys())
    edges1 = dict((kk, []) for kk in names.keys())
    drawsToDo =[ ('m', 'dr'), ('m', 'pt'), ('pt', 'dr'),
                ('m', 'p'), ('p', 'dr'),
                #('m','phi'), ('m', 'cosT'), ('phi', 'cosT'),
                ]

    for kk in edges1:
        edges2[kk] = [sta.scoreatpercentile(wPairVals[kk], ii * 100.0/nBins2)
                    for ii in range(nBins2 + 1)]
        edges1[kk] = [edges2[kk][0] + (edges2[kk][-1]-edges2[kk][0])*ii/nBins1
                        for ii in range(nBins1+1)]

    print 'Edges are', edges2
    print 'Edges are', edges1
    for xKey in ['m', 'dr', 'pt', 'p']:
        Hw, xedges = np.histogram(wPairVals[xKey], bins=edges1[xKey])
        Hg, xedges = np.histogram(gPairVals[xKey], bins=edges1[xKey])
        Hw = np.maximum(Hw, .1)
        Hg = np.maximum(Hg, .1)
        Hw /= np.sum(Hw)
        Hg /= np.sum(Hg)
        bins = (xedges[:-1] + xedges[1:])/2
        plt.semilogy(bins, Hw, 'b', nonposy='clip', linewidth=5, label='W-Jets')
        plt.semilogy(bins, Hg, 'r', nonposy='clip', linewidth=5, label='G-Jets')
        plt.legend()
        plt.ylabel('Probability')
        plt.xlabel(names[xKey])
        plt.savefig(pp, format='pdf')
        plt.clf()


    for xKey, yKey in drawsToDo:
        pdfs = {}
        H, xedges, yedges = np.histogram2d( wPairVals[xKey], wPairVals[yKey],
                                            bins=[edges2[xKey], edges2[yKey]]
                                          )
        yTickLoc, yTickNames   = scaledTicks(yedges)
        xTickLoc, xTickNames   = scaledTicks(xedges)

        for doQCD in [0, 1]:
            pv = gPairVals if doQCD else wPairVals

            H, xedges, yedges = np.histogram2d( pv[xKey], pv[yKey],
                                                bins=[edges2[xKey], edges2[yKey]],
                                                )
            H /= np.sum(H)
            preAspect = (xedges[-1] - xedges[0])/(yedges[-1] - yedges[0])
            aspect = preAspect/1.2
        
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            pdfs[doQCD] = H
            plt.imshow(H.T, extent=extent, aspect=aspect, interpolation='nearest', origin='lower')

            plt.xticks( xTickLoc, xTickNames, rotation=45)
            plt.yticks( yTickLoc, yTickNames)
            plt.ylabel(names[yKey])
            plt.xlabel(names[xKey])
            plt.colorbar()

            plt.title('Q/G' if doQCD == 1 else 'W')
            plt.savefig(pp, format='pdf')
            plt.clf()

        plt.imshow((pdfs[1] - pdfs[0]).T, extent=extent, aspect=aspect, 
                    interpolation='nearest', origin='lower')
        plt.xticks( xTickLoc, xTickNames, rotation=45)
        plt.yticks( yTickLoc, yTickNames)
        plt.ylabel(names[yKey])
        plt.xlabel(names[xKey])
        plt.colorbar()
        plt.title('P(g/q) - P(W)')
        plt.savefig(pp, format='pdf')
        plt.clf()

        plt.imshow(((pdfs[1] - pdfs[0])/(pdfs[1]+pdfs[0])).T, extent=extent, aspect=aspect, 
                    interpolation='nearest', origin='lower')
        plt.xticks( xTickLoc, xTickNames, rotation=45)
        plt.yticks( yTickLoc, yTickNames)
        plt.ylabel(names[yKey])
        plt.xlabel(names[xKey])
        plt.colorbar()
        plt.title('P(g/q) - P(W) / P(g/q) + P(W)')
        plt.savefig(pp, format='pdf')
        plt.clf()

        plt.imshow((pdfs[0]/pdfs[1]).T, extent=extent, aspect=aspect, 
                    interpolation='nearest', origin='lower')
        plt.xticks( xTickLoc, xTickNames, rotation=45)
        plt.yticks( yTickLoc, yTickNames)
        plt.ylabel(names[yKey])
        plt.xlabel(names[xKey])
        plt.colorbar()
        plt.title('P(W)/P(g/q)')
        plt.savefig(pp, format='pdf')
        plt.clf()

    pp.close()
    sys.exit(0)

if __name__ == '__main__':
    from matplotlib.backends.backend_pdf import PdfPages

    incoming = ic.InputCleanerGen(sys.argv[1:], sys.stdin, search=r'^#G(.*)$', repl=r'\1')
#    low, high = 25, 26
#    print 'Beginning group production'
#    recGroups = makeGroups(sys.argv[1:], low, high)
#
    pp = PdfPages('ed.caloJets.pdf')
    makeNCaloImages(pp, [rr for rr in incoming.getGen(15)])
#    cuts = []
#    fTrue = lambda xx: True
#    for ii in range(low, high):
#
#        cuts.append(['W-Jets with %d Particles'     % ii, 'w', ii, None])
#        cuts.append(['Gluon Jets with %d Particles' % ii, 'g', ii, None])
#
#    for title, flavor, nParts, cutFunc in cuts:
#        print 'Working on', title
##        print len(recGroups[flavor][nParts])
##        makeNImages(pp, recGroups[flavor][nParts], title, None, doBoost=False)
#        print len(recGroups[flavor])
#        makeNImages(pp, recGroups[flavor], title, None, doBoost=False)

    pp.close()
