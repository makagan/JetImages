#!/usr/bin/env python
import sys
import math
import numpy                as np
import numpy.linalg as linalg
import matplotlib           as mpl
import matplotlib.ticker    as ticker
from matplotlib.colors import LogNorm
import matplotlib.pyplot    as plt
import matplotlib.cm        as cm
import scipy.interpolate    as inp
from   scipy.ndimage.interpolation import geometric_transform
from matplotlib import rc
rc('font', **{'family':'serif'})#, 'serif':['Times']})
rc('text', usetex=True)
rc('axes', linewidth=2)
rc('xtick.major', width=2)
rc('ytick.major', width=2)

np.set_printoptions(precision=3, suppress=True)

def handleCLI():
    from optparse import OptionParser
    usage = '%prog [opts] -i myInput.txt -o myOuput.{npy,pdf}'
    pp = OptionParser(usage=usage)
    pp.add_option('-n', '--njets', dest='njets', type=int, default=-1,
                    help='Specify number of jets (or images) to rotate.'
                        ' Specifying -1 runs on the entire sample.'
                        ' Currently the run is tacitly truncated at 10K')
    pp.add_option('-j', '--subjets', dest='subjets', action='store_true',
                    help='Use subjets to find the transformation parameters,'
                        ' otherwise, the default uses prinicple axes')

    pp.add_option('-t', '--target', dest='target', type=float, default=None,
                    help='Specify the target (along the axis after rotation)'
                        ' of the leading subjet.  Defaults to None when it is ignored')

    pp.add_option('-l', '--length', dest='length', type=float, default=None,
                    help='Specify the distance between the two subjets as a'
                        ' fraction of the total image length.  -1 to ignore.')

    pp.add_option('-i', '--input', dest='input', type=str, default=None,
                    help='Specify input file of unrotated text-jets')

    pp.add_option('-o', '--output', dest='output', type=str, default=None,
                    help='Specify output file. If the extension is .npy'
                        ' then a numpy file is of rotated jets is created.'
                        ' If a .pdf extension is used then some plots are made'
                        ' instead.')

#    pp.add_option('-s', '--static', dest='static', action='store_true', 
#                    default=False,
#                    help='Do no rotation and put on plot on a page')

    opts, args = pp.parse_args()
    if None in (opts.output, opts.input):
        pp.print_usage()
        sys.exit(1)

    #print "Rotating Jets"
    #print opts, args
    return opts, args

def deltaR(aa, bb):
    dphi = abs(aa[1] - bb[1])
    dphi = (2*math.pi - dphi) if dphi > math.pi else dphi
    return math.sqrt((aa[0] - bb[0])**2 + dphi**2)


def getPrincipleAxes(theMat):

    if theMat.shape != (2,2):
        print "ERROR: getPrincipleAxes(theMat), theMat size is not 2x2. DYING!"
        sys.exit(1)

    laS, laV = linalg.eigh(theMat)
    return -1*laV[::-1], laS[::-1]

def makeXYScatterMatrix(record, pow_forScat=1, pow_forMean=1):
    jet_rad = record['radius']
    #eta_off = record['eta'] - record['jet_offset'][0]
    #phi_off = record['phi'] - record['jet_offset'][1]

    cell_values = np.array([cell[0]   for cell in record['cells']])
    cell_x      = np.array([cell[1]   for cell in record['cells']])
    cell_y      = np.array([cell[2]   for cell in record['cells']])
    cell_values = np.array(cell_values)


    eta_off = (np.max(cell_x) + np.min(cell_x))/2
    phi_off = (np.max(cell_y) + np.min(cell_y))/2

    mask = ( (np.square(cell_x - eta_off) 
            + np.square(cell_y - phi_off)) < jet_rad**2)

    etot = np.sum((cell_values>0) * np.power(cell_values, pow_forMean))
    if etot == 0:
        print 'Found a jet with no energy.  DYING!'
        sys.exit(1)

    x_1  = np.sum((cell_values>0) * np.power(cell_values, pow_forMean) * cell_x)/etot
    y_1  = np.sum((cell_values>0) * np.power(cell_values, pow_forMean) * cell_y)/etot
    x_2  = np.sum(mask * (cell_values>0) * np.power(cell_values, pow_forScat) * np.square(cell_x -x_1))
    y_2  = np.sum(mask * (cell_values>0) * np.power(cell_values, pow_forScat) * np.square(cell_y -y_1))
    xy   = np.sum(mask * (cell_values>0) * np.power(cell_values, pow_forScat) * (cell_x - x_1)*(cell_y -y_1))

    ScatM = np.array([[x_2, xy],[xy, y_2]])
    MeanV = np.array([x_1, y_1])

    return ScatM, MeanV

def getMaxPartialSum(record, n_side_partial=2):
    cell_values = np.array([cell[0]   for cell in record['cells']])
    cell_x      = np.array([cell[1]   for cell in record['cells']])
    cell_y      = np.array([cell[2]   for cell in record['cells']])
    cell_values = np.array(cell_values)

    n_side = int(np.sqrt(cell_x.shape[0]))
    cell_values = cell_values.reshape(n_side, n_side)
    cell_x = cell_x.reshape(n_side, n_side)
    cell_y = cell_y.reshape(n_side, n_side)


    max_partial_sum = -1
    max_x_y = np.array([-1.0,-1.0])
    for xbin in range(n_side - n_side_partial):
        for ybin in range(n_side - n_side_partial):
            partial_sum = np.sum(cell_values[ybin:ybin+2, xbin:xbin+2])

            if partial_sum <= 0:
                continue

            if partial_sum > max_partial_sum:
                max_partial_sum = partial_sum
                max_x_y[0] = np.sum( cell_x[ybin:ybin+2, xbin:xbin+2] * cell_values[ybin:ybin+2, xbin:xbin+2] ) / partial_sum
                max_x_y[1] = np.sum( cell_y[ybin:ybin+2, xbin:xbin+2] * cell_values[ybin:ybin+2, xbin:xbin+2] ) / partial_sum
               
    return max_partial_sum, max_x_y


def etaPhiToCellCoords(ptEtaPhiTriplet):#, startLocation):

    return (ptEtaPhiTriplet['eta'], ptEtaPhiTriplet['phi'])
    return ptEtaPhiTriplet[1:3]
    dEta = ptEtaPhiTriplet[1] - startLocation[0]
    dPhi = ptEtaPhiTriplet[2] - startLocation[1]

    return dEta, dPhi

#def geoBoostCells(record, do_principles, length_frac, y_target, use_pt):
#
#    sideLength = int(math.sqrt(record['n_cells']))
#
#    cell_x      = np.array(record['cells']['eta'])
#    cell_y      = np.array(record['cells']['phi'])
#    cell_values = np.array(record['cells']['E'  ])
#    if use_pt is True:
#        cell_values /= np.cosh(cell_x)
#    cell_points = np.column_stack((cell_x, cell_y))


def geoTransCells(record, do_principles, length_frac, y_target, ):

    sideLength = int(math.sqrt(record['n_cells']))

    cell_x      = np.array(record['cells']['eta'])
    cell_y      = np.array(record['cells']['phi'])
    cell_values = np.array(record['cells']['E'  ])
    cell_points = np.column_stack((cell_x, cell_y))

    ll=1.0
    ctheta=0.0
    stheta=0.0
    centerX=0.0
    centerY=0.0
    imageLength = np.max(cell_y) - np.min(cell_y)
    partial_sum = 0
    partial_coords = np.array([-1.0,-1.0])

    #we know we have at least 1 subjet
    bigJet   = etaPhiToCellCoords( record['rel_subjets'][0] )#, record['jet_offset'])
    smallJet = etaPhiToCellCoords( record['rel_subjets'][1] )#, record['jet_offset'])
    jet = [record['eta']-record['jet_offset'][0], record['phi']-record['jet_offset'][1]]
    if jet[1] < 0:
        jet[1] += 2 * math.pi

    
    if do_principles:
        scat, cent = makeXYScatterMatrix(record)
        paxes, pvars = getPrincipleAxes(scat)
        ll = pvars[0]/sideLength

        partial_sum, partial_coords = getMaxPartialSum(record)
        #check sign of rotation is correct for lead jet on top
        paxes[0,:] = paxes[0,:] * np.sign( np.dot(paxes[0,:], bigJet[0:2] - cent[0:2]) ) 

        stheta = paxes[0,0]
        ctheta = paxes[0,1]
        centerX = cent[0]
        centerY = cent[1]
    else:
        #in this case, we know we have second subjet

        ll = math.sqrt((bigJet[0] - smallJet[0])**2 +(bigJet[1] - smallJet[1])**2)
        ctheta = (bigJet[1] - smallJet[1])/ll
        stheta = (bigJet[0] - smallJet[0])/ll
        #centerX = (record['rel_subjets'][0]['E']*bigJet[0] + record['rel_subjets'][1]['E']*smallJet[0])/np.sum(record['rel_subjets'][:2]['E'])
        #centerY = (record['rel_subjets'][0]['E']*bigJet[1] + record['rel_subjets'][1]['E']*smallJet[1])/np.sum(record['rel_subjets'][:2]['E'])
        #print centerX, centerY
        centerX = (bigJet[0] + smallJet[0])/2
        centerY = (bigJet[1] + smallJet[1])/2
        #print centerX, centerY

    wCenterX = (np.max(cell_x) + np.min(cell_x))/2
    wCenterY = (np.max(cell_y) + np.min(cell_y))/2

    #wLength = sideLength*0.5
    #wLength = ll * 0.2 / deltaR(bigJet, jet)
    #print wLength/ll

    if length_frac > 0:
        wLength = length_frac * imageLength
    else:
        wLength = ll

    if y_target is not None:
        if do_principles:
            wCenterY = y_target + wLength/2
        else:
            wCenterY = y_target - wLength/2
    

    R  = np.matrix([[ctheta, -stheta],[stheta,ctheta]])
    S  = np.matrix([[1,0],[0,wLength/ll]])
    CH = np.matrix([[centerX],[centerY]])
    CW = np.matrix([[wCenterX],[wCenterY]])
    try:
        Rinv = linalg.inv(R)
        Sinv = linalg.inv(S)
    except np.linalg.linalg.LinAlgError, e:
        print 'Failed on eventNumber %d with seed %d' % (record['eventNumber'], record['seed'])
        raise e
        


    def shift_func(coords):
        # used in passive rotations, which is used for rotating the grid
        C = np.matrix(coords).getT()

        A = (Rinv * (Sinv * (C - CW)) + CH)
        A = A.getA1()
        return A

    def shift_func_inv(coords):
        # used in active rotations, which is used for rotating the individual points for display
        # and used for calculating the offsets when targets are specified
        C = np.matrix(coords).getT()

        A = S * (R* (C - CH)) + CW
        A = A.getA1()
        return A

    test_x = np.linspace(min(cell_x), max(cell_x), sideLength)
    test_y = np.linspace(min(cell_y), max(cell_y), sideLength)

    test_grid = np.zeros((test_x.shape[0], test_y.shape[0], 2))

    orig_grid_x, orig_grid_y = np.meshgrid(test_x, test_y)
    test_grid[:,:,0] = orig_grid_x
    test_grid[:,:,1] = orig_grid_y
    test_points = test_grid.reshape((test_grid.shape[0] * test_grid.shape[1], 2))
    test_points = np.array([shift_func(point) for point in test_points])

    test_values = inp.griddata(cell_points, cell_values, test_points, fill_value=0)

    orig_values = listToGrid(cell_values)
    rot_values  = listToGrid(test_values)
    return orig_grid_x, orig_grid_y,                \
           orig_values, rot_values ,                \
           np.array(shift_func_inv(bigJet  [0:2])), \
           np.array(shift_func_inv(smallJet[0:2]))


def listToGrid(ra):
    oldShape = ra.shape
    newLen   = int(math.sqrt(oldShape[0]))
    newShape = (newLen, newLen) + ra.shape[1:]

    return ra.reshape(newShape)

def reflect(cells, axis=1):

    if axis not in [0,1]:
        print "Reflect Error: axis must be either 0 or 1!"
        return
    
    lx,ly = cells.shape

    if axis == 1:
        # flip over y-axis 
        #bounds to exclude the central stripe, needs odd number of bins
        left  = np.sum(cells[:, ly/2+1:    ])
        right = np.sum(cells[:,       :ly/2])

        if left > right:
            return cells, False

        return cells[:,::-1], True

    # if here, axis must be 0
    # flip over x-axis
    #bounds to exclude the central stripe, needs odd number of bins
    top =    np.sum(cells[lx/2+1:    ,:])
    bottom = np.sum(cells[      :lx/2,:])
    print top, bottom

    if top > bottom:
        return cells, False

    return cells[::-1,:], True

def plotStuff(pp, orig_values, rot_values, orig_subjets, rot_subjets):

    fig  = plt.figure(figsize=(5, 5))
    fig.subplots_adjust(top=.97, right=.95, 
                        left=0.13, bottom=.13)
    size = .80
    vmax = 2 * 10**3
    vmin = 5e-4

    ll = int(orig_values.shape[0])
    low = -.1
    high = ll/10.0 + .1

    for do_rot in (False, True, ):
        subjets = rot_subjets if do_rot else orig_subjets
        for do_legend in (False, True, ):
            fig.clf()

            ax = fig.add_axes((0.15, 0.12, size, size))
            ret = ax.imshow(orig_values if not do_rot else rot_values,
                            cmap=cm.jet,
                            norm=LogNorm(),
                            interpolation='nearest',
                            origin='lower',
                            extent=[low, high, low, high],
                            vmin=vmin,
                            vmax=vmax,
                            )
            if subjets is not None and len(subjets) > 1:
                for ii, li in enumerate([('Leading Subjet', 'o'), 
                                         ('2$^{\\mathbf{nd}}$ Subjet', 'v')]):

                    ax.scatter( subjets[ii][0], subjets[ii][1],
                                s=200, edgecolor='#a61700', facecolor='None', 
                                linewidths=3, marker=li[1],label=li[0])
            if do_rot:
                ax.set_ylabel('Q$_{\\mathbf{1}}$', fontsize='x-large', labelpad=9)
                ax.set_xlabel('Q$_{\\mathbf{2}}$', fontsize='x-large', labelpad=12)
            else:
                ax.set_ylabel('Relative $\\mathbf{\phi}$', fontsize='x-large', labelpad=9)
                ax.set_xlabel('Relative $\\mathbf{\eta}$', fontsize='x-large', labelpad=12)

            _ticks = np.linspace(0, ll/10.0, 6)
            ax_fmt = ticker.FormatStrFormatter('$\\mathbf{%s}$')
            ax.tick_params(axis='both', labelsize='large')
            ax.xaxis.set_major_formatter(ax_fmt)
            ax.yaxis.set_major_formatter(ax_fmt)
            ax.set_xticks(_ticks)
            ax.set_yticks(_ticks)

            if do_legend:
                ax_cbar = fig.add_axes((0.78,0.25,0.05,0.5))
                seq = ['$\\mathbf{10^{%d}}$' % ii for ii in range(
                        int(math.ceil(math.log10(vmin))),
                        int(math.log10(vmax)) + 1)
                        ]
                cb_fmt = ticker.FixedFormatter(seq)
                ax_cbar.set_xlabel('$\\mathbf{p_T}$ [GeV]', labelpad=9)
                ax_cbar.tick_params(axis='both', labelsize='large')
                cbar = fig.colorbar(ret, ax_cbar, format=cb_fmt)
                leg = ax.legend(loc='upper right', frameon=True, scatterpoints=1,
                                #fontsize='large')
                                )
                leg.get_frame().set_edgecolor('white')
                leg.get_frame().set_facecolor('white')

            plt.savefig(pp, format='pdf')


def rotationStuff(in_name, out_name, n_jets, do_principles, y_target, length_frac, Use_Y_Reflection = True):

    if out_name.endswith('pdf'):
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(out_name)
        #print 'Writing pdf'
    import InputCleaner as ic
    incoming = ic.InputCleanerGen([in_name],sys.stdin, search='^#G(.*)$')

    whichJet = 0

    import text2bin
    ra = text2bin.makeNumpyArray(incoming, n_jets)
    if n_jets == -1:
        n_jets = len(ra)
    ra = ra[:n_jets]
    out_array = np.array(ra)

    for rr in out_array:
        n_subjets = np.count_nonzero(rr['rel_subjets']['E'])
        n_cells   = np.count_nonzero(rr['cells']['E'])
        if n_subjets <1:
            continue
        if not do_principles and n_subjets <2:
            continue
        if n_cells <2:
            continue

        orig_grid_x, orig_grid_y, orig_values,\
        rot_values, shiftBigJet, shiftSmallJet = geoTransCells(rr, do_principles, length_frac, y_target, )

        #reflect over y-axis.  But not if there is a target, because then axis is already shifted and reflection can throw things out of grid
        if do_principles and Use_Y_Reflection and y_target == None:
            rot_values, isrefY = reflect(rot_values, axis=0)
            print 'reflecting y', isrefY
            if isrefY:
                shiftBigJet[1] = np.max(orig_grid_y[:,0])+np.min(orig_grid_y[:,0]) - shiftBigJet[1]
                shiftSmallJet[1] = np.max(orig_grid_y[:,0])+np.min(orig_grid_y[:,0]) - shiftSmallJet[1]


        #obligatory reflection along x axis
        rot_values, isrefX = reflect(rot_values, axis=1)
        print 'reflecting x', isrefX
        if isrefX:
            shiftBigJet[0] = np.max(orig_grid_x[0])+np.min(orig_grid_x[0]) - shiftBigJet[0]
            shiftSmallJet[0] = np.max(orig_grid_x[0])+np.min(orig_grid_x[0]) - shiftSmallJet[0]


        out_array['cells'][whichJet] = rot_values.reshape(1, rot_values.shape[0]*rot_values.shape[1])

#used to make no rotation set
#        orig_values = rr['cells']['E']
#        out_array['cells'][whichJet] = orig_values#.reshape(1, orig_values.shape[0]*orig_values.shape[1])

        whichJet += 1
        if not out_name.endswith('pdf'):
            continue

        orig_subjets = [(rr['rel_subjets'][0][1], rr['rel_subjets'][0][2]),
                        (rr['rel_subjets'][1][1], rr['rel_subjets'][1][2])]
        plotStuff(pp, orig_values, rot_values, orig_subjets, [shiftBigJet, shiftSmallJet])

    try:
        pp.close()
        return
    except NameError:
        pass

    out_array = out_array[out_array['m']>-1]
    fh = open(out_name, 'w')
    np.save(fh, out_array)
    fh.close()

    print 'Saved %d jets' % whichJet

    return ra, out_array


def main():
    opts, args = handleCLI()

    rotationStuff(opts.input, opts.output, opts.njets,
                    not opts.subjets, opts.target, opts.length)


if __name__ == '__main__':
    main()
