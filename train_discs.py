#!/usr/bin/env python
'''Main application for training many discriminants'''

__author__ = 'Josh Cogan (jcogan@cern.ch)'

import re
import os
import time
import numpy as np
import samples
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
rc('font', **{'family':'serif'})#, 'serif':['Times']})
rc('text', usetex=True)
rc('axes', linewidth=2)
rc('xtick.major', width=2, size=5)
rc('ytick.major', width=2, size=5)
rc('xtick.minor', width=2, size=5)
rc('ytick.minor', width=2, size=5)

_tick_fmt_int   = ticker.FormatStrFormatter('$\\mathbf{%d}$')
_tick_fmt_float = ticker.FormatStrFormatter('$\\mathbf{%s}$')

def moving_average(aa, nn):
#$    cum = np.cumsum(aa, dtype=float)
#$    div = 1.0/np.cumsum(np.ones(aa.shape, dtype=float))
#$    ret[nn-1:(cum[nn - 1:] - cum[:1 - nn]) / nn
    window = np.ones(int(nn))/float(nn)
    return np.convolve(aa, window, 'same')

def last_local_pickle():

    files = filter(os.path.isfile, os.listdir('.'))
    files = filter(lambda ff: os.path.splitext(ff)[1]=='.pickle', files)
#vetos the pickle files that contains hists
    files = filter(lambda ff: '-test-' not in ff, files)
    files.sort(key=lambda x: os.path.getmtime(x))
    try:
        return files[-1]
    except IndexError:
        return None

def handleCLI():
    from optparse import OptionParser

    pp = OptionParser()
    pp.add_option('--gen_a', dest='gen_a', default='p',
                    help='p for pythia or h for herwig for class a')
    pp.add_option('--mu_a', dest='mu_a', type='int', default=0,
                    help='Integer number of avg interactions / crossing for class a')
    pp.add_option('--me_a', dest='me_a', type='str', default='v',
                    help='ME type for class a. Same mapping as pythiaJets'
                         ' i.e. v->Wj->lnuj->, w->WW->lnujj')
    pp.add_option('--pt_a', dest='pt_a', type='int', default=200,
                    help='Pt bin for class a')
    pp.add_option('--trim_a', dest='trim_a', type='int', default=1,
                    help='Trimming on (1) or off (0)')

    pp.add_option('--smear_a', dest='smear_a', type='int', default=0,
                    help='Smearing on (1) or off (0)')    

    pp.add_option('--pdg_a',  dest='pdg_a', default=True, action='store_false',
                    help='Disable filtering on pdgid')

    pp.add_option('--key_a', dest='key_a', type='str', default=None,
                    help='This commands overrides all other'
                        ' dataset specifiers (if used) and looks for the'
                        ' exact key provided.')
    
    pp.add_option('--gen_b', dest='gen_b', default=None,
                    help='p for pythia or h for herwig for class b')
    pp.add_option('--mu_b', dest='mu_b', type='int', default=None,
                    help='Integer number of avg interactions / crossing for class b')
    pp.add_option('--me_b', dest='me_b', type='str', default=None,
                    help='ME type for class b. Same mapping as pythiaJets'
                         ' i.e. v->Wj->lnuj->, w->WW->lnujj')
    pp.add_option('--pt_b', dest='pt_b', type='int', default=None,
                    help='Pt bin for class b')
    pp.add_option('--trim_b', dest='trim_b', type='int', default=1,
                    help='Trimming on (1) or off (0)')

    pp.add_option('--smear_b', dest='smear_b', type='int', default=None,
                    help='Smearing on (1) or off (0)')    
    pp.add_option('--pdg_b',  dest='pdg_b', default=None, action='store_true',
                    help='Disable filtering on pdgid')

    pp.add_option('--key_b', dest='key_b', type='str', default=None,
                    help='This commands overrides all other'
                        ' dataset specifiers (if used) and looks for the'
                        ' exact key provided.')

    pp.add_option('-r', '--rad', dest='rad', type='float', default='1.2',
                    help='Radius of jet.  Unlike the other dataset specifiers'
                         ' this one must be identical for both classes (and'
                         ' testing!)')

    pp.add_option('-s', '--start', dest='start', type='float', default=0,
                    help='Which index to start the training slice in for input'
                        ' dataset.  If abs(start) < 1 then start*=num_entries.'
                        ' Training is performed on all_input[start:end].')
    pp.add_option('-e', '--end', dest='end', type='float', default=-1,
                    help='Which index to start the training slice in for input'
                        ' dataset.  If abs(start) < 1 then start*=num_entries.'
                        ' Training is performed on all_input[start:end].')

    pp.add_option('--min_m',  dest='min_m', type='float', default=-1,
                  help='The minimum jet mass.')
    pp.add_option('--max_m',  dest='max_m', type='float', default=-1,
                  help='The maximum jet mass.')

    pp.add_option('--min_pt',  dest='min_pt', type='float', default=-1,
                  help='The minimum jet pT.')
    pp.add_option('--max_pt',  dest='max_pt', type='float', default=-1,
                  help='The maximum jet pT.')

    pp.add_option('--min_bdrs',  dest='min_bdrs', type='float', default=-1,
                  help='The minimum jet bdrs mass.')
    pp.add_option('--max_bdrs',  dest='max_bdrs', type='float', default=-1,
                  help='The maximum jet bdrs mass.')
    
    pp.add_option('--min_dr',  dest='min_dr', type='float', default=-1,
                  help='The minimum delta R between subjets')
    pp.add_option('--max_dr',  dest='max_dr', type='float', default=-1,
                  help='The minimum delta R between subjets')

    pp.add_option('--n_b',  dest='n_b', type='int', default=-1,
                  help='The exact number of B hadrons per jet')
    pp.add_option('--subjet_e',  dest='subjet_e', type='float', default=-1,
                  help='The minimum energy found in at least one subjet')

    

    pp.add_option('-d', '--re_discs', dest='re_discs', type='str', default=None,
                    help='Specify a regex that the list of picked up '
                        'discriminant names must match.  ie "PCA.*Comp" will match'
                        ' all the individual pca componenets but not fisher.'
                        '  Only for testing cycles')

    pp.add_option('-o', '--output', dest='output', type='str', default=None,
                    help='Location of output pickle file for discriminants')
    pp.add_option('-i', '--input', dest='input', type='str', default=None,
                    help='Locations of input pickle file for discriminants')
    pp.add_option('--varbins', dest='varbins', default=False, action='store_true',
                    help='Use variable width binning for the likelihood curves')

    pp.add_option('--tables', dest='tables', default=False, action='store_true',
                  help='Print acceptance and rejection values.')

    opts, args = pp.parse_args()
    if opts.input == '-' or opts.input == '--':
        opts.input = last_local_pickle()
        print 'Setting input to: ', opts.input
    if opts.input is not None and not os.path.isfile(opts.input):
        print 'Input specified but not a file!! ', opts.input
        sys.exit(1)

    for kk, vv in opts.__dict__.iteritems():
        if kk[-2:] == '_b' and vv is None and kk[:-2] != 'key':
            opts.__dict__[kk] = opts.__dict__[kk[:-2] + '_a']

    if len(args) > 0:
        print "There where unprocessed arguments: ", args

    return opts, args

def make_filename(ds_a, ds_b, rad, start, end, min_m, max_m, min_bdrs, max_bdrs, min_dr, max_dr):
    s_m = ('%.1fm%.1f' % (min_m , max_m )) if (min_m  > -1 or max_m  > -1) else ''
    s_b = ('%.1fbdrs%.1f' % (min_bdrs , max_bdrs )) if (min_bdrs  > -1 or max_bdrs  > -1) else ''
    s_d = ('%.2fd%.2f' % (min_dr, max_dr)) if (min_dr > -1 or max_dr > -1) else ''
    s_s = ('%.3f' % start) if 0 < abs(start) < 1 else '%d' % start
    s_e = ('%.3f' % end  ) if 0 < abs(end  ) < 1 else '%d' % end
    return '%s-vs-%s-%s_%s_%s_%s_%s' % (ds_a.u_name, ds_b.u_name,
                                        s_m, s_b, s_d, s_s, s_e)

def train_main(opts, ds_a, ds_b):
    print 'No input specified, will train'

    output = op.output
    if output is None:
        output = make_filename(ds_a, ds_b, op.rad, 
                                op.start, op.end, 
                                op.min_m, op.max_m,
                                op.min_bdrs, op.max_bdrs,
                                op.min_dr, op.max_dr)
        output += '.pickle'

#    for ds in [ds_a, ds_b]:
#        ds.get_array()
#        ds._ra = ds._ra[ds._ra['rel_subjets']['E'][:,1] >= 20]
    print 'Training on %s and %s from [%s, %s)' % (ds_a.u_name, ds_b.u_name, op.start, op.end)
    print 'Will save discriminants to ', output

    import discriminants as discs
    discs.train_all(ds_a, ds_b)

    discs.pickle_all(output)

def make_histograms(var_a, var_b, nBins=100, varbins=False, normed=False):

    if type(nBins) != type(1) and type(nBins) != type(1.):
        bins = np.array(nBins)
        h_a, bin_edges = np.histogram(var_a, bins, density=normed)
        h_b, bin_edges = np.histogram(var_b, bins, density=normed)
    elif varbins:
        bins = make_variable_bins(var_a, var_b, 1. / nBins)
        h_a, bin_edges = np.histogram(var_a, bins, density=normed)
        h_b, bin_edges = np.histogram(var_b, bins, density=normed)
    else:
        import scipy.stats as sps
        low  = min( sps.scoreatpercentile(var_a,  1), 
                    sps.scoreatpercentile(var_b,  1))
        high = max( sps.scoreatpercentile(var_a, 99), 
                    sps.scoreatpercentile(var_b, 99))
    
        h_a, bin_edges = np.histogram(var_a, nBins, (low, high))
        h_b, bin_edges = np.histogram(var_b, nBins, (low, high))
        
    return (h_a, h_b, bin_edges)


def plot_most_splendid_histogram(bin_centers, xlabel, h_a, name_a, h_b, name_b,
                                do_title):

    cc, li_diffs = samples.diff_prop_strs([name_a, name_b])
    lab_c = samples.make_verbose_title(cc)
    lab_a = samples.make_verbose_title(li_diffs[0])
    lab_b = samples.make_verbose_title(li_diffs[1])

    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)
    if True:
        #pretty but debatable
        avg_a = moving_average(h_a, 5)
        avg_b = moving_average(h_b, 5)
    else:
        #ugly but "fair"
        avg_a = moving_average(h_a, 1)
        avg_b = moving_average(h_b, 1)
    #do we really need error bars anyway??
    #ax.fill_between(bin_centers, avg_a - np.sqrt(avg_a), avg_a + np.sqrt(avg_a),
    #                color='g', alpha=0.5,)
    #ax.fill_between(bin_centers, avg_b - np.sqrt(avg_b), avg_b + np.sqrt(avg_b),
    #                color='b', alpha=0.5,)
    ax.plot(bin_centers, h_a, label=lab_a, color='g', linewidth=3)
    ax.plot(bin_centers, h_b, label=lab_b, color='b', linewidth=3)

    ax.set_xlabel(xlabel, size='x-large')
    ax.set_ylabel('Occupancy', size='x-large')
    ax.set_ylim(0,1.1*ax.get_ylim()[1])
    if do_title:
        ax.set_title(lab_c, size='x-large')

    ax.xaxis.set_major_formatter(_tick_fmt_float)
    ax.yaxis.set_major_formatter(_tick_fmt_int)
    ax.tick_params(axis='both', labelsize='large')
    ax.tick_params(which='both', axis='both', labelsize='large')

    ax.grid(lw=2.0, alpha=0.5)
    leg = plt.legend(
                    frameon=True, prop={'size':'15'},
                    loc=1,
                    )

    if leg is not None:
#this happens when label parsing above fails
        leg.get_frame().set_edgecolor('white')
        leg.get_frame().set_facecolor('white')

def ll_train(h_s, h_b):
    h_s = np.array(h_s.flatten(), dtype='f8')/np.sum(h_s)
    h_b = np.array(h_b.flatten(), dtype='f8')/np.sum(h_b)

    best_indices = np.argsort(h_s/(h_b + 1.0/1000))[::-1]

    return best_indices

def ll_test(h_s, h_b, best_indices):
    simulated_ll  = np.zeros((best_indices.shape[0],2))
    simulated_ll[:,0] = np.cumsum(h_s[best_indices])
    simulated_ll[:,1] = np.cumsum(h_b[best_indices])

    return simulated_ll

def eff_eff_pair(h_s, h_b):
    is_1d = len(h_s.shape) == 1

    h_s = np.array(h_s.flatten(), dtype='f8')/np.sum(h_s)
    h_b = np.array(h_b.flatten(), dtype='f8')/np.sum(h_b)

    best_indices = np.argsort(h_s/(h_b + 1.0/1000))[::-1]

    simulated_ll  = np.zeros((best_indices.shape[0],2))
    simulated_ll[:,0] = np.cumsum(h_s[best_indices])
    simulated_ll[:,1] = np.cumsum(h_b[best_indices])

    if not is_1d:
        return simulated_ll

    fixed_at_low  = np.zeros((h_s.shape[0],2))
    fixed_at_high = np.zeros((h_s.shape[0],2))
    fixed_at_low[:,0] = np.cumsum(h_s)
    fixed_at_low[:,1] = np.cumsum(h_b)
    fixed_at_high[:,0] = np.cumsum(h_s[::-1])
    fixed_at_high[:,1] = np.cumsum(h_b[::-1])

    return simulated_ll, fixed_at_low, fixed_at_high

def make_most_splendid_roc_curves(train_hists, test_hists = None):
    colors = ['#ff0000',  '#009999', '#00cc00',
              '#a60000', '#1d7373','#ff7400', '#008500',]
    fig = plt.figure(figsize=(7,5))
    ax = fig.add_subplot(111)

    for ii, pair in enumerate(sorted(train_hists.keys())):
        which_disc, is_linear = pair
        h_a, h_b, bin_edges = train_hists[pair]
        
        h_a = np.array(h_a, dtype='f8')/np.sum(h_a)
        h_b = np.array(h_b, dtype='f8')/np.sum(h_b)

        best_indices = ll_train(h_a, h_b)

        if test_hists:
            h_a, h_b, bin_edges = test_hists[pair]
            h_a = np.array(h_a, dtype='f8')/np.sum(h_a)
            h_b = np.array(h_b, dtype='f8')/np.sum(h_b)

        simulated_ll = ll_test(h_a, h_b, best_indices)

        simulated_ll = simulated_ll[simulated_ll[:,1] != 0]
        simulated_ll = simulated_ll[simulated_ll[:,0] != 1]

        marker = 's' if is_linear else 'o'

        cc = colors[ii%len(colors)]
        ax.semilogy(100*simulated_ll [:,0], 1.0/(.0001+simulated_ll [:,1]), 
                    marker=marker, label=which_disc,
                    mec=cc,
                    mfc=cc,
                    ls='None',
                    alpha=0.99)

        #Handy way to get some working points out
        #print '-'*50,'\n',which_disc,'sig_eff, then bkg_rej'
        #print 100*simulated_ll [:,0]
        #print 1.0/(.0001+simulated_ll [:,1])

    ax.set_ylabel('Background Rejection', size='x-large')
    ax.set_xlabel('Signal Efficiency [$\\%$]', size='x-large')

#    guess_x = np.linspace(0,1,1000)
#    ax.semilogy(100*guess_x, 1.0/(.000001 + guess_x), 'k--', 
#                label='Worst Case Guessing',
#                lw = 3
#                )

    #seq = ['$\\mathbf{10^{%d}}$' % ii for ii in range(0, 5)]
    ax.xaxis.set_major_formatter(_tick_fmt_int)
    #ax.yaxis.set_major_formatter(ticker.FixedFormatter(seq))
    ax.yaxis.set_major_formatter(_tick_fmt_int)
    _rng = [1,3,6]
    _rng += [10*xx for xx in _rng] + [100*xx for xx in _rng]
    ax.set_yticks(_rng)
    ax.set_ylim(.9, 250)
    #ax.set_yticks([10**i for i in range(0,3)])
#    ax.set_ylim(.9, 1.1e4)
#    ax.set_yticks([10**i for i in range(0,5)])
    ax.set_xlim(-5, 105)
    ax.set_xticks(np.arange(0,110,10))

    ax.tick_params(axis='both', labelsize='large')
    ax.tick_params(which='both', axis='both', labelsize='large')

    ax.grid(which='major', axis='both', lw=1.0, ls='--', color='gray')
    #ax.grid(which='minor', axis='both', lw=1.0, ls='--', color='gray')
    leg = plt.legend(
                    frameon=True, prop={'size':'15'},
                    loc=1, #upper right
                    scatterpoints=1, #show only one marker
                    )

    leg.get_frame().set_edgecolor('white')
    leg.get_frame().set_facecolor('white')
    
def make_most_splendid_plots(pp, order, ra_a, name_a, ra_b, name_b, varbins=False):
    import discriminants

    train_hists = {}
    test_hists = {}
    for dd in order:
        kk = dd.name()
        is_linear = isinstance(dd, discriminants.LinearDisc)
        h_a, h_b, bin_centers = make_histograms(ra_a[kk], ra_b[kk], 100, varbins)
        bin_centers = (bin_centers[:-1] + bin_centers[1:]) / 2.

        plot_most_splendid_histogram(bin_centers, dd.title(), 
                                    h_a, name_a, h_b, name_b,
                                    do_title=True)
        plt.savefig(pp, format='pdf')
        plot_most_splendid_histogram(bin_centers, dd.title(), 
                                    h_a, name_a, h_b, name_b,
                                    do_title=False)
        plt.savefig(pp, format='pdf')

        size_a, size_b = len(ra_a[kk])/2, len(ra_b[kk])/2
        test_hists[(kk, is_linear)] = make_histograms(ra_a[kk][:size_a], ra_b[kk][:size_b], 20, varbins)
        train_hists[(kk, is_linear)]  = make_histograms(ra_a[kk][size_a:], ra_b[kk][size_b:], test_hists[(kk, is_linear)][2], varbins)

    make_most_splendid_roc_curves(train_hists, test_hists)
    plt.savefig(pp, format='pdf')

    return train_hists, test_hists

def get_working_points(llh, li_wps=[2,10,20,100], fixed_bkg_rejection=True):
    import scipy.interpolate as inp
    import scipy.optimize    as opt

    if fixed_bkg_rejection:
        xx = llh[:,0]
        yy = 1.0/(.0001+llh[:,1])
    else:
        xx = 1.0/(.0001+llh[:,1])
        yy = llh[:,0]
        xx = xx[::-1]
        yy = yy[::-1]

    roc = inp.interp1d(xx, yy, kind='linear')
    low, high = np.min(xx), np.max(xx)

    lines = ['Sig Eff @ Bkg Rej']
    for y_wp in li_wps:
        try:
            x_wp = opt.brentq(lambda x: roc(x) - y_wp, low, high)
        except ValueError:
            x_wp = -1

        if fixed_bkg_rejection:
            di = {'bkg': str(y_wp),
                  'sig': ('% 3d' % int(100*x_wp)) if x_wp > -1 else ' <1'}
        else:
            di = {'bkg': ('%.1f' % x_wp) if x_wp > 01 else '??', 'sig':'% 3d' % int(100*y_wp)}

        

        lines.append('%s$\\%%$ @ x%s' % (di['sig'], di['bkg']))
    return '\n'.join(lines)

def make_variable_bins(var_a, var_b, bin_size = 0.05):
    quantiles = np.arange(0, 1.0 + bin_size, bin_size)
    a = np.sort(var_a)
    n = a.shape[0]
    quantiles = np.array(quantiles * n, dtype=int)
    if quantiles[-1] != n-1:
        quantiles[-1] = n-1
    if quantiles[0] != 0:
        quantiles[0] = 0
    output = a[quantiles]
    output = np.unique(output)
    output[0] = min(min(var_a), min(var_b))
    output[-1] = max(max(var_a), max(var_b))
    return np.unique(output)

def make_most_splendid_2dplots(pp, ra_a, ra_b, val_a, val_b, u_name_a, u_name_b,
                               title_a = None, title_b = None,
                                gridsize = 20, varbins = False, tables = False):

    print val_a, val_b, title_a, title_b
    from scipy.stats import pearsonr
    import samples
    test_cc , test_diffs  = samples.diff_prop_strs([u_name_a, u_name_b])
    polys = []

    fig = plt.figure()
    fig.suptitle('%s' % (samples.make_verbose_title(test_cc)))
    left  = min([np.min(rr[val_a]) for rr in [ra_a,ra_b]])
    right = max([np.max(rr[val_a]) for rr in [ra_a,ra_b]])
    down  = min([np.min(rr[val_b]) for rr in [ra_a,ra_b]])
    up    = max([np.max(rr[val_b]) for rr in [ra_a,ra_b]])

    size_a, size_b = len(ra_a[val_a])/2, len(ra_b[val_b])/2
    if varbins:
        xBreaks = make_variable_bins(ra_a[val_a], ra_b[val_a], 1. / gridsize)
        yBreaks = make_variable_bins(ra_a[val_b], ra_b[val_b], 1. / gridsize)

    split_counts = ([],[])
    for ii, pair in enumerate([(ra_a, test_diffs[0]), (ra_b, test_diffs[1])]):
        sample, label = pair
        ax = plt.subplot2grid((2,2), (0,ii))

        if not varbins:
            for s in range(2):
                split_counts[s].append(
                    ax.hexbin(sample[val_a][s*size_a:(s+1)*size_a], sample[val_b][s*size_b:(s+1)*size_b],
                              gridsize=gridsize, extent=(left, right, down, up), visible=False).get_array()
                    )
            ax.cla()
            poly = ax.hexbin(sample[val_a], sample[val_b], gridsize=gridsize, 
                             extent=(left, right, down, up))
            polys.append(poly)
        else:
            for s in range(2):
                split_counts[s].append(
                    np.histogram2d(sample[val_a][s*size_a:(s+1)*size_a], sample[val_b][s*size_b:(s+1)*size_b], bins=(xBreaks, yBreaks), normed=True)[0].reshape(1,-1)[0]
                    )

            H, xedges, yedges = np.histogram2d(sample[val_a], sample[val_b], bins=(xBreaks, yBreaks), normed=True)
            edges = (xedges, yedges)
            xgrid, ygrid = np.meshgrid(xedges, yedges)
            poly = ax.pcolormesh(xgrid, ygrid, H)
            ax.set_ylim(down, up)
            ax.set_xlim(left, right)
            polys.append(poly)

        ax.set_ylim(down, up)
        ax.set_xlim(left, right)
        if ii == 0:
            ax.set_ylabel(val_b)
        ax.set_xlabel(val_a)
        corr = pearsonr(sample[val_a], sample[val_b])
        ax.set_title("%s (corr = %.2f)" % (samples.make_verbose_title(label), corr[0]), size='medium')


    best_indices = ll_train(split_counts[0][0], split_counts[0][1])
    simulated_ll = ll_test(split_counts[1][0], split_counts[1][1], best_indices)

    #print simulated_ll
    #print val_a, val_b, u_name_a, u_name_b,title_a, title_b 
    simulated_ll[:,0]/= simulated_ll[-1,0]
    simulated_ll[:,1]/= simulated_ll[-1,1]
    if not varbins:
        simulated_ll = simulated_ll[simulated_ll[:,1] != 0]
        simulated_ll = simulated_ll[simulated_ll[:,0] != 1]

    ax = plt.subplot2grid((2,2), (1,0), colspan=2)
    ax.semilogy(100*simulated_ll [:,0], 1.0/(.0001+simulated_ll [:,1]), 
                'go', label='Simulated LHood', mew=0, mec='g')

    guess_x = np.linspace(0,1,1000)
    ax.semilogy(100*guess_x, 1.0/(.000001 + guess_x), 'k-', 
                label='Worst Case Guessing')

    ss_bkg_text = get_working_points(simulated_ll, [2, 10, 20, 100], True)
    ss_sig_text = get_working_points(simulated_ll, [.95, .90, .75, .50], False)

    if val_a == val_b and tables:
        print '\n\n %s' % (val_a)
        for a,b in (('Background Rejection:', ss_bkg_text), ('Signal Acceptance:', ss_sig_text)):
            print '  ', a
            print '\\hline\n' + b.replace('Eff @', 'Eff [%] &').replace('$\\%%$ @x', '& ').replace('\n', ' //\n').replace('$\\%%$', '').replace('Rej //', 'Rej //\n\hline') + ' // \n\\hline'

    
    #print ss_bkg_text

    ax.text(70, 100, ss_bkg_text)
    ax.text(40, 100, ss_sig_text)
    ax.set_ylabel('Background Rejection')
    ax.set_xlabel('Signal Efficiency ($\\%$)')
    ax.set_ylim(.9, 1e4)#1.1*ax.get_ylim()[1])
    ax.set_xlim(-5, 105)
    ax.set_xticks(np.arange(-10,110,10))
    ax.grid()

    plt.savefig(pp, format='pdf')

def cover_page(pp, opts, test_a, test_b):
    import time

    di_input = samples.parse_pickle_name(opts.input)
    train_cc, train_diffs = samples.diff_prop_strs([di_input['ds_a'], di_input['ds_b']])
    test_cc , test_diffs  = samples.diff_prop_strs([test_a, test_b])
    lines = [
    'Trained on: %s' % ' vs '.join(samples.make_verbose_title(ts) for ts in train_diffs),
    '$%(min_m)s < m_j < %(max_m)s, %(start)s < n_j < %(end)s$' % di_input,
    '$%(min_bdrs)s < bdrs_{m_j} < %(max_bdrs)s$' % di_input,
    '$%(min_dr)s < DR < %(max_dr)s$' % di_input,
    samples.make_verbose_title(train_cc),
    '','',
    'Tested on: %s' % ' vs '.join(samples.make_verbose_title(ts) for ts in test_diffs),
    '$%(min_m)s < m_j < %(max_m)s, %(start)s < n_j < %(end)s$' % eval(str(opts)),
    '$%(min_bdrs)s < bdrs_{m_j} < %(max_bdrs)s$' % eval(str(opts)),
    '$%(min_dr)s < DR < %(max_dr)s$' % eval(str(opts)),
    samples.make_verbose_title(test_cc),
    '','',time.strftime('%d/%m/%Y %H:%M %Z'), os.environ['USER']
            ]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    tt = '\n'.join((ll+' ') for ll in lines)
    ax.text(0.1, 0.4, tt)
    ax.set_axis_off()
    plt.savefig(pp, format='pdf')

def test_main(opts, ds_a, ds_b):
    print 'Input specified, will test'

    import cPickle
    from matplotlib.backends.backend_pdf import PdfPages

    output = op.output
    if output is None:
        output = make_filename(ds_a, ds_b, op.rad, op.start, op.end, 
                                op.min_m, op.max_m,
                                op.min_bdrs, op.max_bdrs,
                                op.min_dr, op.max_dr)
        output = '%s-test-%s.pdf' % (os.path.splitext(op.input)[0], output)

    name_a = ds_a.name
    name_b = ds_b.name
    di = cPickle.load(open(opts.input))

    print 'Considering ', [dd.name() for dd in di['discs']]
    discs = [dd for dd in di['discs'] if opts.re_discs is None 
                or re.search(opts.re_discs, dd.name()) is not None]
    print 'Testing on ', [dd.name() for dd in discs]

    dType = [(dd.name(), 'f8') for dd in discs]
    ra_a  = np.zeros(ds_a.get_array().shape, dType)
    ra_b  = np.zeros(ds_b.get_array().shape, dType)
    for dd in discs:
        t_s = time.time()
        ra_a[dd.name()] = dd.get_values(ds_a)
        ra_b[dd.name()] = dd.get_values(ds_b)
        print dd.name(), 'took', time.time() - t_s


    pp = PdfPages(output)
    cover_page(pp, opts, ds_a.u_name, ds_b.u_name)

    order = [dd for dd in discs]
    hists = make_most_splendid_plots(pp, order, ra_a, ds_a.u_name, ra_b, ds_b.u_name, varbins=opts.varbins)
    cPickle.dump(hists, open(output.replace('.pdf', '.pickle'), 'w'))
    for dd in discs:
        if not dd.plottable():
            continue

        dd.plot_from_scratch(do_title=True)
        plt.savefig(pp, format='pdf')
        dd.plot_from_scratch(do_title=False)
        plt.savefig(pp, format='pdf')

    musts = set(('SVC', 'Logistic', 'Fisher', 'NSubjet2Over1', 'Mass', 'BDRS'))
    musts = [dd.name() for dd in order if any(ss in dd.name() for ss in musts)]

    for name in musts:
        make_most_splendid_2dplots(pp, ra_a, ra_b, name, name,
                                   ds_a.u_name, ds_b.u_name,
                                   varbins=opts.varbins,
                                   tables=opts.tables
                                   )
    for ii in xrange(len(order)):
        for kk in xrange(ii+1,len(order)):
            make_most_splendid_2dplots(pp, ra_a, ra_b, 
                                        order[ii].name(), order[kk].name(),
                                        ds_a.u_name, ds_b.u_name,
                                        order[ii].title(), order[kk].title(),
                                        varbins=opts.varbins)

    pp.close()
    #for kk in ds_a.get_array().dtype.names:
    #    if not re.search(r'wtag_[0-9]+', kk):
    #        continue
    #    sig_eff = np.count_nonzero(ds_a.get_array()[kk]) * 100.0/ds_a.get_array().size
    #    bkg_eff = np.count_nonzero(ds_b.get_array()[kk]) * 100.0/ds_b.get_array().size

    #    print kk, sig_eff, bkg_eff, 100./bkg_eff

if __name__ == '__main__':
    import sys
    print sys.argv

    op, args = handleCLI()

    gd  = samples.get_dataset
    args = [op.start, op.end, op.min_m, op.max_m, op.min_bdrs, op.max_bdrs, op.min_dr, op.max_dr, op.min_pt, op.max_pt, op.n_b, op.subjet_e]
    if op.key_a is None:
        ds_a = gd(op.gen_a, op.me_a, op.pt_a, op.mu_a, 
                  op.rad, op.trim_a,op.smear_a, op.pdg_a,
                  *args)
    else:
        ds_a = samples.get_exact_dataset(op.key_a, *args)

    if op.key_b is None:
        ds_b = gd(op.gen_b, op.me_b, op.pt_b, op.mu_b, 
                  op.rad, op.trim_b, op.smear_b, op.pdg_b,
                  *args)
    else:
        ds_b = samples.get_exact_dataset(op.key_b, *args)


    s_a = ds_a.get_array().size
    s_b = ds_b.get_array().size
    print 'Dataset %s has %d jets' % (ds_a.u_name, s_a)
    print 'Dataset %s has %d jets' % (ds_b.u_name, s_b)
    if s_a * s_b == 0:
        print 'Both datasets must have jets in them!'
        sys.exit(1)
    if op.input is None: train_main(op, ds_a, ds_b)
    else:                test_main (op, ds_a, ds_b)
