#!/usr/bin/env python
import time
import sys
import getopt
import math
import numpy                as np
import scipy.interpolate    as inp
import matplotlib.pyplot    as plt
import matplotlib.cm        as cm
np.set_printoptions(precision=3, suppress=True)

def makeTitle(index, nComps, doCombinedFittting):
    src = 'W' if index < nComps else 'G'
    if doCombinedFitting: 
        if index == nComps:
            return  'Dot Prod W-Avg'
        if index == nComps+1:
            return  'Dot Prod G-Avg'

        comp = index % nComps
        src='Combined'
        return '%dth %s Component' % (comp, src)


def shiftAndScale(in_array, low=0, high=1):
    old_min = np.min(in_array)
    old_max = np.max(in_array)

    out_array = (in_array-old_min)*(high-low)/(old_max - old_min)
    out_array += low

    return out_array

#def condition(*largs):
#
#    avg = None
#    outputs = []
#    for arg in largs:
#        tt = np.array(arg).transpose()
#        ss = np.sqrt(np.sum(tt*tt, axis=0))
#
#        tt /= ss
#        tt = tt.transpose()
#
#        if avg is None:
#            avg = np.average(tt, axis=0)
#
#        tt -= avg
#        outputs.append(tt)
#
#    outputs.append(avg)
#    return outputs

def unit_norm(arg):
    tt = np.array(arg).transpose()
    ss = np.sqrt(np.sum(tt*tt, axis=0))

    tt /= ss
    tt = tt.transpose()
    return tt

def average(arg):
    return np.average(tt.transpose(), axis=0)

def plot_components(comps):
    #in case last one is zero
    trunc_comps = comps[:-1]
    fig = plt.figure()
    ax_top = fig.add_subplot(211)
    ax_top.semilogy(trunc_comps)
    ax_top.set_ylabel("Explained Variance")
    ax_top.grid(which='major', linestyle='-', color='0.75')
    ax_top.grid(which='minor', linestyle=':', color='0.50')

    cum_comps = np.cumsum(trunc_comps)
    ax_bot = fig.add_subplot(212)
    ax_bot.semilogy(1-cum_comps)
    ax_bot.set_ylabel("1-Cumulative Expl. Var.")
    ax_bot.set_xlabel("PCA Component")
    ylim = ax_bot.get_ylim()
    ax_bot.set_ylim(ylim[0]*.8, max(ylim[1], 1.2))
    ax_bot.grid(which='major', linestyle='-', color='0.75')
    ax_bot.grid(which='minor', linestyle=':', color='0.50')



def plot_histogram(trans_w, trans_g, index, nComps=0, doCombinedFitting=True):
    import scipy.stats as sps

    nBins = 100
    if isinstance(index, str):
        ww = trans_w[index]
        gg = trans_g[index]
    else:
        ww = trans_w[:, index]
        gg = trans_g[:, index]
#    low  = min(np.min(gg), np.min(ww))
#    high = max(np.max(gg), np.max(ww))
    low  = min(sps.scoreatpercentile(gg,  1), sps.scoreatpercentile(ww,  1))
    high = max(sps.scoreatpercentile(gg, 99), sps.scoreatpercentile(ww, 99))

    h_w, bin_edges = np.histogram(ww, nBins, (low, high))
    h_g, bin_edges = np.histogram(gg, nBins, (low, high))
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ebkw = {'linewidth':1,}
    ax.errorbar(bin_centers, h_w, np.sqrt(h_w),label=w_label   ,color='g', **ebkw)
    ax.errorbar(bin_centers, h_g, np.sqrt(h_g),label=g_label   ,color='b', **ebkw)
    if type(index) == str:
        ax.set_xlabel(index, size='x-large')
    else:
        ax.set_xlabel(makeTitle(index, nComps, doCombinedFitting), size='x-large')
    ax.set_ylabel('Occupancy', size='x-large')
    ax.set_ylim(0,ax.get_ylim()[1])
    ax.grid()
    plt.legend()


def plot_scatter(trans_2d_w, trans_2d_g, xax_index, yax_index,
                 nComps, doCombinedFitting):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(trans_2d_w[:,xax_index], trans_2d_w[:,yax_index], c='g', lw=0)
    ax.scatter(trans_2d_g[:,xax_index], trans_2d_g[:,yax_index], c='b', lw=0)

    ax.set_xlabel(makeTitle(xax_index, nComps, doCombinedFitting))
    ax.set_ylabel(makeTitle(yax_index, nComps, doCombinedFitting))

def plot_roc(trans_w, trans_g, index):
    import scipy.stats as sps
    if isinstance(index, str):
        ww = trans_w[index]
        gg = trans_g[index]
    else:
        ww = trans_w[:, index]
        gg = trans_g[:, index]

#    low  = min(np.min(gg), np.min(ww))
#    high = max(np.max(gg), np.max(ww))
    low  = min(sps.scoreatpercentile(gg,  1), sps.scoreatpercentile(ww,  1))
    high = max(sps.scoreatpercentile(gg, 99), sps.scoreatpercentile(ww, 99))

    njets = ww.shape[0]
    nbins = max(10, min(30, njets/100))
    nbins = 20

    h_w, bin_edges = np.histogram(ww, bins=nbins, range=(low,high))
    h_g, bin_edges = np.histogram(gg, bins=nbins, range=(low,high))

    h_w = np.array(h_w, dtype='f8')/np.sum(h_w)
    h_g = np.array(h_g, dtype='f8')/np.sum(h_g)

    simulated_ll, fixed_at_low, fixed_at_high = eff_eff_pair(h_w, h_g)

    fig = plt.figure()
    ax = fig.add_subplot(111)


    #ax.plot(fixed_at_low [:,1], fixed_at_low [:,0], 'bo', label='Fixed At Low End')
    #ax.plot(fixed_at_high[:,1], fixed_at_high[:,0], 'ro', label='Fixed At High End')
    ax.semilogy(simulated_ll [:,0], 1.0/(.0001+simulated_ll [:,1]), 'go', label='Simulated LHood')

    ax.set_ylabel('Background Rejection', size='x-large')
    ax.set_xlabel('Signal Efficiency', size='x-large')
    #ax.set_ylim(0, 1.1*np.max(1.0/(.0001+simulated_ll[:,1])))
    ax.set_ylim(.9, 1.1*ax.get_ylim()[1])
    ax.set_xlim(-.1, 1.1)

    guess_x = np.linspace(0,1,1000)
    ax.semilogy(guess_x, 1.0/(.000001 + guess_x), 'k-', label='Worst Case Guessing')

    plt.suptitle('ROC for Index=%s' % str(index))
    ax.grid()
    plt.legend()

def get_vals(struct, index):
    if isinstance(index, str):
        return struct[index]
    return struct[:, index]


def plot_roc2(trans_w, trans_g, xax_index, yax_index, nComps=0, doCombinedFitting=True):
    wx = get_vals(trans_w, xax_index)
    gx = get_vals(trans_g, xax_index)
    wy = get_vals(trans_w, yax_index)
    gy = get_vals(trans_g, yax_index)

    low_x  = min(np.min(wx), np.min(gx))
    low_y  = min(np.min(wy), np.min(gy))
    high_x = max(np.max(wx), np.max(gx))
    high_y = max(np.max(wy), np.max(gy))
    nbins = max(5, min(30,int(math.sqrt(gy.shape[0]/100.0))))

    h_w, xedges, yedges = np.histogram2d(wx, wy, bins=(nbins, nbins), range=((low_x, high_x),(low_y, high_y)))
    h_g, xedges, yedges = np.histogram2d(gx, gy, bins=(nbins, nbins), range=((low_x, high_x),(low_y, high_y)))

    simulated_ll = eff_eff_pair(h_w, h_g)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(np.linspace(0,1,10), np.linspace(0,1,10), 'k-')

    ax.plot(simulated_ll [:,1], simulated_ll [:,0], 'go', label='Simulated LHood')

    ax.set_xlabel('Background Efficiency', size='x-large')
    ax.set_ylabel('Signal Efficiency', size='x-large')
    ax.set_ylim(-.1, 1.4)
    ax.set_xlim(-.1, 1.1)
    plt.suptitle('ROC for Index=%s & %s' % \
        (makeTitle(xax_index,nComps,doCombinedFitting),
         makeTitle(yax_index,nComps,doCombinedFitting)))
    ax.grid()
    plt.legend()

def divide_input(all_jets, n_max, divisions=(0.6, 0.2, 0.2)):

    slices = []
    divs = np.cumsum(n_max * np.array([0]+list(divisions)))
    for ii in range(divs.shape[0] - 1):
        slices.append(all_jets[divs[ii]:divs[ii+1]])

    return slices

#def spawn_avg_transformer(*largs):
#    avgs = []
#    for arg in args:

def eff_eff_pair(h_s, h_b):
    is_1d = len(h_s.shape) == 1

    h_s = np.array(h_s.flatten(), dtype='f8')/np.sum(h_s)
    h_b = np.array(h_b.flatten(), dtype='f8')/np.sum(h_b)

    best_indices = np.argsort(h_s/(h_b + 1.0/1000))[::-1]
    print best_indices

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

class LinearAndMVA:
    def __init__(self, n_pca_comps=5):
        self.n_pca_comps = n_pca_comps
        self.changers = []

    def fit(self, class_a, class_b):
        from sklearn.decomposition import RandomizedPCA
        from fisher import Fisher

        print 'Fitting on %d of class A and %d of class B' % (class_a.shape[0], class_b.shape[0])
        #cond_a, avg_a = condition(class_a)
        #cond_b, avg_b = condition(class_b)
        cond_a = unit_norm(class_a)
        cond_b = unit_norm(class_b)
        avg_a = np.average(cond_a, axis=0)
        avg_b = np.average(cond_b, axis=0)
        #print cond_a.shape
        #print avg_a.shape
        #sys.exit(0)

        avg_a /= np.sqrt(np.sum(avg_a*avg_a))
        avg_b /= np.sqrt(np.sum(avg_b*avg_b))

        self.avg_a = avg_a
        self.avg_b = avg_b

        class_all = np.concatenate((class_a, class_b), axis=0)
        cond_all = unit_norm(class_all)
        avg_all  = np.average(cond_all, axis=0)

        #not needed?? appears the pca method does this itself
        #cond_all -= avg_all

        labels = np.concatenate((np.ones (class_a.shape[0]),
                                 np.zeros(class_b.shape[0])), axis=0)

        self.pca = RandomizedPCA(self.n_pca_comps)
        ts = time.time()
        self.pca.fit(cond_all)
        _pca_comps = RandomizedPCA(100)
        _pca_comps.fit(cond_all)
        self.pca_explained_variance_ratio_ = _pca_comps.explained_variance_ratio_
        print 'Fitting PCA took:', time.time() - ts

        self.fish = Fisher(n_components=1)
        ts = time.time()
        self.fish.fit(cond_all, labels, tol=1e-3, do_smooth_reg=False)
        print 'Fitting Fish took:', time.time() - ts

        self.changers.extend( [
            (lambda xx: unit_norm(xx), '!'), #overwrite incoming trans values
            (lambda xx: xx-avg_all      , '!'),
            (lambda xx: np.sum(xx*avg_a,axis=1).reshape(xx.shape[0], 1) , 'avg_a'),
            (lambda xx: np.sum(xx*avg_b,axis=1).reshape(xx.shape[0], 1) , 'avg_b'),
            (lambda xx: self.pca .transform(xx), 'pca'),
            (lambda xx: self.fish.transform(xx), 'fish'),
            ])

        self.generate_out_format()

    def generate_out_format(self):
        dtype = []
        dtype = [loc for ff, loc in self.changers if loc not in ('!', 'pca')]
        dtype.extend(['pca_%d' % ii for ii in range(self.pca.components_.shape[0])])
        self.dtype = [(dd, float) for dd in dtype]


    def transform(self, incoming):
        tt = np.array(incoming, dtype='f8')
        out = np.zeros(shape=tt.shape[0], dtype=self.dtype)

        for ff, loc in self.changers:
            vals = ff(tt)

            if loc[0] == '!':
                tt = vals

            elif loc == 'pca':
                for ii in range(vals.shape[1]):
                    out['pca_%d' % ii] = vals[:,ii]

            else:
                out[loc] = vals[:,0]

        return out


def tree_trans(trainer, labels, test_w, test_g):
    from sklearn import tree 
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier

    #clf = DecisionTreeClassifier(max_depth=3)
    clf = RandomForestClassifier(n_estimators=700,max_depth=6, min_samples_split=10, min_samples_leaf=10)
    ts = time.time()
    clf.fit(trainer, labels)
    print 'Trees training: ', (time.time() - ts)

    ts = time.time()
    prob_w = np.zeros(test_w.shape[0], dtype=[('Forest Output Prob', 'f8')])
    prob_g = np.zeros(test_w.shape[0], dtype=[('Forest Output Prob', 'f8')])
    print prob_w.shape
    print test_w.shape
    print clf.predict_proba(test_w).shape
    print clf.predict_proba(test_w)[:,0].shape
    prob_w['Forest Output Prob'] = clf.predict_proba(test_w)[:,0]
    prob_g['Forest Output Prob'] = clf.predict_proba(test_g)[:,0]
    print 'Trees transforming: ', (time.time() - ts)
    return prob_w, prob_g

    h_w, bin_edges = np.histogram(prob_w, 20, (0,1))
    h_g, bin_edges = np.histogram(prob_g, 20, (0,1))
    bin_centers = (bin_edges[0:-1] + bin_edges[1:])/2

    fig = plt.figure()
    ebkw = {'linewidth':1,}
    ax = fig.add_subplot(111)
    ax.errorbar(bin_centers, h_w, np.sqrt(h_w),label=w_label  ,color='g', **ebkw)
    ax.errorbar(bin_centers, h_g, np.sqrt(h_g),label=g_label ,color='b', **ebkw)
    ax.set_xlabel('Decision Tree Ouput Prob', size='x-large')
    ax.set_ylabel('Occupancy', size='x-large')
    plt.legend()


def flatten_rec(rec_array):
    out = np.zeros((rec_array.shape[0], len(rec_array.dtype)))
    for ii,dd in enumerate(sorted(rec_array.dtype.names)):
        out[:,ii] = rec_array[dd]

    return out

def retrieve_arrays(globbable_path):
    import glob
    import os
    path = os.path.expanduser(globbable_path)
    file_list = glob.glob(path)

    output_ra = None
    ra_list = [np.load(open(ff)) for ff in file_list]
    return np.concatenate(ra_list)
#    for ff in file_list:
#        this_ra = np.load(open(ff))
#        if output_ra is None:
#            output_ra = this_ra
#        else:
#            output_ra = np.concatenate((output_ra, this_ra))
#
#    return output_ra

def prepare_jet_arrays(whole_jets, parton_test, mass_range, max_jets = -1):
    if whole_jets == None:
        return None

    whole_jets  = whole_jets [(parton_test(np.abs(whole_jets ['pdgIDHardParton']))) & (whole_jets['m'] > mass_range[0]) & (whole_jets['m'] < mass_range[1]) ]
    if max_jets > 0:
        whole_jets = whole_jets[:max_jets]
    whole_jets   = whole_jets ['cells'].reshape(whole_jets .shape[0], 625)
    ##below line should require at least a cell > 0 GeV
    whole_jets   = whole_jets [np.max(whole_jets , axis=1) > 0]

    return whole_jets

def main(argv):

    #usePUSample = True
    #usePUSample = False
    #output_filename = 'pca_mva_train_500_0_test_30.pdf'
    
    trainWithPU = False
    testWithPU = False
    ptbin = 200

    try:
        opts, args = getopt.getopt(argv,"hp:r:t:",["ptbin=","trainPU=","testPU="])
    except getopt.GetoptError:
        print 'pca.py -p <pt bin> -r <train on pileup> -t <test on pileup>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'pca.py -p <pt bin> -r <train on pileup> -t <test on pileup>'
            sys.exit()
        elif opt in ("-p", "--ptbin"):
            ptbin = arg
        elif opt in ("-r", "--trainPU"):
            trainWithPU = bool(arg)
        elif opt in ("-t", "--testPU"):
            testWithPU = bool(arg)

    base_dir = '/u/eb/joshgc/mynfs/CSJets/logs/'
    samples = {"pythia_g" : (base_dir + "bsub_v_", "Pythia Light Jets", lambda x : x != 24),
               "pythia_w" : (base_dir + "bsub_w_", "Pythia W Jets"    , lambda x: x == 24),
               "herwig_g" : (base_dir + "bsub_v_", "Herwig Light Jets", lambda x: x != 24),
               "herwig_w" : (base_dir + "bsub_w_", "Herwig W Jets"    , lambda x: x == 24),
               }

    mass_range = (60, 100)

    max_training_jets = 10000
    max_training_jets = 1000
    max_testing_jets  = 30000

    global g_basename, g_label, g_parton_test
    global w_basename, w_label, w_parton_test

    training_sample_a, training_sample_b = ("pythia_g", "pythia_w") if len(args) < 2 else args[0:2]
    testing_sample_a, testing_sample_b = (training_sample_a, training_sample_b) if len(args) < 4 else args[2:4]

    g_basename, g_label, g_parton_test = samples[training_sample_a]
    w_basename, w_label, w_parton_test = samples[training_sample_b]
    g_test_basename = samples[testing_sample_a][0]
    w_test_basename = samples[testing_sample_b][0]
     

    output_filename = 'pca_mva_pt_'+str(ptbin)+'_train_'+str(30 if trainWithPU else 0)+'_test_'+str(30 if testWithPU else 0) + '_' + training_sample_a + '_vs_' + training_sample_b + '_tested_' + testing_sample_a + '_vs_' + testing_sample_b + '.pdf'
    print "---------------------------------------------------------------------"
    print "Running pca.py, output file is: ", output_filename
    print "---------------------------------------------------------------------"


    print "Read sample A from %s" % g_basename
    print "Read sample B from %s" % w_basename
    print "Test sample A from %s" % g_test_basename
    print "Test sample B from %s" % w_test_basename
    

    whole_g_0_jets  = retrieve_arrays(g_basename + str(ptbin) + "_0_1.2_*.npy" )
    whole_w_0_jets  = retrieve_arrays(w_basename + str(ptbin) + "_0_1.2_*.npy" )
    test_g_0_jets  = None if g_basename == g_test_basename else retrieve_arrays(g_test_basename + str(ptbin) + "_0_1.2_*.npy" )
    test_w_0_jets  = None if w_basename == w_test_basename else retrieve_arrays(w_test_basename + str(ptbin) + "_0_1.2_*.npy" )
#     whole_g_30_jets = retrieve_arrays("../npy_data/merge_bsub_v_"+str(ptbin)+"_30_1.2*.npy")
#     whole_w_30_jets = retrieve_arrays("../npy_data/merge_bsub_w_"+str(ptbin)+"_30_1.2*.npy")

    print "Number of jets:"
    print "Training sample A: %d" % len(whole_g_0_jets)
    print "Training sample B: %d" % len(whole_w_0_jets)
    if test_g_0_jets != None:
        print "Testing sample A: %d" % len(test_g_0_jets)
    if test_g_0_jets != None:
        print "Testing sample B: %d" % len(test_w_0_jets)

    many_w_0_jets  = prepare_jet_arrays(whole_w_0_jets, w_parton_test, mass_range, max_training_jets + max_testing_jets)
    many_g_0_jets  = prepare_jet_arrays(whole_g_0_jets, g_parton_test, mass_range, max_training_jets + max_testing_jets)
    test_many_w_0_jets  = prepare_jet_arrays(test_w_0_jets, w_parton_test, mass_range, max_testing_jets)
    test_many_g_0_jets  = prepare_jet_arrays(test_g_0_jets, g_parton_test, mass_range, max_testing_jets)

    if not trainWithPU and not testWithPU:
        n_test_remain = 0
        n_remain = min(many_w_0_jets.shape[0], many_g_0_jets.shape[0])
        f_train = min(0.5, float(max_training_jets)/n_remain)
        f_test  = float(min(max_testing_jets, n_remain - max_training_jets))/n_remain
        train_w, test_w, validate_w = divide_input(many_w_0_jets,n_remain, divisions=(f_train, f_test, 0.0))
        train_g, test_g, validate_g = divide_input(many_g_0_jets,n_remain, divisions=(f_train, f_test, 0.0))
        n_train_remain = train_w.shape[0]
         
        if test_many_w_0_jets != None:
            print "Get other W testing sample"
            test_w = test_many_w_0_jets[:max_testing_jets]
        if test_many_g_0_jets != None:
            print "Get other g testing sample"
            test_g = test_many_g_0_jets [:max_testing_jets]
        n_train_remain = train_w.shape[0]
        
#         del whole_w_0_jets
#         del whole_g_0_jets
        del many_w_0_jets
        del many_g_0_jets
         
    elif trainWithPU and testWithPU:
        n_test_remain = 0
        n_remain = min(many_w_30_jets.shape[0], many_g_30_jets.shape[0])
        train_w, test_w, validate_w = divide_input(many_w_30_jets,n_remain)
        train_g, test_g, validate_g = divide_input(many_g_30_jets,n_remain)
        n_train_remain = train_w.shape[0]
        
    elif not trainWithPU and testWithPU:
        n_train_remain = min(many_w_0_jets.shape[0] , many_g_0_jets.shape[0])
        n_test_remain  = min(many_w_30_jets.shape[0], many_g_30_jets.shape[0])
        #n_train_remain = min(5000, n_train_remain)
        train_w = many_w_0_jets [:n_train_remain]
        test_w  = many_w_30_jets[:n_test_remain]
        train_g = many_g_0_jets [:n_train_remain]
        test_g  = many_g_30_jets[:n_test_remain]

    elif trainWithPU and not testWithPU:
        n_train_remain = min(many_w_30_jets.shape[0] , many_g_30_jets.shape[0])
        n_test_remain  = min(many_w_0_jets.shape[0], many_g_0_jets.shape[0])
        #n_train_remain = min(5000, n_train_remain)
        train_w = many_w_30_jets [:n_train_remain]
        test_w  = many_w_0_jets[:n_test_remain]
        train_g = many_g_30_jets [:n_train_remain]
        test_g  = many_g_0_jets[:n_test_remain]


    print 'Working with a total input of ', n_train_remain + n_test_remain
    print 'Training W size: ', train_w.shape
    print 'Training G size: ', train_g.shape
    print 'Testing W size: ' , test_w .shape
    print 'Testing G size: ' , test_g .shape
    n_pca_comps = 12
    doCombinedFitting = True
    lmva = LinearAndMVA(n_pca_comps)
    lmva.fit(train_w, train_g)
    print "Fisher singular values", lmva.fish.singular_vals

    ts = time.time()
    trans_train_w = lmva.transform(train_w)
    trans_train_g = lmva.transform(train_g)
    trans_test_w  = lmva.transform(test_w)
    trans_test_g  = lmva.transform(test_g)
    print 'Transforming LMVA took:', time.time() - ts


    print 'Done fitting.  Drawing.'

    trans_2d_w = trans_test_w[:100,]
    trans_2d_g = trans_test_g[:100,]


    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(output_filename)

    flat_trans_train_w = flatten_rec(trans_train_w)
    flat_trans_train_g = flatten_rec(trans_train_g)
    flat_trans_test_w  = flatten_rec(trans_test_w)
    flat_trans_test_g  = flatten_rec(trans_test_g)

    flat_trans_train_all = np.concatenate((flat_trans_train_w, flat_trans_train_g), axis=0)
    labels = np.concatenate((np.ones (flat_trans_train_w.shape[0]),
                             np.zeros(flat_trans_train_g.shape[0])), axis=0)

    prob_w, prob_g = tree_trans(flat_trans_train_all, labels, flat_trans_test_w, flat_trans_test_g)
    print prob_w.shape
    print prob_g.shape
    plot_histogram(prob_w, prob_g, 'Forest Output Prob')
    plt.savefig(pp, format='pdf')

    plot_roc(prob_w, prob_g, 'Forest Output Prob')
    plt.savefig(pp, format='pdf')

    plot_components(lmva.pca_explained_variance_ratio_)
    plt.savefig(pp, format='pdf')


    print 'Plotting 1d'
    for xax_index in trans_test_w.dtype.names:
        plot_histogram(trans_test_w, trans_test_g, xax_index, n_pca_comps, doCombinedFitting)
        plt.savefig(pp, format='pdf')

        plot_roc(trans_test_w, trans_test_g, xax_index)
        plt.savefig(pp, format='pdf')

    print 'Done Plotting 1d'
#    for xax_index in range(trans_train_w.shape[1]):
#        for yax_index in range(xax_index+1, trans_train_w.shape[1]):
#            plot_scatter(trans_train_2d_w, trans_train_2d_g, xax_index, yax_index,
#                         n_pca_comps, doCombinedFitting)
#            plt.savefig(pp, format='pdf')
#
#            plot_roc2(trans_train_w, trans_train_g, xax_index, yax_index,
#                         n_pca_comps, doCombinedFitting)
#            plt.savefig(pp, format='pdf')
#


    elem_w = shiftAndScale(lmva.avg_a.reshape(25,25), 0, 1)
    elem_g = shiftAndScale(lmva.avg_b.reshape(25,25), 0, 1)

    fig = plt.figure()
    fig.suptitle('Average Calo Image')
    ax_w= fig.add_subplot(121, aspect=1)
    ax_g= fig.add_subplot(122, aspect=1)
    ax_w.set_xlabel(w_label, size='x-large')
    ax_g.set_xlabel(g_label, size='x-large')
    ax_w.imshow(elem_w, interpolation='nearest')
    ax_g.imshow(elem_g, interpolation='nearest')
    #fig.colorbar()
    plt.savefig(pp, format='pdf')

    elem_all = lmva.fish.w_[0].reshape(25,25)

    #if np.sum(elem_all[0:5,0:5]) < 0:  elem_all *= -1
    #elem_all = shiftAndScale(elem_all)

    fig = plt.figure()
    fig.suptitle('FisherJet')
    ax_all= fig.add_subplot(111, aspect=1)
    a = ax_all.imshow(elem_all, interpolation='nearest', cmap=cm.Blues)
    fig.colorbar(a)
    plt.savefig(pp, format='pdf')

    for iComp in range(n_pca_comps):

        elem_all = lmva.pca.components_[iComp].reshape(25,25)
        var_all  = lmva.pca.explained_variance_ratio_[iComp]

        #if np.sum(elem_all[0:5,0:5]) < 0:  elem_all *= -1
        #elem_all = shiftAndScale(elem_all)

        fig = plt.figure()
        fig.suptitle('%dth PCA Component' % iComp)
        ax_all= fig.add_subplot(111, aspect=1)
        ax_all.set_title('Expl. Var. = %.3f' % (var_all))
        a = ax_all.imshow(elem_all, interpolation='nearest', cmap=cm.Blues)
        fig.colorbar(a)
        plt.savefig(pp, format='pdf')

    whole_w_jets = whole_w_30_jets if testWithPU else whole_w_0_jets
    whole_g_jets = whole_g_30_jets if testWithPU else whole_g_0_jets
    #whole_w_jets = whole_w_30_jets
    #whole_g_jets = whole_g_30_jets
    for kk in ['wtag_50', 'preTrimMass', 'subjet_dr', 'm', 'e', 'eta', 'pt', 'phi', 'hsf', 'npv']:
        try:
            whole_w_jets[kk]
        except ValueError:
            continue
        plot_histogram(whole_w_jets[:n_train_remain], whole_g_jets[:n_train_remain], kk)
        plt.savefig(pp, format='pdf')

        plot_roc(whole_w_jets[:n_train_remain], whole_g_jets[:n_train_remain], kk)
        plt.savefig(pp, format='pdf')

    pp.close()


if __name__ == '__main__':
    main(sys.argv[1:])
