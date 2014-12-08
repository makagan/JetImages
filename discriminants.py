'''Library for calculating discriminants for jets images
'''

__author__ = 'Josh Cogan jcogan@cern.ch'

import math
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm     as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker

n_pca_comps = 2
#cm_jw = colors.LinearSegmentedColormap.from_list('jv', ('white', '#0c5aa6'))
#cm_jv = colors.LinearSegmentedColormap.from_list('jw', ('white', '#ff9700'))
#cm_bi = colors.LinearSegmentedColormap.from_list('bi', ('white', '#ff9700'))

class Disc(object):
    def __init__(self, title = None):
        self.trained_on = None
        self.can_be_drawn = True
        if title:
            self._title = title

    @staticmethod
    def cardinal_number(ii):
        if ii==1: return '1st'
        if ii==2: return '2nd'
        else    : return '%dth' % ii

    @staticmethod
    def get_cell_values(*args, **kwargs):
        ras = []
        for arg in args:
            try:
                ra = arg.get_array()
            except AttributeError:
                ra = arg

            for kk, li in kwargs.iteritems():
                low, high = li
                if low  is not None:
                    ra = ra[ra[kk] >= low ]
                if high is not None:
                    ra = ra[ra[kk] <  high]

            ras.append(ra)

        max_allowed = min([xx.size for xx in ras])
        for ii, ra in enumerate(ras):
            if isinstance(ra, np.ndarray):
                ra = ra[:max_allowed]
            try:
                ra = ra['cells']['E']
            except (ValueError, IndexError):
                pass

            ras[ii] = ra

        return ras

    @staticmethod
    def _unit_norm(arg, power=2):
        if power == 0:
            return np.array(arg)
        tt = np.array(arg).transpose()
        ss = np.sum(np.power(tt, power), axis=0)
        if power > 0:
            ss = np.power(ss, 1.0/power)

        tt /= ss
        tt = tt.transpose()
        return tt

    @staticmethod
    def unit_norm(arg):
#### IF THIS ISN"T 2 YOU HAVE  PROBLEM
        return Disc._unit_norm(arg, 2)
#### /IF THIS ISN"T 2 YOU HAVE  PROBLEM

    def name(self):
        try:
            return self._name
        except AttributeError:
            return self.__class__.__name__

    def title(self):
        try:
            return self._title
        except AttributeError:
            return self.name()


    @classmethod
    def plottable(cls):
        return False

    def train(self, dataset_a, dataset_b):
        try:
            self.trained_on = (dataset_a.u_name, dataset_b.u_name)
        except AttributeError:
            self.trained_on = ('Unnamed Array', 'Unnamed Array')


    def get_values(self, ra_input):
        return np.zeros(ra_input.shape, dtype=float)

class LinearDisc(Disc):

    def __init__(self, title = None):
        self.mean = 0
        if title:
            self._title = title

    def name(self):
        return self._name

    @classmethod
    def plottable(cls):
        return True

    def scale_func(self, ds_a, ds_b):
        vals_a = self.get_values(ds_a)
        vals_b = self.get_values(ds_b)
        scale = max(np.max(np.abs(vals_a)), np.max(np.abs(vals_b)))
        if scale == 0:
            return 1

        scale = 1.0/scale

        if np.average(vals_a) > np.average(vals_b):
            scale *= -1

        return scale

    def plot_from_scratch(self, do_title):
        fig = plt.figure(figsize=(6, 5))
        size = .80
        ax = fig.add_axes((0.08, 0.12, size, size))
        a = self.plot_on_axes(ax, do_title)

        ax_cbar = fig.add_axes((0.84,0.12,.05, size))
        cb_fmt = ticker.FormatStrFormatter('$\\mathbf{%s}$')
        cbar = fig.colorbar(a, ax_cbar)#, format=cb_fmt)
        ax_cbar.set_xlabel('$\\mathbf{Cell}$\n$\\mathbf{Coefficient}$', 
                            labelpad=9)
        ax_cbar.tick_params(axis='both', labelsize='large')
        ax_cbar.tick_params(which='both', axis='both', labelsize='large')

    def plot_on_axes(self, ax, do_title=False):
        import math

        xx = self.comp
        ll = int(math.sqrt(xx.shape[0]))
        elem = xx.reshape(ll,ll)
        low = -.1
        high = ll/10.0 + .1

        #elem = shiftAndScale(elem, 0, 1)
        vmin = np.min(elem)
        vmax = np.max(elem)
        if vmin >= 0:
            ret = ax.imshow(elem,
                           cmap=cm.jet,
                           norm=colors.LogNorm(),
                           interpolation='nearest',
                           origin='lower',
                           extent=[low, high, low, high],
                           )
        else:
            cm_bi = colors.LinearSegmentedColormap.from_list('bi', 
                    [(0,'red'), (abs(vmin)/(vmax-vmin), 'white'),(1,'blue')])
            ret = ax.imshow(elem,
                           cmap=cm_bi,
                           interpolation='nearest',
                           origin='lower',
                           extent=[low, high, low, high],
                           )
        ax.set_ylabel('Q$_1$', fontsize='x-large', labelpad=9)
        ax.set_xlabel('Q$_2$', fontsize='x-large', labelpad=12)

        _ticks = np.linspace(0, ll/10.0, 6)
        ax_fmt = ticker.FormatStrFormatter('$\\mathbf{%s}$')
        ax.tick_params(axis='both', labelsize='large')
        ax.tick_params(which='both', axis='both', labelsize='large')
        ax.xaxis.set_major_formatter(ax_fmt)
        ax.yaxis.set_major_formatter(ax_fmt)
        ax.set_xticks(_ticks)
        ax.set_yticks(_ticks)
        if do_title:
            ax.set_title(self.name(), size='xx-large')

        return ret

    def get_values(self, ds_input):

        cells = Disc.get_cell_values(ds_input)[0]
        cells = Disc.unit_norm(cells)

        cells = cells - self.mean
        return np.dot(cells, self.comp.T)

class Dummy(LinearDisc):
    _name = ""
    _title = ""
    comp = np.array([1]*1)
    
    def get_values(self, jet):
        return 0

    def train(self, *args):
        pass

        

class DRLinear(Disc):
    '''Bins a LinearDiscriminant in DR'''
    _name = ' DR Binned '

    def __init__(self, LinDisc_inst, bins=[None,None]):
        self._nbins = len(bins) - 1
        self._bins = bins
        self._LinDisc_inst = LinDisc_inst
        self._name = str(self._nbins) + DRLinear._name + LinDisc_inst._name
        self._discs = [copy.deepcopy(LinDisc_inst) for ii in range(self._nbins)]

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)

        for i_bin in range(self._nbins):
            rng = self._bins[i_bin:i_bin+2]
            ras = Disc.get_cell_values(dataset_a, dataset_b, subjet_dr=rng)
            if ras[0].size * ras[1].size == 0:
                print 'NOT ENOUGH STATS'
                print i_bin, rng
                self._discs[i_bin] = Dummy()
                
            self._discs[i_bin].train(ras[0], ras[1])

    def get_bin(self, val):
#this assumes sorted!
        for i_bin,low_edge in enumerate(self._bins[:-1]):
            high_edge = self._bins[i_bin+1]
            if val < low_edge:
                return -1
            if val < high_edge or high_edge is None:
                return i_bin

        return -1

    def _get_values(self, which_bin, array):
        theout = self._discs[which_bin].get_values(array)
        return theout

    def get_values(self, ds_input):
#        gb = np.vectorize(self.get_bin)
#        gv = np.vectorize(self._get_values)
#
#        which_bin = gb(ds_input.get_array()['subjet_dr'])
        vals = np.zeros(ds_input.get_array().shape, dtype=float)
        #god damnit.
        #vals = np.where(which_bin == np.nan, np.nan, 
        #                gv(which_bin, ds_input.get_array()))
        #for i_bin in range(self._nbins):
        #    rng = self._bins[i_bin:i_bin+2]
        #    ras = Disc.get_cell_values(dataset_a, dataset_b, subjet_dr=rng)
        #    self._discs[i_bin].train(ras[0], ras[1])
        ra = ds_input.get_array()
        for ii, jet in enumerate(ds_input.get_array()):
            which_bin = self.get_bin(jet['subjet_dr'])
            if which_bin < 0: 
                vals[ii] = -10
            else:
                #assumes each disc maps a jet to -1 to 1 (OR SMALLER)
                vals[ii] = (1+2*which_bin)+\
                            self._discs[which_bin].get_values(jet)

#        for ii, i_bin in enumerate(which_bin):
#            if i_bin >= 0:
#                vals[ii] = self._get_values(i_bin, ra[ii])
#            else:
#                vals[ii] = np.nan

        return vals

    @classmethod
    def plottable(cls): return True

    def make_subtitle(self, i_fig):
        var = '\Delta R'
        if self._bins[i_fig] is not None:
            var = '%.1f \\leq %s' % (self._bins[i_fig], var)
        if self._bins[i_fig+1] is not None:
            var = '%s < %.1f'  % (var, self._bins[i_fig+1])

        return '$\\mathbf{%s}$' % var

    def plot_from_scratch(self, do_title):
        ncols = int(math.sqrt(self._nbins)) + 1
        nrows = int(math.ceil(self._nbins/float(ncols)))
        fig = plt.figure()
        if do_title:
            fig.suptitle(self._discs[0]._title, fontsize='xx-large')

        for i_fig in range(self._nbins):
            #print self._nbins, nrows, ncols, i_fig
            ax  = fig.add_subplot(nrows, ncols, i_fig+1, aspect=1)
            self._discs[i_fig].plot_on_axes(ax, False)
            ax.set_title(self.make_subtitle(i_fig))
            on_left   = (i_fig % ncols) == 0
            on_bottom = (i_fig / ncols) == (nrows - 1)

            if not on_left:
                ax.set_ylabel('')
                ax.set_yticklabels([])
            if not on_bottom:
                ax.set_xlabel('')


class VarDisc(Disc):
    def __init__(self, var, name = None, title = None):
        self._name = name if name else 'Var_%s' % var
        self._title = title if title else self._name
        self._var = var
    def get_values(self, ds_input):
        try:
            return ds_input.get_array()[self._var]
        except ValueError:
            return None

class MassDisc(Disc):
    def get_values(self, ds_input):
        return ds_input.get_array()['m']

class BDRSDisc(Disc):
    def get_values(self, ds_input):
        ra = ds_input.get_array()
        try:
            return ra['bdrs_mass']
        except ValueError:
            return  None

class SubjetDR(Disc):
    def get_values(self, ds_input):
        return ds_input.get_array()['subjet_dr']

class NSubjet2Over1(Disc):
    _title = r'N-Subjetiness ($\mathbf{\tau_2 / \tau_1}$)'
    def get_values(self, ds_input):
        return ds_input.get_array()['tau_2']/ds_input.get_array()['tau_1']

class NSubjet3Over2(Disc):
    _title = r'N-Subjetiness ($\mathbf{\tau_3 / \tau_2}$)'
    def get_values(self, ds_input):
        return ds_input.get_array()['tau_3']/ds_input.get_array()['tau_2']

class YBalance(Disc):
    def get_values(self, ds_input):
        return ds_input.get_array()['y_balance']

class LogBalance(Disc):
    def get_values(self, ds_input):
        return ds_input.get_array()['log_balance']

class WTagger(Disc):
    def get_values(self, ds_input):
        ra = ds_input.get_array()
        try:
            return ra['wtag_opt']
        except ValueError:
            return  None

def shiftAndScale(in_array, low=0, high=1):
    old_min = np.min(in_array)
    old_max = np.max(in_array)

    out_array = (in_array-old_min)*(high-low)/(old_max - old_min)
    out_array += low

    return out_array

class AllPCA(LinearDisc):
    _name = 'AllPCA'

    def __init__(self):
        self.pca = None

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)

        from sklearn.decomposition import RandomizedPCA
        ra_a, ra_b = Disc.get_cell_values(dataset_a, dataset_b)
#        ra_a = dataset_a.get_array()
#        ra_b = dataset_b.get_array()
#
#        ra_a = ra_a[:min(ra_a.size, ra_b.size)]['cells']['E']
#        ra_b = ra_b[:min(ra_a.size, ra_b.size)]['cells']['E']

        ra_a = Disc.unit_norm(ra_a)
        ra_b = Disc.unit_norm(ra_b)
        avg_a = np.average(ra_a, axis=0)
        avg_b = np.average(ra_b, axis=0)

        avg_a /= np.sqrt(np.sum(avg_a*avg_a))
        avg_b /= np.sqrt(np.sum(avg_b*avg_b))

#        self.avg_a = avg_a
#        self.avg_b = avg_b
        self.avg_1 = avg_a
        self.avg_2 = avg_b

        ra_all       = np.concatenate((ra_a, ra_b), axis=0)
        cond_all     = Disc.unit_norm(ra_all)
        self.avg_all = np.average(cond_all, axis=0)
        self.mean    = self.avg_all

        #not needed?? appears the pca method does this itself
        #cond_all -= self.avg_all

        import time
        self.pca = RandomizedPCA(n_pca_comps)
        ts = time.time()
        self.pca.fit(cond_all)
        print 'Fitting %s took %d seconds' % (self._name, time.time() - ts)

    def get_copy_slice(self, i_comp):
        if self.pca is None:
            return None

        return np.array(self.pca.components_[i_comp])

    def get_values(self, ds_input):
        cells = Disc.get_cell_values(ds_input)[0]
        cells = Disc.unit_norm(cells)
        cells -= self.avg_all

        return self.pca.transform(cells)


class PCAComp(LinearDisc):
    def __init__(self, i_comp, AllPCA_inst):
        LinearDisc.__init__(self)

        self._name  = '%s PCA Component' % Disc.cardinal_number(i_comp)
        self.i_comp = i_comp
        self.mean  = AllPCA_inst.mean
        self.comp  = AllPCA_inst.get_copy_slice(i_comp)
        self.trained_on = tuple(AllPCA_inst.trained_on)
        self.comp = Disc.unit_norm(self.comp)

class AvgComp(LinearDisc):
    def __init__(self, i_class, AllPCA_inst):
        LinearDisc.__init__(self)

        self._name  = '%s Class Average' % Disc.cardinal_number(i_class)
        self.i_class = i_class
        if i_class > 0:
            self.comp  = np.array(getattr(AllPCA_inst, 'avg_%d' % i_class)) 
        else:
            self.comp = np.array(AllPCA_inst.avg_1) - np.array(AllPCA_inst.avg_2)
            self._name  = 'Class Average A-B'

        self.comp = Disc.unit_norm(self.comp)

        self.trained_on = tuple(AllPCA_inst.trained_on)

class Width(LinearDisc):
    _name = 'JetWidth'
    _title = 'Jet Width'

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)
        ex = dataset_a.get_array()[0]['cells']['E']
        ncells = int(math.sqrt(ex.size))

        rr = np.mgrid[0:ncells,0:ncells]
        rr = (rr - ncells/2)/10.0
#because each bin is .1 not 1
        self.comp = np.sqrt(np.sum(np.square(rr), axis=0)).ravel()
        self.mean = 0

    def get_values(self, ds_input):

        cells = Disc.get_cell_values(ds_input)[0]
#this power = 1 means this function can't be rolled into std lineardisc
        cells = Disc._unit_norm(cells, power=1)

        cells = cells - self.mean
        return np.dot(cells, self.comp.T)


class LinearSVC(LinearDisc):
    _name = 'LinearSVC'
    _title = 'Linear Support Vector'

    def __init__(self, penalty='l2', loss='l2', dual=False, tol=1, C=1.0):
        self.mean = 0
    #    self._name += '_%s_%s_%s_%.1e_%.1e' % (penalty, loss, dual,tol, C)
        self.penalty = penalty
        self.loss = loss
        self.dual = dual
        self.tol  = tol
        self.C    = C

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)

        ras = Disc.get_cell_values(dataset_a, dataset_b)

        ra_all = np.concatenate((ras[0], ras[1]), axis=0)
        cond_all = Disc.unit_norm(ra_all)

        labels = np.concatenate((np.ones (ras[0].shape[0]),
                                 np.zeros(ras[1].shape[0])), axis=0)

        import time
        ts = time.time()

        from sklearn.svm import LinearSVC
        self._disc = LinearSVC(penalty=self.penalty, loss=self.loss, tol=self.tol, C=self.C,
                                dual=self.dual, fit_intercept=False)
        self._disc.fit(cond_all, labels)
        self.mean = 0
        self.comp = self._disc.coef_
        self.comp = self.comp.reshape((self.comp.shape[0] * self.comp.shape[1],))

        print 'Fitting %s took %d seconds' % (self._name, time.time() - ts)
        self.comp = self.comp * self.scale_func(dataset_a, dataset_b)

#    def get_values(self, ds_input):
#        cells = Disc.get_cell_values(ds_input)[0]
#        cells = Disc.unit_norm(cells)
#
#        print self._disc.transform(cells).shape
#        out = self._disc.predict(cells)
#
#        print out.shape
#        return out

class LogisticRegression(LinearDisc):
    _name = 'Logistic'
    _title = 'Logistic Regression'

    #these really are pretty optimal!
    def __init__(self, penalty='l2', tol=1e-3, reg=1.0):
        LinearDisc.__init__(self)
        #self._name = '%s_%.1f_%.1f' % (self.name, reg, tol)
        self.penalty = penalty
        self.tol = tol
        self.reg = reg

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)

        ras = Disc.get_cell_values(dataset_a, dataset_b)
            
        ra_all = np.concatenate((ras[0], ras[1]), axis=0)
        cond_all = Disc.unit_norm(ra_all)

        labels = np.concatenate((np.ones (ras[0].shape[0]),
                                 np.zeros(ras[1].shape[0])), axis=0)
        import time
        ts = time.time()

        from sklearn.linear_model import LogisticRegression
        self.logistic = LogisticRegression(penalty=self.penalty, tol=self.tol, C=self.reg)
        self.logistic.fit(cond_all, labels)
        self.comp = self.logistic.coef_
        self.comp = self.comp.reshape((self.comp.shape[0] * self.comp.shape[1],))

        print 'Fitting %s took %d seconds' % (self._name, time.time() - ts)
        self.comp = self.comp * self.scale_func(dataset_a, dataset_b)
        self.mean = 0

#    def get_values(self, ds_input):
#        cells = Disc.get_cell_values(ds_input)[0]
#        cells = Disc.unit_norm(cells)
#
#        cells = cells - self.mean
#        out = self.logistic.predict_proba(cells)
#        out = out[:,0].reshape((out.shape[0],))
#        return out



class AllFisher(LinearDisc):
    _name = 'Fisher'
    _title = 'Fisher-Jet'

    def __init__(self, tol=1e-3, do_smooth_reg=False, cov_class=None, cov_power=1):
        LinearDisc.__init__(self)
        self.tol = tol
        self.do_smooth_reg = do_smooth_reg
        self.cov_class = cov_class
        self.cov_power = cov_power

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)

        ras = Disc.get_cell_values(dataset_a, dataset_b)
        ra_all = np.concatenate((ras[0], ras[1]), axis=0)
        cond_all = Disc.unit_norm(ra_all)

        labels = np.concatenate((np.ones (ras[0].shape[0]),
                                 np.zeros(ras[1].shape[0])), axis=0)

        import time
        ts = time.time()

        from fisher import Fisher
        self.fish = Fisher(n_components=1)
        self.fish.fit(cond_all, labels, tol=self.tol, 
                      do_smooth_reg=self.do_smooth_reg, 
                      cov_class=self.cov_class,
                      cov_power=self.cov_power, store_covariance=True)
        self.comp = self.fish.w_[0]
        self.comp = Disc.unit_norm(self.comp)

        print 'Fitting %s took %d seconds' % (self._name, time.time() - ts)
        self.comp = self.comp * self.scale_func(dataset_a, dataset_b)

class Ellipse(LinearDisc):
    _name = 'Ellipse'

    def __init__(self, i_class, AllFish_inst, reg=1e-4):
        Disc.__init__(self)
        self._name  = '%s Class SqrDist_%.1e' % (i_class, reg)
        self.i_class = i_class
        self.means = np.array(AllFish_inst.fish.means_[0])
        self.covs  = np.array(AllFish_inst.fish.covs_ [0])
        self.trained_on = tuple(AllFish_inst.trained_on)
        self.AllFish_inst = AllFish_inst

        #print self.covs.shape
        self.invs  = np.empty(shape=self.covs.shape) #= np.linalg.pinv(self.covs, rcond=reg)
        self.comps = np.empty(shape=self.means.shape) #= np.linalg.pinv(self.covs, rcond=reg)
        for ii in range(self.covs.shape[0]):
            self.invs[ii]  = np.linalg.pinv(self.covs[ii], rcond=reg)
            self.comps[ii] = Disc.unit_norm(np.dot(self.invs[ii], self.means[ii].T))

        if self.i_class == 'A':
            self.comp= self.comps[0]
        elif self.i_class == 'B':
            self.comp= self.comps[1]
        elif self.i_class == 'r':
            self.comp= (self.comps[0] - self.comps[1])

    def train(self, *args):
        pass

    def get_values(self, ds_input):
        cells_per_event = Disc.get_cell_values(ds_input)[0]
        cells_per_event = Disc.unit_norm(cells_per_event)

        
        #dot goes last, second to last
        class_sqrdists = np.zeros(shape=(2, cells_per_event.shape[0]), dtype='f8')
        for i_class in range(2):
            #out = np.zeros(cells_per_event.shape[0],dtype='f8')
            for ii, cells in enumerate(cells_per_event-self.means[i_class]):
                tt = np.dot(self.invs[i_class], cells.T)
                yy = np.dot(cells, tt)
                class_sqrdists[i_class][ii] = yy
    
        if self.i_class == 'A':
            #out = class_sqrdists[0] + class_sqrdists[1]
            out = np.log10(class_sqrdists[0])
        elif self.i_class == 'B':
            out = np.log10(class_sqrdists[1])
        elif self.i_class == 'r':
            #out = (class_sqrdists[0] - class_sqrdists[1])/(class_sqrdists[0] + class_sqrdists[1])
            out = np.sqrt((class_sqrdists[1])/(class_sqrdists[0] + class_sqrdists[1]))
        return out


class KernelFisher(Disc):
    _name = 'KernelFisher'

    def __init__(self):
        self.mean = 0
        self.kfish = None

    def name(self):
        return self._name

    @classmethod
    def plottable(cls):
        return False

    @staticmethod
    def cardinal_number(ii):
        if ii==1: return '1st'
        if ii==2: return '2nd'
        else    : return '%dth' % ii

    def train(self, dataset_a, dataset_b, kernel="poly", gamma=1.0, degree=2.0, coef0=-1.0, sigma_sqrd = 1e-4, tol=1e-5, use_total_scatter=True, print_timing=True):
        Disc.train(self, dataset_a, dataset_b)

        ras = Disc.get_cell_values(dataset_a, dataset_b)
        #ras = {'a':dataset_a.get_array(), 'b':dataset_b.get_array()}
        #for kk, ra in ras.iteritems():
        #    ra = ra[:min([xx.size for xx in ras.values()])]
        #    ra = ra['cells']['E']
        #    #ra = ra.reshape(ra.shape[0], ra.shape[2])
        #    ras[kk] = ra
            
        ra_all = np.concatenate((ras[0], ras[1]), axis=0)
        cond_all = Disc.unit_norm(ra_all)

        labels = np.concatenate((np.ones (ras[0].shape[0]),
                                 np.zeros(ras[1].shape[0])), axis=0)

        import time
        ts = time.time()

        from fisher import KernelFisher
        if gamma is None:
            gamma = 1 / (1.0*cond_all.shape[1])
        self.kfish = KernelFisher(kernel=kernel, gamma=gamma, degree=degree, coef0=coef0,
                                  use_total_scatter=use_total_scatter, sigma_sqrd =sigma_sqrd, tol=tol, print_timing=print_timing)
        print self.kfish
        self.kfish.fit(cond_all, labels)
        print 'Fitting %s took %d seconds' % (self._name, time.time() - ts)


    def get_values(self, ds_input):
        if self.kfish is None:
            print "ERROR in KernelFisher:  tried to get values before training!"
            return 0
        
        ra = ds_input.get_array()
        cells = ra['cells']['E']
        #cells = cells.reshape(cells.shape[0], cells.shape[2])
        cells = Disc.unit_norm(cells)

        #cells = cells - self.mean
        return self.kfish.transform(cells)[:,0]

class LinearTree(Disc):
    '''Has a LinearDisc'''

    def __init__(self, nlayers=2):
        Disc.__init__(self)
        self._nlayers = nlayers
        self._ndiscs = (2**nlayers) - 1
        self._discs    = [[]]*(self._ndiscs)

    @staticmethod
    def point_of_max_entropy(vals_a, vals_b):
        from math import log
        low  = min(np.min(vals_a), np.min(vals_b))
        high = max(np.max(vals_a), np.max(vals_b))

        nbins = 1000
        h_a, edges = np.histogram(vals_a, bins=nbins, range=(low,high))
        h_b, edges = np.histogram(vals_b, bins=nbins, range=(low,high))

        epsilon = 1e-7
        pur_a = np.cumsum(h_a, dtype=float)/vals_a.size + epsilon
        pur_b = np.cumsum(h_b, dtype=float)/vals_b.size + epsilon
        entropy  = -1 * np.where(pur_a <= epsilon, 0, pur_a * np.log(pur_a))
        entropy += -1 * np.where(pur_b <= epsilon, 0, pur_b * np.log(pur_b))
        best_bin = np.argmax(entropy)
        return edges[best_bin + 1]

    def _get_next_indices(self, index):
        '''storing binary tree in array trick'''
        if index*2 + 2 >= self._ndiscs:
            return None, None

        return index*2+1, index*2+2

    def _rec_train(self, index, ra_a, ra_b):
        ld = Fisher()
        ld.train(ra_a, ra_b)
        vals_a = ld.get_values(ra_a)
        vals_b = ld.get_values(ra_b)

        split_here = LinearTree.point_of_max_entropy(vals_a, vals_b)
        go_left_a  = ra_a[vals_a <  split_here]
        go_right_a = ra_a[vals_a >= split_here]

        go_left_b  = ra_b[vals_b <  split_here]
        go_right_b = ra_b[vals_b >= split_here]

        #left has higher sig purity
        if go_left_a.size / go_left_b.size < go_right_a.size/go_right_b.size:
            go_left_a, go_right_a = go_right_a, go_left_a
            go_left_b, go_right_b = go_right_b, go_left_b

        left_index, right_index = self._get_next_indices(index)
        if left_index is None: return

        self._rec_train(left_index, go_left_a, go_left_b)
        self._rec_train(right_index, go_right_a, go_right_b)

    def train(self, dataset_a, dataset_b):
        Disc.train(self, dataset_a, dataset_b)
        ra_a = dataset_a.get_array()
        ra_b = dataset_b.get_array()

        self._rec_train(0, ra_a, ra_b)



_all_discs_ =  [
                VarDisc('m', 'MassDisc', 'Jet Invariant Mass'),
                VarDisc('bdrs_mass', 'BDRSDisc', 'Jet BDRS Mass'),
                VarDisc('subjet_dr', 'SubjetDr', r'Subjet $\Delta$R'),
                VarDisc('pt', 'pT', r'Jet $p_T$'),
                VarDisc('e', 'e', r'Jet $E$'),
                Width(),
                NSubjet2Over1(),
                #LogisticRegression(reg=.1),
                LinearSVC(),
#                DRLinear(AllFisher(), [None, 0.4, 0.6, 0.8, 1.0, 1.2,None]),
#                DRLinear(LogisticRegression(reg=0.0001), [None, 0.4, 0.6, 0.8, 1.0, 1.2,None]),
#                DRLinear(LogisticRegression(), [None, 0.4, 0.6, 0.7, 0.8, 0.9, 1.2,None]),
                VarDisc('y_balance', 'YBalance', 
                        r'Rapidity Balance ($p_T^{1} \times \Delta R_{1,2} / m$)'),
                #VarDisc('log_balance', 'LogBalance', 
                #        r'Log Balance ($- \log \left( 1-p_T^{1}/p_T \right)$)'),
                #VarDisc('wtag_opt', 'Cui,Han,Schwartz', 
                #        r'Cui, Han and Schwartz BDT'),
                #KernelFisher(),
                ]

def train_all(ds_a, ds_b):
    [dd.train(ds_a, ds_b) for dd in _all_discs_]
    if n_pca_comps > 0:
        allpca = AllPCA()
        allpca.train(ds_a, ds_b)
        _all_discs_.extend([PCAComp(ii, allpca) for ii in range(n_pca_comps)])

    _all_discs_.append(AvgComp(1 , allpca))
    _all_discs_.append(AvgComp(2 , allpca))
    _all_discs_.append(AvgComp(-1, allpca))

    allfish = AllFisher()
    allfish.train(ds_a, ds_b)
    _all_discs_.append(allfish)
##    _all_discs_.append(Ellipse('A', allfish, 1e-4))
##    _all_discs_.append(Ellipse('B', allfish, 1e-4))
#    _all_discs_.append(Ellipse('r', allfish, 1e-4))

def pickle_all(fname):
    import cPickle
    import sys
    import time

    fh = open(fname, 'w')
    tt = {'discs': list(_all_discs_)}
    tt['args'] = ' '.join(sys.argv)
    tt['time'] = time.strftime('%b-%d-%Y-%H-%M-%S')


    print 'Pickling %s' % ', '.join(dd.name() for dd in _all_discs_)
    cPickle.dump(tt, fh, 2)
    fh.close()

if __name__ == '__main__':
    #just some tests cases
    pass

    import samples

    aa = samples.get_dataset('pythia', 'w', 200, 0, 1.2, 0, start=0, end=100)
    bb = samples.get_dataset('pythia', 'v', 200, 0, 1.2, 0, start=0, end=100)
    cc = samples.get_dataset('pythia', 'v', 200, 0, 1.2, 0, start=1000, end=1010)

    lt = LinearSVC()
    lt.train(aa, bb)
