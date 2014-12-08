'''Library for retriving the dataset arrays of jets
Users should only need to use get_dataset, developers might need to expand
class dataset'''

__author__ = 'Josh Cogan jcogan@cern.ch'

###Priming the pump
_di_pdg_rng_ = {}
_di_pdg_rng_['v'] = set([      -3,-2,-1,1,2,3    ,21])  #u d s and later g
_di_pdg_rng_['q'] = set([-5,-4,-3,-2,-1,1,2,3,4,5])     #u d s c b
_di_pdg_rng_['c'] = set([-4,4])
_di_pdg_rng_['b'] = set([-5,5])
_di_pdg_rng_['t'] = set([-6,6])
_di_pdg_rng_['g'] = set([21])
_di_pdg_rng_['w'] = set([-24,24])
_di_pdg_rng_['h'] = set([-25,25])
###/Priming the pump

###Pretty names
pretty = {'gen'   :{'p': 'Pythia8', 'h':'Herwig', 'm':'MadGraph'},
          'me'    :{'w': 'W', 'v': 'Light', 'q':'Quark', 'g':'Gluon', 'h':'Higgs'},
          'trim'  :{'0': 'Untrimmed', '1': 'Trimmed'},
          'smear' :{'0': 'No smearing', '1': 'Radial Smearing'},
          'mu'    :{'0': 'No PU', '30': 'w/ PU'},
          'ptbin' :{'200': 'Low pT', '500': 'High pT'},
          'rad'   :{'1.2': 'R=1.2', '.4': 'R=0.4'},
}


###/Pretty names

class CutEfficiency:
    def __init__(self, label = ""):
        self._count = []
        self._name = []
        self._idx = -1
        self._label = label

    def begin(self):
        self._idx = 0

    def passed(self, name, count = 1):
        if name != "" and self._idx == len(self._name):
            self._name.append(name)
            self._count.append(0)
        if name == "" or name == self._name[self._idx]:
            self._count[self._idx] += count
            self._idx += 1
        
    def summary(self):
        print ""
        print "Cut summary: %s" % self._label
        eff = [1.] + [float(self._count[i+1]) / float(self._count[i]) for i in range(len(self._count)-1)]
        print "\\hline"
        print "  Cut       & Counts (Eff [%]) //"
        print "\\hline"
        for n, c, e in zip(self._name, self._count, eff):
            print "  %s%s & %d (%.3f) //" % (n, " "*(10 - len(n)), c, e)
        print "\\hline"

class dataset:
    def __init__(self, path, u_name='', name='', s_pdgid=None, l_f_cuts=None):
        self.path     = path
        self.u_name   = u_name
        self.name     = name
        self.s_pdgid  = set() if s_pdgid  is None else s_pdgid
        self.l_f_cuts = []    if l_f_cuts is None else l_f_cuts
        self._ra = None
        self.start  = 0
        self.end    = -1
        self.min_m  = -1
        self.max_m  = -1
        self.min_bdrs  = -1
        self.max_bdrs  = -1
        self.min_dr = -1
        self.max_dr = -1
        self.min_pt = -1
        self.max_pt = -1

    def set_mass_range(self, min_m, max_m):
        self.min_m = min_m
        self.max_m = max_m

    def set_bdrs_range(self, min_bdrs, max_bdrs):
        self.min_bdrs = min_bdrs
        self.max_bdrs = max_bdrs

    def set_dr_range(self, min_dr, max_dr):
        self.min_dr = min_dr
        self.max_dr = max_dr

    def set_pt_range(self, min_pt, max_pt):
        self.min_pt = min_pt
        self.max_pt = max_pt

    def set_n_b(self, n_b):
        self.n_b = n_b

    def set_subjet_e(self, subjet_e):
        self.subjet_e = subjet_e

    def set_range(self, start, end):
        self.start = start
        self.end   = end 
        self._do_range_math()

    def _do_range_math(self):
        if self._ra is None:
            return

        ll = self._ra.shape[0]
        if abs(self.start) < 1:
            self.start = int(ll * self.start)
        if abs(self.end) < 1:
            self.end = int(ll * self.end)

    @staticmethod
    def apply_pdgid_cut(ra, pdgs):
        if len(pdgs) == 0:
            return ra
        import numpy as np
        ff = np.vectorize(lambda pdg: pdg in pdgs)
        return ra[ff(ra['pdgIDHardParton'])]

    @staticmethod
    def apply_mass_cut(ra, min_m, max_m):
        if min_m >= 0:
            ra = ra[ra['m'] >= min_m]

        if max_m >= 0:
            ra = ra[ra['m'] <= max_m]

        return ra

    @staticmethod
    def apply_bdrs_cut(ra, min_bdrs, max_bdrs):
        if min_bdrs >= 0:
            ra = ra[ra['bdrs_mass'] >= min_bdrs]

        if max_bdrs >= 0:
            ra = ra[ra['bdrs_mass'] <= max_bdrs]

        return ra

    @staticmethod
    def apply_dr_cut(ra, min_dr, max_dr):
        if min_dr >= 0:
            ra = ra[ra['subjet_dr'] >= min_dr]

        if max_dr >= 0:
            ra = ra[ra['subjet_dr'] <= max_dr]

        return ra


    @staticmethod
    def apply_pt_cut(ra, min_pt, max_pt):
        if min_pt >= 0:
            ra = ra[ra['pt'] >= min_pt]

        if max_pt >= 0:
            ra = ra[ra['pt'] <= max_pt]

        return ra


    def clear(self):
        self._ra = None

    def get_array(self):
        #don't recalculate if already found
        if self._ra is not None: return self._ra

        import glob
        import os
        path = os.path.expanduser(self.path)
        file_list = glob.glob(path)

        import numpy as np
        ra_list = []
        self._orig_size = 0
        size_so_far     = 0
        dtype = None
        cut_eff = CutEfficiency(self.name)
        if len(file_list) == 0:
            print 'No files were found matching your expression: %s' % path
            print 'Could be an AFS PROBLEM!!'
            print 'Meh, I should exit here but wheres the fun in that?'
        for ff in file_list:
            try:
                ra_temp = np.load(open(ff))
            except ValueError:
                print ff, "Bad file: ValueError, continue"
                continue
            except IOError:
                print ff, "Bad file: IOError, continue"
                continue
            if dtype == None:
                dtype = ra_temp.dtype
            if ra_temp.dtype != dtype:
                ra_temp = np.array(ra_temp, dtype=dtype)
                print ff, "Bad dtype, overwrite"

            self._orig_size += ra_temp.shape[0]

            cut_eff.begin()
            cut_eff.passed('Initial', ra_temp.shape[0])

            #pdgid cuts
            if len(self.s_pdgid) > 0:
                ra_temp = dataset.apply_pdgid_cut(ra_temp, self.s_pdgid)
                cut_eff.passed('PDG ID', ra_temp.shape[0])
            if self.n_b >= 0:
                ra_temp = ra_temp[ra_temp['n_b'] == self.n_b]
                cut_eff.passed('N b', ra_temp.shape[0])
            #at least one cell is > 0
            ra_temp = ra_temp[np.max(ra_temp['cells']['E'], axis=1) > 0]
            cut_eff.passed('Cell E', ra_temp.shape[0])
            if self.subjet_e >= 0:
                ra_temp = ra_temp[ra_temp['subjets']['E'][:,1] > self.subjet_e]
                cut_eff.passed('Subjet E', ra_temp.shape[0])
            if self.min_dr >= 0 or self.max_dr >= 0:
                ra_temp = dataset.apply_dr_cut  (ra_temp, self.min_dr, self.max_dr)
                cut_eff.passed('Delta R', ra_temp.shape[0])
            if self.min_pt >= 0 or self.max_pt >= 0:
                ra_temp = dataset.apply_pt_cut(ra_temp, self.min_pt, self.max_pt)
                cut_eff.passed('pT', ra_temp.shape[0])
            if self.min_m >= 0 or self.max_m >= 0:
                ra_temp = dataset.apply_mass_cut(ra_temp, self.min_m , self.max_m )
                cut_eff.passed('Mass', ra_temp.shape[0])
            if self.min_bdrs >= 0 or self.max_bdrs >= 0:
                ra_temp = dataset.apply_bdrs_cut(ra_temp, self.min_bdrs , self.max_bdrs )
                cut_eff.passed('BDRS Mass', ra_temp.shape[0])


            #user cuts
            for f_cut in self.l_f_cuts:
                ra_temp = f_cut(ra_temp)

            if len(self.l_f_cuts) > 0:
                cut_eff.passed('User Cuts', ra_temp.shape[0])                

            ra_list.append(ra_temp)
            size_so_far += ra_temp.shape[0]

            if size_so_far >= self.end >= 1:
                break

        self._ra = np.concatenate(ra_list)
        self._do_range_math()
        self._ra = self._ra[self.start:self.end]

        cut_eff.summary()
        print 'Returning %d elements' % self._ra.shape[0]
        return self._ra



_jbase_ = '/u/eb/joshgc/mynfs/CSJets/logs/'
_ebase_ = '/a/sulky51/atlaswork.u1/e/estrauss/WTagger/npy/'
_ebase2_ = '/u/at/estrauss/atlint02/npy_pythia/'
_ebase_smear_ = '/a/sulky51/atlaswork.u1/e/estrauss/WTagger/npy_smear/'
_ebase_madgraph_h = '/u/at/estrauss/atlint02/npy_madgraph_from_scratch_v2/'
_ebase_madgraph_g = '/u/at/estrauss/atlint02/npy_madgraph_from_scratch_v2/'

register_us = [
    (_jbase_ + 'test_w_jets_no_rot.npy', 'test_w_jets_no_rot', 'test_w_jets_no_rot', 'w'),
    (_jbase_ + 'test_w_jets_no_ref_subjets.npy', 'test_w_jets_no_ref_subjets', 'test_w_jets_no_ref_subjets', 'w'),
    (_jbase_ + 'test_w_jets_no_ref_pa.npy', 'test_w_jets_no_ref_pa', 'test_w_jets_no_ref_pa', 'w'),

    (_jbase_ + 'tmp_fixed_w_0.npy', 'test_fixed_w_0', 'test_fixed_w_0', 'w'),
    (_jbase_ + 'tmp_fixed_v_0.npy', 'test_fixed_v_0', 'test_fixed_v_0', 'v'),
    (_jbase_ + 'tmp_fixed_w_1.npy', 'test_fixed_w_1', 'test_fixed_w_1', 'w'),
    (_jbase_ + 'tmp_fixed_v_1.npy', 'test_fixed_v_1', 'test_fixed_v_1', 'v'),
    (_jbase_ + 'tmp_fixed_w_2.npy', 'test_fixed_w_2', 'test_fixed_w_2', 'w'),
    (_jbase_ + 'tmp_fixed_v_2.npy', 'test_fixed_v_2', 'test_fixed_v_2', 'v'),

    (_jbase_ + 'tmp_fixed_w_0_c.npy','test_fixed_w_0_c','test_fixed_w_0_c','w'),
    (_jbase_ + 'tmp_fixed_v_0_c.npy','test_fixed_v_0_c','test_fixed_v_0_c','v'),

    (_jbase_ + 'tmp_fixed_w_l_c.npy','test_fixed_w_l_c','test_fixed_w_l_c','w'),
    (_jbase_ + 'tmp_fixed_v_l_c.npy','test_fixed_v_l_c','test_fixed_v_l_c','v'),


    (_jbase_ + 'tmp_fix_w_j_sj0.npy','test_fixed_w_j_sj0','test_fixed_w_j_sj0','w'),
    (_jbase_ + 'tmp_fix_v_j_sj0.npy','test_fixed_v_j_sj0','test_fixed_v_j_sj0','v'),

    (_jbase_ + 'test_std_w_pt.npy','test_std_w_pt','test_std_w_pt','w'),
    (_jbase_ + 'test_std_v_pt.npy','test_std_v_pt','test_std_v_pt','v'),

    (_ebase_smear_+'bsub_[0-9]*_p_v_1.2_200_0_T1_*.npy'  ,  'p_v_1.2_200_0_T1_S1'  , 'Pythia Light Smeared w/o PU' , 'v') ,
    (_ebase_smear_+'bsub_[0-9]*_p_v_1.2_200_30_T1_*.npy' ,  'p_v_1.2_200_30_T1_S1' , 'Pythia Light Smeared w/ PU'  , 'v') ,
    (_ebase_smear_+'bsub_[0-9]*_p_v_1.2_500_0_T1_*.npy'  ,  'p_v_1.2_500_0_T1_S1'  , 'Pythia Light Smeared w/o PU' , 'v') ,
    (_ebase_smear_+'bsub_[0-9]*_p_v_1.2_500_30_T1_*.npy' ,  'p_v_1.2_500_30_T1_S1' , 'Pythia Light Smeared w/ PU'  , 'v') ,

    (_ebase_smear_+'bsub_[0-9]*_p_w_1.2_200_0_T1_*.npy'  ,  'p_w_1.2_200_0_T1_S1'  , 'Pythia W Smeared w/o PU'     , 'w') ,
    (_ebase_smear_+'bsub_[0-9]*_p_w_1.2_200_30_T1_*.npy' ,  'p_w_1.2_200_30_T1_S1' , 'Pythia W Smeared w/ PU'      , 'w') ,
    (_ebase_smear_+'bsub_[0-9]*_p_w_1.2_500_0_T1_*.npy'  ,  'p_w_1.2_500_0_T1_S1'  , 'Pythia W Smeared w/o PU'     , 'w') ,
    (_ebase_smear_+'bsub_[0-9]*_p_w_1.2_500_30_T1_*.npy' ,  'p_w_1.2_500_30_T1_S1' , 'Pythia W Smeared w/ PU'      , 'w') ,
    
    (_ebase_+'bsub_[0-9]*_p_v_1.2_200_0_T0_*.npy' ,  'p_v_1.2_200_0_T0' , 'Pythia Light w/o PU (No Trimming)', 'v'),
    (_ebase_+'bsub_[0-9]*_p_v_1.2_500_0_T0_*.npy' ,  'p_v_1.2_500_0_T0' , 'Pythia Light w/o PU (No Trimming)', 'v'),
    (_ebase_+'bsub_[0-9]*_p_v_1.2_200_30_T0_*.npy',  'p_v_1.2_200_30_T0', 'Pythia Light w/ PU (No Trimming)' , 'v'),
    (_ebase_+'bsub_[0-9]*_p_v_1.2_500_30_T0_*.npy',  'p_v_1.2_500_30_T0', 'Pythia Light w/ PU (No Trimming)' , 'v'),

    (_ebase_+'bsub_[0-9]*_p_w_1.2_200_0_T0_*.npy' ,  'p_w_1.2_200_0_T0' , 'Pythia W w/o PU (No Trimming)', 'w'),
    (_ebase_+'bsub_[0-9]*_p_w_1.2_500_0_T0_*.npy' ,  'p_w_1.2_500_0_T0' , 'Pythia W w/o PU (No Trimming)', 'w'),
    (_ebase_+'bsub_[0-9]*_p_w_1.2_200_30_T0_*.npy',  'p_w_1.2_200_30_T0', 'Pythia W w/ PU (No Trimming)' , 'w'),
    (_ebase_+'bsub_[0-9]*_p_w_1.2_500_30_T0_*.npy',  'p_w_1.2_500_30_T0', 'Pythia W w/ PU (No Trimming)' , 'w'),

    (_ebase_+'bsub_[0-9]*_p_v_1.2_200_0_T1_*.npy' ,  'p_v_1.2_200_0_T1' , 'Pythia Light w/o PU', 'v'),
    (_ebase2_+'bsub_[0-9]*_p_v_1.2_250_0_T1_*.npy' ,  'p_v_1.2_250_0_T1' , 'Pythia Light w/o PU', 'v'),
    (_ebase2_+'bsub_[0-9]*_p_v_1.2_300_0_T1_*.npy' ,  'p_v_1.2_300_0_T1' , 'Pythia Light w/o PU', 'v'),
    (_ebase2_+'bsub_[0-9]*_p_v_1.2_350_0_T1_*.npy' ,  'p_v_1.2_350_0_T1' , 'Pythia Light w/o PU', 'v'),
    (_ebase2_+'bsub_[0-9]*_p_v_1.2_400_0_T1_*.npy' ,  'p_v_1.2_400_0_T1' , 'Pythia Light w/o PU', 'v'),
    (_ebase2_+'bsub_[0-9]*_p_v_1.2_450_0_T1_*.npy' ,  'p_v_1.2_450_0_T1' , 'Pythia Light w/o PU', 'v'),
    (_ebase_+'bsub_[0-9]*_p_v_1.2_500_0_T1_*.npy' ,  'p_v_1.2_500_0_T1' , 'Pythia Light w/o PU', 'v'),
    
    (_ebase_+'bsub_[0-9]*_p_v_1.2_200_30_T1_*.npy',  'p_v_1.2_200_30_T1', 'Pythia Light w/ PU' , 'v'),
    (_ebase_+'bsub_[0-9]*_p_v_1.2_500_30_T1_*.npy',  'p_v_1.2_500_30_T1', 'Pythia Light w/ PU' , 'v'),

    (_ebase_+'bsub_[0-9]*_p_w_1.2_200_0_T1_*.npy' ,  'p_w_1.2_200_0_T1' , 'Pythia W w/o PU', 'w'),
    (_ebase2_+'bsub_[0-9]*_p_w_1.2_250_0_T1_*.npy' ,  'p_w_1.2_250_0_T1' , 'Pythia W w/o PU', 'w'),
    (_ebase2_+'bsub_[0-9]*_p_w_1.2_300_0_T1_*.npy' ,  'p_w_1.2_300_0_T1' , 'Pythia W w/o PU', 'w'),
    (_ebase2_+'bsub_[0-9]*_p_w_1.2_350_0_T1_*.npy' ,  'p_w_1.2_350_0_T1' , 'Pythia W w/o PU', 'w'),
    (_ebase2_+'bsub_[0-9]*_p_w_1.2_400_0_T1_*.npy' ,  'p_w_1.2_400_0_T1' , 'Pythia W w/o PU', 'w'),
    (_ebase2_+'bsub_[0-9]*_p_w_1.2_450_0_T1_*.npy' ,  'p_w_1.2_450_0_T1' , 'Pythia W w/o PU', 'w'),
    (_ebase_+'bsub_[0-9]*_p_w_1.2_500_0_T1_*.npy' ,  'p_w_1.2_500_0_T1' , 'Pythia W w/o PU', 'w'),

    (_ebase_+'bsub_[0-9]*_p_w_1.2_200_30_T1_*.npy',  'p_w_1.2_200_30_T1', 'Pythia W w/ PU' , 'w'),
    (_ebase_+'bsub_[0-9]*_p_w_1.2_500_30_T1_*.npy',  'p_w_1.2_500_30_T1', 'Pythia W w/ PU' , 'w'),



    (_ebase_+'bsub_[0-9]*_h_v_1.2_200_0_T0_*.npy' ,  'h_v_1.2_200_0_T0' , 'Herwig Light w/o PU (No Trimming)', 'v'),
    (_ebase_+'bsub_[0-9]*_h_v_1.2_500_0_T0_*.npy' ,  'h_v_1.2_500_0_T0' , 'Herwig Light w/o PU (No Trimming)', 'v'),

    (_ebase_+'bsub_[0-9]*_h_w_1.2_200_0_T0_*.npy' ,  'h_w_1.2_200_0_T0' , 'Herwig W w/o PU (No Trimming)', 'w'),
    (_ebase_+'bsub_[0-9]*_h_w_1.2_500_0_T0_*.npy' ,  'h_w_1.2_500_0_T0' , 'Herwig W w/o PU (No Trimming)', 'w'),

    (_ebase_+'bsub_[0-9]*_h_v_1.2_200_0_T1_*.npy' ,  'h_v_1.2_200_0_T1' , 'Herwig Light w/o PU', 'v'),
    (_ebase_+'bsub_[0-9]*_h_v_1.2_500_0_T1_*.npy' ,  'h_v_1.2_500_0_T1' , 'Herwig Light w/o PU', 'v'),

    (_ebase_+'bsub_[0-9]*_h_w_1.2_200_0_T1_*.npy' ,  'h_w_1.2_200_0_T1' , 'Herwig W w/o PU', 'w'),
    (_ebase_+'bsub_[0-9]*_h_w_1.2_500_0_T1_*.npy' ,  'h_w_1.2_500_0_T1' , 'Herwig W w/o PU', 'w'),


    (_ebase_madgraph_g+'bsub_[0-9]*_m_g_1.2_300_0_T1_*.npy',   'm_g_1.2_300_0_T1' , 'Madgraph Light w/o PU', 'g'),
    (_ebase_madgraph_h+'bsub_[0-9]*_m_h_1.2_300_0_T1_*.npy',   'm_h_1.2_300_0_T1' , 'Madgraph Higgs w/o PU', 'h'),
    (_ebase_madgraph_g+'bsub_[0-9]*_m_g_1.2_500_0_T1_*.npy',   'm_g_1.2_500_0_T1' , 'Madgraph Light w/o PU', 'g'),
    (_ebase_madgraph_h+'bsub_[0-9]*_m_h_1.2_500_0_T1_*.npy',   'm_h_1.2_500_0_T1' , 'Madgraph Higgs w/o PU', 'h'),


    (_ebase2_+'bsub_[0-9]*_p_g_1.2_200_0_T1_*.npy' ,  'p_g_1.2_200_0_T1' , 'Pythia Gluons w/o PU', 'g'),
    (_ebase2_+'bsub_[0-9]*_p_q_1.2_200_0_T1_*.npy' ,  'p_q_1.2_200_0_T1' , 'Pythia Quarks w/o PU', 'q'),
    (_ebase2_+'bsub_[0-9]*_p_g_0.4_200_0_T1_*.npy' ,  'p_g_0.4_200_0_T1' , 'Pythia Gluons w/o PU', 'g'),
    (_ebase2_+'bsub_[0-9]*_p_q_0.4_200_0_T1_*.npy' ,  'p_q_0.4_200_0_T1' , 'Pythia Quarks w/o PU', 'q'),


            ]



_ds_ = {}
for path, key, name, me in register_us:
    _ds_[key] = dataset(path, key, name, _di_pdg_rng_[me])

def make_prop_str(generator, metype, rad, ptbin, mu, trim, smear):
    key = '%s_%s_%.1f_%d_%d_T%d' % (generator[0], metype, rad, ptbin, mu, trim)
    if smear > 0:
        key += '_S%d' % smear

    return key

_ps_re_  = r'(?P<gen>[^_]+)_(?P<me>[^_]+)_(?P<rad>[0-9.e+-]+)_(?P<ptbin>[0-9.e+-]+)'
_ps_re_ += r'_(?P<mu>[0-9e+-]+)_T(?P<trim>[0-9e+-]+)(?P<smear>_S[0-9.e+-]+)?'
def parse_prop_str(ss):
    import re
    mm = re.search(_ps_re_, ss)
    if mm is None:
        return {'misc': ss}

    di = mm.groupdict()
    if di['smear'] is None:
        di['smear'] = '0'
    else:
        di['smear'] = di['smear'][2:]

    di['misc'] = ''
    return di


_pk_re_ = r'(?P<ds_a>[^-]+)-vs-(?P<ds_b>[^-]+)-(?P<extra_rad>[0-9.]+_)?'
_pk_re_+= r'((?P<min_m>[0-9.]+)m(?P<max_m>[0-9.]+))?_'
_pk_re_+= r'(((?P<min_bdrs>[0-9.-]+)bdrs(?P<max_bdrs>[0-9.-]+))?_)?'
_pk_re_+= r'((?P<min_dr>[0-9.]+)d(?P<max_dr>[0-9.]+))?_'
_pk_re_+= r'(?P<start>[0-9.-]+)_(?P<end>[0-9.-]+)$'
def parse_pickle_name(inname):
    import os
    import re

    ss = os.path.splitext(inname)[0]
    mm = re.search(_pk_re_, ss)
    if mm is None:
        return {'misc': ss}

    di = mm.groupdict()
    del di['extra_rad']
    di['misc'] = ''
    return di


def diff_prop_strs(nested_list):
    if all([isinstance(ob, str) for ob in nested_list]):
        return diff_dicts([parse_prop_str(ob) for ob in nested_list])

def make_verbose_title(prop_di):
    if isinstance(prop_di, str):
        prop_di = parse_prop_str(prop_di)

    for kk in ['gen', 'rad', 'trim', 'me', 'mu', 'smear', 'ptbin']:
        if kk not in prop_di:
            prop_di[kk] = ''
        else:
            try:
                prop_di[kk] = pretty[kk][prop_di[kk]]
            except KeyError:
                pass

    jet_desc = ['ptbin', 'rad', 'trim', 'me']
    env_desc = ['gen', 'mu', 'smear']
    jet_desc = ', '.join(prop_di[kk] for kk in jet_desc if prop_di[kk] != '')
    env_desc = ', '.join(prop_di[kk] for kk in env_desc if prop_di[kk] != '')

    if jet_desc == '' and env_desc == '' and 'misc' in prop_di:
        return prop_di['misc'].replace('_', '\\_')

    if jet_desc == '': return env_desc
    if env_desc == '': return jet_desc + ' Jets'
    return jet_desc + ' Jets: ' + env_desc


def diff_dicts(all_dicts):
    cc = dict()
    li_dicts = [{} for ss in all_dicts]
    comm_keys = set(all_dicts[0].keys())
    for di in all_dicts[1:]:
        comm_keys = comm_keys & set(di.keys())
    for kk in comm_keys:
        vv = all_dicts[0][kk]
        if all(di[kk] == vv for di in all_dicts):
            cc[kk] = vv
        else:
            for ii in range(len(all_dicts)):
                li_dicts[ii][kk] = all_dicts[ii][kk]

    for ii, di in enumerate(all_dicts):
        for kk, vv in di.iteritems():
            if kk not in comm_keys:
                li_dicts[ii][kk] = vv

    return cc, li_dicts

def get_dataset(generator, metype, ptbin, mu, rad, trim, smear=0, pdg=True,
                start=0, end=-1, 
                min_m=-1, max_m=-1, 
                min_bdrs=-1, max_bdrs=-1, 
                min_dr=-1, max_dr=-1,
                min_pt=-1, max_pt=-1,
                n_b=-1, subjet_e=-1
                ):
    key = make_prop_str(generator, metype, rad, ptbin, mu, trim, smear)
    return get_exact_dataset(key, start, end, min_m, max_m, min_bdrs, max_bdrs, min_dr, max_dr, min_pt, max_pt, pdg, n_b, subjet_e)

def get_exact_dataset(key, start=0, end=-1, min_m=-1, max_m=-1, min_bdrs=-1, max_bdrs=-1, min_dr=-1, max_dr=-1, min_pt=-1, max_pt=-1, pdg=True, n_b=-1, subjet_e=-1):
    import copy
    ds = copy.copy(_ds_[key])
    ds.set_range(start, end)
    ds.set_mass_range(min_m , max_m )
    ds.set_bdrs_range(min_bdrs , max_bdrs )
    ds.set_dr_range  (min_dr, max_dr)
    ds.set_pt_range(min_pt, max_pt)
    ds.set_n_b(n_b)
    ds.set_subjet_e(subjet_e)
    if not pdg:
        ds.s_pdgid = set()
    return ds

if __name__ == '__main__':
    #just toy code for testing

    keyl = make_prop_str('p', 'v', 1.2, 200, 30, 1, 1)
    keyr = make_prop_str('p', 'w', 1.2, 200, 30, 1, 1)
    cc, li_diffs = diff_prop_strs([keyl, keyr])
    print make_verbose_title(cc)
    print make_verbose_title(li_diffs[0])
    print make_verbose_title(li_diffs[1])
    #ds = get_exact_dataset('pythia_w_500_30_1.2', 0, 100).get_array()
    #print ds.shape
    #print ds['cells'].shape
    #print ds['cells'].reshape(ds.shape[0], 625).shape
