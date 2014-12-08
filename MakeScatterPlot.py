import sys
from math import sqrt, cos, pi
import numpy as np
import matplotlib.spines
import matplotlib.pyplot    as plt
import matplotlib.patches   as ptc
from sklearn.decomposition import RandomizedPCA
from fisher import Fisher

#so that the data dont' keep shifting call to call
np.random.seed(123456789)


m1 = np.array([-0.5, 1])
m2 = np.array([0.5 ,-1])
cov1 = np.array([[1, 1.5],
                 [1.5, 2]])

xx = cos(pi/7)
rot  = np.array([[xx        , -1*sqrt(1-xx)],
                 [sqrt(1-xx), xx        ]])

m1 = np.dot(m1, rot)
m2 = np.dot(m2, rot)
cov1 = np.dot(cov1, rot)

d1 = np.random.multivariate_normal(m1, cov1, 100)
d1mean = np.sum(d1, axis=0) / (1.0*d1.shape[0])

#cov2 = np.array([[1, 1.5],[1.5,2]])
#d2 = np.random.multivariate_normal(m2, cov2, 100)
d2 = np.random.multivariate_normal(m2, cov1, 100)
d2mean = np.sum(d2, axis=0) / (1.0*d2.shape[0])

x_data_max = max(np.max(d1[:,0]), np.max(d2[:,0]))
y_data_max = max(np.max(d1[:,1]), np.max(d2[:,1]))
x_data_min = min(np.min(d1[:,0]), np.min(d2[:,0]))
y_data_min = min(np.min(d1[:,1]), np.min(d2[:,1]))


classmeandiff = d1mean-d2mean


alldata = np.concatenate( (d1, d2), axis=0)
alllabels = np.concatenate( ( np.array([0 for i in range(len(d1))]), np.array([1 for i in range(len(d2))]) ), axis=0)
meanall = np.sum(alldata, axis=0) / (1.0*alldata.shape[0])

pca = RandomizedPCA(2)
pca.fit( alldata )
d1pca = pca.transform(d1)
d2pca = pca.transform(d2)
mean1pca = pca.transform(d1mean)
pcasign = np.array([1.0, 1.0])
if mean1pca[0]<0: pcasign[0] = -1.0
if mean1pca[1]<0: pcasign[1] = -1.0


fish = Fisher(norm_covariance=True)
fish.fit_multiclass(alldata, alllabels, use_total_scatter=False, solution_norm="N", sigma_sqrd=1e-8, tol=5e-3, print_timing=False)
d1fish = fish.transform(d1)
d2fish = fish.transform(d2)
fishsign=1.0
mean1fish = fish.transform(d1mean)
if mean1fish[0]<0: fishsign = -1.0


c_sig = 'g'
c_bkg = 'r'
to_do = [(classmeandiff     , 'Class Mean Difference', '#04819e'),
         (pca.components_[0], 'First PCA Component'  , '#a65200'),
         (pca.components_[1], 'Second PCA Component' , '#ff7f00'),
         (fish.w_[0]        , 'Fisher Discriminant'  , '#1435ad'),
        ]

#for doVector in (False,):
for doVector in (True, False):
    for n_subs in range(3,4):
        fig = plt.figure(figsize=(3.7*2,3*2))
        ax = plt.subplot2grid((4,5), (0,0), colspan=4, rowspan=4)
        ax.scatter(d1[:,0], d1[:,1], c=c_sig, edgecolor=c_sig, lw=0)
        ax.scatter(d2[:,0], d2[:,1], c=c_bkg, edgecolor=c_bkg, lw=0)

        
        patches = [plt.Rectangle((1,1), 1, 1, ec='w', fc='w')]*4
        names   = ['']*4
        x = meanall[0]
        y = meanall[1]
        ll = 3

        #to be fair i should use both datasets.
        #Lifes not fair and neither am I
        max_dx = np.max(np.abs(d1[:,0] - meanall[0]))
        max_dy = np.max(np.abs(d1[:,1] - meanall[0]))
        for ii, li, in enumerate(to_do[:n_subs+1]):
            ra, name, color = li
            dx = ll * ra[0]/sqrt(ra[0] ** 2 + ra[1]**2)
            dy = ll * ra[1]/sqrt(ra[0] ** 2 + ra[1]**2)
            if doVector:
                if dy < 0:
                    dx *= -1
                    dy *= -1
                obj = ax.arrow(x, y, dx, dy, color=color,
                         head_width=0.20, head_length=0.2,
                         linewidth=3, label=name)
            else:
                dx, dy = -dy, dx
                factor = min(max_dy/abs(dy), max_dx/abs(dx))
                dx *= factor
                dy *= factor
                obj = ax.arrow(x-dx, y-dy, 2*dx, 2*dy, color=color,
                         head_width=0.001, head_length=0.001,
                         linewidth=3, label=name, linestyle='solid')

            patches[ii] = obj
            names  [ii] = name

        ax.legend(patches, names, loc=4, fontsize=14, 
                    frameon=False)
        ax.get_legend().set_title(
            'Discriminant Directions' if doVector else 'Discriminanting Planes',
            prop={'size':18})


        ax.set_title("Toy Data", fontsize=18)
        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_xlim(-6, 6)
        ax.set_ylim(-7, 5)
        #ax.set_axis_off()
        
        _axes = []
        _to_hist = [
                    (np.dot(d1, classmeandiff), np.dot(d2, classmeandiff), 'Class Means'),
                    (pcasign[0]*d1pca[:,0], pcasign[0]*d2pca[:,0], 'First PCA'),
                    (pcasign[1]*d1pca[:,1], pcasign[1]*d2pca[:,1], 'Second PCA'),
                    (fishsign*d1fish[:,0], fishsign*d2fish[:,0]  , 'Fisher'),
                    ]
        
        for ii, li in enumerate(_to_hist[:n_subs+1]):
            _axes.append(plt.subplot2grid((4,5), (ii, 4)))
            if ii==0:
                plt.title('Projections', fontsize=18)
            _axes[-1].hist(li[0], 10, normed=1, facecolor=c_sig, alpha=0.5)
            _axes[-1].hist(li[1], 10, normed=1, facecolor=c_bkg, alpha=0.5)
            _axes[-1].set_xticks([])
            _axes[-1].set_yticks([])
            _axes[-1].set_xlabel(li[2], size='medium')
            for child in _axes[-1].get_children():
                if isinstance(child, matplotlib.spines.Spine):
                    child.set_color(to_do[ii][2])

        plt.show()
        #plt.savefig('demo_linear_%s_%d.png' % ('directions' if doVector else 'planes', n_subs))
