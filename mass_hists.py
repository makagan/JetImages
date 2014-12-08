import numpy as np
import matplotlib.pyplot as plt


def draw(nameAndLists, limits, fname):
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    for name, li in nameAndLists:
        if '30' in name:
            if 'Untrim' in name:
                cc = 'b00000'
            else:
                cc = 'ff7272'
        else:
            if 'Untrim' in name:
                cc = '0000b0'
            else:
                cc = '7272ff'
        cc = '#'+cc
        hh, edges = np.histogram(li, limits[0], limits[1:3])
        xs = (edges[0:-1] + edges[1:])/2.0
        #ax.plot(xs, hh, 'o', label=name)
        ax.errorbar(xs, hh, yerr=np.sqrt(hh), markeredgecolor=cc,markerfacecolor=cc,ecolor=cc, fmt='o', label=name)


    ax.set_xlabel("Jet Mass (GeV)")
    ax.set_ylabel("#Entries per %.1f GeV" % ((limits[2] - limits[1]+0.0)/limits[0]))
    ax.set_ylim(0, 1.2*ax.get_ylim()[1])
    plt.legend()
    plt.savefig(fname)
    #plt.show()



if __name__ == '__main__':
    inpath='/u/eb/joshgc/mynfs/CSJets/'
    nameAndLocations = [('Trimmed t mu=0'   , 't_0_mass_trim_1.mass' ),
                        ('Trimmed t mu=30'  , 't_30_mass_trim_1.mass'),
                        ('Untrimmed t mu=0' , 't_0_mass_trim_0.mass'),
                        ('Untrimmed t mu=30', 't_30_mass_trim_0.mass'),
                        ]
    outname = 't_trimming.png'
    limits = (40,50,250)    

    #nameAndLocations = [('Trimmed W mu=0'   , 'w_0_mass_trim_1.mass' ),
    #                    ('Trimmed W mu=30'  , 'w_30_mass_trim_1.mass'),
    #                    ('Untrimmed W mu=0' , 'w_0_mass_trim_0.mass'),
    #                    ('Untrimmed W mu=30', 'w_30_mass_trim_0.mass'),
    #                    ]

    #outname = 'w_trimming.png'
    #limits = (50, 0, 150)

    nameAndLocations = [('Trimmed g mu=0'   , 'g_0_mass_trim_1.mass' ),
                        ('Trimmed g mu=30'  , 'g_30_mass_trim_1.mass'),
                        ('Untrimmed g mu=0' , 'g_0_mass_trim_0.mass'),
                        ('Untrimmed g mu=30', 'g_30_mass_trim_0.mass'),
                        ]

    outname = 'g_trimming.png'
    limits = (50, 0, 150)

    nameAndLists = [(name, list(float(xx) for xx in open(inpath+loc)))
                    for name,loc in nameAndLocations]


    draw(nameAndLists, limits,outname)
