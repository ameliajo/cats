import numpy as np
ENDF_y = [0, .0066, .0264, .0594, .1055, .1649, .2374, .3232, .4221, .5342, .6595, .7980, .9497, 1.1146, 1.2927, 1.4839, 1.6884, 2.0169, 2.4373, 2.9366, 3.6133, 4.6775, 7.1346, 7.3650, 7.5156, 7.6733, 7.8309, 8.0740, 8.4419, 9.0595, 9.6773, 7.3645, 6.2674, 5.1965, 4.7958, 4.8024, 4.6841, 4.4673, 4.1914, 3.8169, 3.3439, 2.7855, 3.2782, 5.3082, 8.5930, 12.3377, 8.4616, 5.6695, 4.1585, 2.6081, 0.0] 
ENDF_x = [0.0008*i for i in range(len(ENDF_y))]
invArea = 1.0/(np.trapz(ENDF_y,x=ENDF_x))
ENDF_y = [invArea * y for y in ENDF_y]

if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    plt.title('ENDF-B VIII.0 013-Al Phonon Distribution')
    plt.plot(ENDF_x,ENDF_y)
    plt.show()
