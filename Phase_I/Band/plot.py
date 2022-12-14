import numpy as np
import matplotlib.pyplot as plt
from plot_band import read_fermi, read_bnd, plot

fermi = read_fermi("Ammonia.out")
bands, x = read_bnd("AC_bands.dat")

plot(bands, x, fermi)
plt.xticks(ticks= [0.0000, 0.5000, 1.2071, 1.5133, 2.0996], labels=['L', '$\Gamma$', 'X', 'U', '$\Gamma$'])
plt.ylabel("Energy (eV)")
plt.show()
