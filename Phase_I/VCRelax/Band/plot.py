import numpy as np
import matplotlib.pyplot as plt
from plot_band import read_fermi, read_bnd, read_Hsym, plot

fermi = read_fermi("Ammonia.out")
bands, x = read_bnd("AC_bands.dat")
Xsym = read_Hsym("CAPI_band.out")

plot(bands, x, fermi)
plt.xticks(ticks= Xsym, labels=['L', '$\Gamma$', 'X', 'U', '$\Gamma$'])
plt.ylabel("Energy (eV)")
plt.show()
