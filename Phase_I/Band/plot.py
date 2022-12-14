import numpy as np
import matplotlib.pyplot as plt
from plot_band import read_fermi, read_bnd, plot

fermi = read_fermi("Ammonia.out")
bands, x = read_bnd("AC_bands.dat")

plot(bands, x, fermi)
plt.show()
