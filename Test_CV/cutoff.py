import numpy as np
from subprocess import run

Input = """&CONTROL
  calculation = 'scf'
  outdir = './out/'
  prefix = 'AC_Ph_1'
  pseudo_dir = './pseudo/'
  tprnfor = .true.
  tstress = .true.
  verbosity = 'low'
  dt = 1.d2
  etot_conv_thr =   {ethr}
  forc_conv_thr =   {fthr}
/
&SYSTEM
  degauss =   {dgauss}
  ecutrho =   {rho}
  ecutwfc =   {wfc}
  ibrav = 0
  nat = 16
  nosym = .false.
  ntyp = 2
  occupations = 'smearing'
  smearing = 'cold'
/
&ELECTRONS
  conv_thr =   {econv}
  electron_maxstep = 80
  mixing_beta =   {beta}
/
ATOMIC_SPECIES
H      1.00794 H.pbe-rrkjus_psl.1.0.0.UPF
N      14.0067 N.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
N            0.2110000000       0.2110000000       0.2110000000
N            0.2890000000       0.7890000000       0.7110000000
N            0.7890000000       0.7110000000       0.2890000000
N            0.7110000000       0.2890000000       0.7890000000
H            0.3900000000       0.2430000000       0.1240000000
H            0.1100000000       0.7570000000       0.6240000000
H            0.6100000000       0.7430000000       0.3760000000
H            0.8900000000       0.2570000000       0.8760000000
H            0.1240000000       0.3900000000       0.2430000000
H            0.6240000000       0.1100000000       0.7570000000
H            0.3760000000       0.6100000000       0.7430000000
H            0.8760000000       0.8900000000       0.2570000000
H            0.2430000000       0.1240000000       0.3900000000
H            0.7570000000       0.6240000000       0.1100000000
H            0.7430000000       0.3760000000       0.6100000000
H            0.2570000000       0.8760000000       0.8900000000
K_POINTS automatic
7 7 7 0 0 0
CELL_PARAMETERS angstrom
      5.1850000000       0.0000000000       0.0000000000
      0.0000000000       5.1850000000       0.0000000000
      0.0000000000       0.0000000000       5.1850000000
"""

vals = {
    "ethr": 1e-5,
    "fthr": 1e-5,
    "dgauss": 1.4699723600e-2,
    "rho": 4.8e2,
    "wfc": 6e1,
    "econv": 1e-9,
    "beta": 0.5,
}


def CV_input(file, **kwargs):
    """Save the input and run the test"""
    for k in kwargs:
        kwargs[k] = f"{kwargs[k]:12.10e}".replace("e", "d")

    with open(file, "w") as f:
        f.write(Input.format(**kwargs))

    run(["pw.x", "-i", file], stdout=file.replace("in", "out"))


if __name__ == "__main__":

    for i, rho in enumerate(np.logspace(1, np.log10(1000))):
        vals["rho"] = rho
        CV_input(f"file_{i}.in", **vals)
