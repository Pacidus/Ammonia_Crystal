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
K_POINTS gamma
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
    "wfc": 1e1,
    "econv": 1e-9,
    "beta": 0.5,
}


def CV_input(file, **kwargs):
    """Save the input and run the test"""
    for k in kwargs:
        kwargs[k] = f"{kwargs[k]:12.10e}".replace("e", "d")

    with open(file, "w") as f:
        f.write(Input.format(**kwargs))

    with open(file.replace("in", "out"), "w") as f:
        run(["pw.x", "-i", file], stdout=f)


dt = [1, 60, 3600]
tr = {ord("h"): " ", ord("m"): " ", ord("s"): ""}


def sec(time):
    times = time.translate(tr).split()[::-1]
    return sum([float(i) * j for i, j in zip(times, dt)])


if __name__ == "__main__":
    Fermi = []
    Energy = []
    TimeCPU = []
    TimWall = []

    Rho = np.logspace(1, 3, 10)
    for i, rho in enumerate(Rho):
        vals["rho"] = rho
        CV_input("input.in", **vals)

        with open("output.out", "r") as f:
            lines = f.readlines()

        for line in lines:
            match line.split():
                case ["the", "Fermi", "energy", "is", fermi, "ev"]:
                    Fermi.append(float(fermi))
                case ["!", "total", "energy", "=", energy, "Ry"]:
                    Energy.append(float(energy))
                case ["PWSCF", ":", tcpu, "CPU", twall, "WALL"]:
                    TimeCPU.append(sec(tcpu))
                    TimWall.append(sec(twall))

    np.savetxt("conv.csv", np.array([Rho, Energy, Fermi, TimWall, TimeCPU]).T)
