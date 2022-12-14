# Phase I
La phase I du cristal d'ammoniaque ce présente comme un cristal cubique.
Sa symétrie est P 21 3 (198), on utilise VESTA pour pouvoir initialiser les position des atomes dans le cristal en utilisant les ressources :

- [Liste des logiciels fonctionnant avec Quantum Espresso](https://www.quantum-espresso.org/auxiliary-software/)
- Amonia_project.pdf
- Simulating_thermal_motion_in_crystalline_phase-I_a.pdf

## PWscf 
On utilise l'outil "qe input generator" disponible en ligne :

- [https://www.materialscloud.org/work/tools/qeinputgenerator](https://www.materialscloud.org/work/tools/qeinputgenerator)

## Band Structure
À l'aide de différentes ressources nous pouvons établir une méthode pour obtenir des structures de bandes :

- [https://www.jappoker.com/blog/2019/band-diagram-QE/](https://www.jappoker.com/blog/2019/band-diagram-QE/)
- [https://pranabdas.github.io/espresso/](https://pranabdas.github.io/espresso/)

## Relaxation
Nous avons réussis à faire différentes relaxations, le paramètre

- [Inputs correspondant à la pression](https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1105)

