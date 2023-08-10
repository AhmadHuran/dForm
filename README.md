# dForm
A cli tool for recasting GTH/HGH pseudo-potentials into diagonal-projection form
$$\sum_{lm} \sum_{ij} |p_i^{lm} \rangle h^l_{ij}\langle p_j^{lm}|\qquad\implies\qquad \sum_{nlm}|\beta_n^{lm} \rangle e^l_{n}\langle \beta_n^{lm}|$$

The radial projectors written eventually to the output file are multiplied by the radial coordinate i.e. $r \beta_n^{l}\(r\)$ are written.

## Installation 
#### Installing cargo
If you install rust from [rustup](https://rustup.rs/), cargo will be installed automatically. 
#### Installing dForm
dForm is installed through cargo's cli
```
cargo install --root <path to installation root> --git https://github.com/AhmadHuran/dForm
```
The binary `dform` will be installed at `<path to installation root>/bin/`
## Usage
#### Help 
```
dform --help
```
#### Example
Download the cp2k GTH potentials database [file](https://github.com/cp2k/cp2k/blob/master/data/GTH_POTENTIALS) from the [cp2k](https://github.com/cp2k/cp2k) repository and execute
```
dform database <path to download>/GTH_POTENTIALS --element=Na --xc=LDA --qion=1 --output=Na_LDA_q1.psp8
```
The file `Na_LDA_q1.psp8` is created and the following is printed to the console
```
[INFO] [database mode]: Attempt reading potential 'Na-LDA-q1' from file: '<path to download>/GTH_POTENTIALS'
[INFO] [database mode]: Success!
[INFO] [database mode]: Attempt writing potential 'Na-LDA-q1' to file: 'Na_LDA_q1.psp8'
[INFO] [database mode]: Success!
```
## Testing
Single-potential files in both [cp2k](https://www.cp2k.org/) and [abinit](https://www.abinit.org/) formats can be downloaded [here](https://htmlpreview.github.io/?https://github.com/cp2k/cp2k-data/blob/master/potentials/Goedecker/index.html).
Abinit can then be used to contrast calculations using the abinit-formatted potential and the `psp8` file created by `dform` from the downloaded potential in cp2k format. The `psp8` file is created in the `single` mode for less verbose input
```
dform single <path to cp2k-formatted potential> --output=some.psp8
```
If necessary the radial grid can be made finer using e.g. `--sep=0.005`. The default separation is 0.01 Bohr.

## References
#### Transformation
- [P. E. Blöchl, Generalized separable potentials for electronic-structure calculations, Phys. Rev. B 41(8), 5414 (1990)](https://www.doi.org/10.1103/PhysRevB.41.5414)
#### GTH/HGH pseudo-potentials
- [S. Goedecker, M. Teter, and J. Hutter, Separable dual-space Gaussian pseudopotentials, Phys. Rev. B 54, 1703 (1996)](https://www.doi.org/10.1103/PhysRevB.54.1703)
- [C. Hartwigsen, S. Goedecker, and J. Hutter, Relativistic separable dual-space Gaussian pseudopotentials from H to Rn, Phys. Rev. B 58, 3641 (1998)](https://www.doi.org/10.1103/PhysRevB.58.3641)
- [M. Krack, Pseudopotentials for H to Kr optimized for gradient-corrected exchange-correlation functionals, Theor. Chem. Acc. 114, 145 (2005)](https://www.doi.org/10.1007/s00214-005-0655-y)
