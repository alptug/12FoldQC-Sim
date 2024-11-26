# Open Boundary 12-fold Symmetric Quasicrystal Simulator

This repository contains the simulation code used in the paper:  
**“Defects Enhance Stability in 12-fold Symmetric Soft-Matter Quasicrystals.”**

## Supported Operating Systems
The code is designed to run on **Linux** and **macOS**.

## Compilation Instructions
To compile the code, use the following command:
```bash
gcc -o rhombi -Wall -std=c99 -O3 rhombi.c -lm
```
## Usage

To run the simulation, use:
```bash
./rhombi <st_D> <st_I> <gamma> <eps_rh>
```

### Parameters

1.	**st_D**: Parameter of the randomized $[D,I]$-Stampfli inflationary tiling used as the initial configuration.
2.	**st_I**: Parameter of the randomized $[D,I]$-Stampfli inflationary tiling.
3.	**gamma**: Line tension value in $\beta\gamma a$.
4.	**eps_rh**: Energy cost of adding a rhombus to the system ($\beta\epsilon_R$).

## Output

The code generates a directory in the same path as the executable, named:
```
data_size_<st_D>_<st_I>_gamma_<gamma>_epsRh_<eps_rh>
```

Inside this directory, two files are created:

1.	state.csv: Contains details of the state point.
2.	measurements.csv: Contains simulation measurements.

### state.csv Format

This file contains the following columns:

1.	Number of vertices
2.	$D$ parameter of Stampfli tiling
3.	$I$ parameter of Stampfli tiling
4.	Line tension ($\beta\gamma a$)
5.	Energy cost of adding a triangle ($\beta\epsilon_{Tr}$)
6.	Energy cost of adding a square ($\beta\epsilon_{Sq}$)
7.	Energy cost of adding a rhombus ($\beta\epsilon_R$)

### measurements.csv Format

This file contains the following columns:

1.	Monte Carlo Sweep
2.	Number of squares
3.	Number of triangles
4.	Number of rhombi
5.	Boundary length
6.	Number of directed edges in each direction ($e_1, e_2, \ldots, e_{12}$)

## Citation

If you use this code, please cite the associated paper:
“Defects Enhance Stability in 12-fold Symmetric Soft-Matter Quasicrystals”
```
@article{ulugol2024defects,
  title={Defects Enhance Stability in 12-fold Symmetric Soft-Matter Quasicrystals},
  author={Ulug{\"o}l, Alptu{\u{g}} and Hardeman, Robert J and Smallenburg, Frank and Filion, Laura},
  journal={arXiv preprint arXiv:2408.08168},
  year={2024}
}
```

## Contact

For questions or issues, please contact:

Alptuğ Ulugöl

[e-mail address is stated in the paper]
