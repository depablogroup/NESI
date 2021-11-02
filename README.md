# NESI: Non-Equilibrium Simulations of Interfaces

NESI is a simple C++ program for conducting nonequilibrium molecular dynamics (NEMD) simulation for liquid-vapor interfaces.
This code was used to generate the data in our forthcoming paper on interfacial thermodynamics [1].

The standard system contains walls at the top and bottom of the simulation cell which are either attractive or reflecting, ensuring
that the adjoining phases are either liquid or vapor, respectively. Near these walls, thermodynamic baths are implemented which
maintain constant temperature and chemical potential difference between components; for now only pure and binary systems are supported.

The code is capable of calculating a wide variety of spatially-resolved physical and thermodynamic properties on-the-fly. 
These include kinetic and potential energies, total and species momenta, chemical potentials, mass and species mass densities, stresses
and pressures, energy and mass fluxes, and more.

## Compiling

To compile for single-core use:
```
make
```

For OpenMP support:
```
make omp
```

## Usage

The program looks for a file called "input" which provides run instructions. An example input file is provided with comments.
To run the code, simple call `./nesi` or `./nesi_omp` depending on how you compiled (see above).

## License

NESI is available under the MIT license, included in this repository.

### References

[1] Rauscher, Oettinger, and de Pablo, "Nonequilibrium Statistical Thermodynamics of Multicomponent Interfaces," Submitted (2021).

### Author

Phil Rauscher, phillip.m.rauscher@gmail.com




