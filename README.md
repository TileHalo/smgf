# SMGF ![build status](https://github.com/TileHalo/fnir/actions/workflows/rust.yml/badge.svg)
SMGF provides stratified (layered) media electromagnetic green's functions for 
microwave and antenna engineering problems.

## Usage
Add the following under `[dependencies]` in your cargo.toml

```toml
smgf = 0.1.0
```


WARNING: This project is currently a vaporware AND not uploaded to [crates.io](https://crates.io)

## Use cases
SMGF is meant to be used in modelling electromagnetic problems arising from
fields in layered media (ex. PCB antennas) and used when building impedance
matrices in Method of Moments applications.

Currently planned Green's functions and dyadics (GF) in spectral domain are
 - [ ] Electric field GF due to electric current
 - [ ] Magnetic field GF due to electric current
 - [ ] Electric field GF due to magnetic current
 - [ ] Magnetic field GF due to magnetic current
 - [ ] MPIE GF:s

The planned evaluations from spectral domain are:
 - [ ] Michalski-Mosig algorithm
 - [ ] DCIM with various techniques (ex. Double DCIM, Taguchi-gradient etc.)

In some distant future the following would also be nice:
 - [ ] Vectorization (i.e. SIMD)
 - [ ] Branch cut integration
 - [ ] Integration along steepest descent path

## Contributing
Feel free to submit bug reports, feature requests, and the best of all,
pull requests to those issues.
