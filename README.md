# Thermogravity
Codes of subroutines to calculate the thermal diffusion and gravitational effects on liquid columns

By: Sugata Tan, Planetary Science Institute, 2022

The codes were written in FORTRAN and had been tested and run for the calculations using Compaq Visual Fortran © 2000 Professional Edition 6.6.0. The routines to calculate the fugacity coefficients and density, which depend on the users’ choice of EOS, as well as that calculates the inverse of a matrix, are not included, so they need to be provided by the users.

These codes are also published in Planetary and Space Science as the Supplementary Materials of this paper: Tan, S.P.; Adidharma, H. On the stability and phase behavior of Titan's subsurface liquid columns, Planetary and Space Science 214 (2022) 105451. https://doi.org/10.1016/j.pss.2022.105451

The codes were applied with PC-SAFT equation of state (not included here) to describe Titan's subsurface liquid columns, which are modeled as ternary mixtures of nitrogen/methane/ethane. With gravity and temperature gradient as input, the outputs are pressure, density, and compositional profiles versus depth. Details can be found in the above paper.
