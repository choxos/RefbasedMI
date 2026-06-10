# RefBasedMI: Reference-Based Imputation for Longitudinal Clinical Trials with Protocol Deviation

Imputation of missing numerical outcomes for a longitudinal trial with
protocol deviations. The package uses distinct treatment arm-based
assumptions for the unobserved data, following the general algorithm of
Carpenter, Roger, and Kenward (2013)
[doi:10.1080/10543406.2013.834911](https://doi.org/10.1080/10543406.2013.834911)
, and the causal model of White, Royes and Best (2020)
[doi:10.1080/10543406.2019.1684308](https://doi.org/10.1080/10543406.2019.1684308)
. Sensitivity analyses to departures from these assumptions can be done
by the Delta method of Roger. The program uses the same algorithm as the
'mimix' 'Stata' package written by Suzie Cro, with additional coding for
the causal model and delta method. The reference-based methods are jump
to reference (J2R), copy increments in reference (CIR), copy reference
(CR), and the causal model, all of which must specify the reference
treatment arm. Other methods are missing at random (MAR) and the last
mean carried forward (LMCF).

## See also

Useful links:

- <https://github.com/choxos/RefbasedMI>

- <https://choxos.github.io/RefbasedMI/>

- Report bugs at <https://github.com/choxos/RefbasedMI/issues>

## Author

**Maintainer**: Matteo Quartagno <m.quartagno@ucl.ac.uk>

Authors:

- Kevin McGrath <kkevinmcgrath@yahoo.co.uk>

- Ahmad Sofi-Mahmudi <a.sofimahmudi@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-6829-0823))
