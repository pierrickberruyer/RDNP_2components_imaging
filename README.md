# RDNP_2components_imaging
Matlab code for imaging the internal composition of core-shell particles by fitting relayed DNP (dynamic nuclear polarization) enhancement build-up curves via pattern search on spin-diffusion PDE solutions. Accompanies Berruyer et al., J. Am. Chem. Soc. 2023 (10.1021/jacs.3c01279); data at 10.5281/zenodo.7756308. License: CC BY-SA 4.0.

## Publication and citation

This code accompanies:

> Berruyer, P.; Cibaka-Ndaya, C.; Pinon, A.; Sanchez, C.; Drisko, G. L.; Emsley, L.
> *Imaging radial distribution functions of complex particles by relayed dynamic nuclear polarization.*
> **J. Am. Chem. Soc.** **2023**. https://doi.org/10.1021/jacs.3c01279

The archived version of record (raw NMR data, extraction scripts, and the
sample-specific fitting scripts `SampleA_3_RDNP_imaging.m` /
`SampleB_3_RDNP_imaging.m`) is available on Zenodo:
https://doi.org/10.5281/zenodo.7756308

**Any use of this code, original or modified, is required to cite the Zenodo
deposit (DOI 10.5281/zenodo.7756308). Citation of the original publication is
invited.**

## License

CC BY-SA 4.0 (Creative Commons Attribution-ShareAlike 4.0 International).
See the [LICENSE](LICENSE) file or
https://creativecommons.org/licenses/by-sa/4.0/.

The code is provided free of charge and comes with no warranty.

## Requirements

- Matlab (R2021a or later recommended)
- Optimization Toolbox (`patternsearch`)
- Partial Differential Equation Toolbox (`pdepe`)
- Parallel Computing Toolbox (recommended; the fit is computationally heavy)

## Usage

1. Download the data from the Zenodo deposit and place
   `SampleA_2_extracted_nmr_data.mat` (or `SampleB_...`) in the working folder.
   If the file is absent, the script runs on synthetic random data (smoke test
   only).
2. Edit the physical constants in the header of `RDNP_2components_github.m`
   (T1 values, 1H concentrations, diffusion coefficient `D0`, MAS rate,
   enhancement `epsilon0`, particle radius `totalpixelsmaterial`) to match your
   sample. The two fitted components are generic (`comp1`/`comp2`); for the
   published samples comp1 = PS and comp2 = SiO2-CTAB.
3. Run `RDNP_2components_github.m`. Runtime is hours to days depending on
   hardware; a high-performance computing (HPC) environment is recommended.
   The optimization stops at convergence or after `maxtimecluster` (24 h by
   default).
4. The fitted composition is saved in `SampleA_4_composition_fitted.mat`,
   which contains the fitted comp1 fraction per pixel of the particle, the full
   `fit_result` from `patternsearch`, and the final residue `sse1`.

## Differences from the Zenodo deposit

This repository provides a generalized two-component version of the
sample-specific scripts published on Zenodo:

- generic `comp1`/`comp2` notation instead of PS/CTAB, so any two solid
  components can be fitted;
- no hardcoded cluster settings: local `parpool` fallback, data-file fallback;
- correction of a concentration weighting in the comp2 signal
  (`fc_comp2` now uses `conc_comp2`; the deposited scripts used `conc_comp1`
  by mistake), so fits may differ slightly from the published results.

## Acknowledgment of prior work

This code is based on the Matlab codes published in:

- Rossini, A. J. *et al.* *J. Am. Chem. Soc.* **2012**, 134, 16899–16908. https://doi.org/10.1021/ja308135r
- Pinon, A. C. *et al.* *J. Phys. Chem. C* **2017**, 121, 15993–16005. https://doi.org/10.1021/acs.jpcc.7b04438
- Berruyer, P. *et al.* *J. Phys. Chem. Lett.* **2020**, 11, 8386–8391. https://doi.org/10.1021/acs.jpclett.0c02493
