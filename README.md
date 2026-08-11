# Kepler: Photometric Transit Modeling & Analysis Suite

Source code for NASA's Kepler Mission photometric transit modeling, light curve detrending, transit discovery, and Markov Chain Monte Carlo (MCMC) parameter estimation.

**Author**: Jason F. Rowe ([jasonfrowe@gmail.com](mailto:jasonfrowe@gmail.com))

---

## Prerequisites

* **Fortran Compiler**: GFortran (GCC 10+) or Intel Fortran (`ifort`/`ifx`) with OpenMP support.
* **C Compiler**: GCC or Clang (`gcc`).
* **Graphics**: X11 / XQuartz (for interactive PGPLOT plotting).
* **Libraries**:
  * **CFITSIO**: For reading MAST FITS light curve files.
  * **PGPLOT**: For interactive light curve, model, and posterior distribution plots.

### 1. Installing Dependencies on macOS (Homebrew)

Install CFITSIO and XQuartz via Homebrew:
```bash
brew install cfitsio
brew install --cask xquartz
```

#### Building PGPLOT from Source on macOS
PGPLOT is compiled locally (e.g. into `~/Software/pgplot/build`):
```bash
mkdir -p ~/Software/pgplot/build
cd ~/Software/pgplot/build
cp ../drivers.list .
```
1. Enable the desired drivers in `drivers.list` (e.g., `/NULL`, `/PNG`, `/TPNG`, `/PS`, `/VPS`, `/CPS`, `/VCPS`, `/XWINDOW`, `/XSERVE`).
2. Generate the Makefile with `makemake` using the Darwin GFortran configuration:
   ```bash
   ~/Software/pgplot/makemake ~/Software/pgplot darwin gfortran_gcc
   ```
3. Build the library:
   ```bash
   make libpgplot.a
   ```

#### PGPLOT Environment Setup
Add the following to your `~/.zshrc` or `~/.bashrc`:
```bash
export PGPLOT_DIR=$HOME/Software/pgplot/build/
export PGPLOT_DEV=/xserve   # Preferred interactive output device
```

### 2. Installing Dependencies on Linux

Install required packages via your distribution's package manager:
```bash
# Debian / Ubuntu
sudo apt install gfortran libcfitsio-dev pgplot5 libpng-dev

# Fedora / RHEL / CentOS
sudo dnf install gcc-gfortran cfitsio-devel pgplot-devel libpng-devel
```

---

## Quick Installation

### Option A: Standard Build (Autoconf) — *Recommended*
```bash
./autogen.sh
./configure
make -j4
```

### Option B: Direct Build (Make)
```bash
make -j4
```

> **Notes**:
> * `./autogen.sh` generates `./configure` from `configure.ac` using `autoreconf`.
> * Compiled binaries are automatically placed in the centralized `./bin/` directory.

---

## Overview of Modules & Executables

### 1. `TRANSITFIT5`
Photometric transit modeling, radial velocities, transit timing variations (TTVs), and Differential Evolution MCMC (deMCMC) posteriors.

| Executable | Description |
| :--- | :--- |
| **`transitfit51`** | Modern Fortran 90 Levenberg-Marquardt optimizer for multi-planet transit models. |
| **`transitfit5`** | Levenberg-Marquardt photometric transit model fitter. |
| **`transitmcmc5`** | Parallelized Differential Evolution MCMC (deMCMC) posterior sampler. |
| **`transitcut5`** | Trims light curve photometry beyond $\pm$ transit durations. |
| **`sigclip`** | Outlier rejection routine via iterative $\sigma$-clipping. |
| **`transitsn5`** | Estimates transit signal-to-noise ratios (S/N) based on transit models. |
| **`transittiming5`** | Fits individual center-of-transit times for TTV analysis. |
| **`transitdepth5`** | Computes transit depths from model parameters. |
| **`transitplot5`** | Interactive PGPLOT visualization of phase-folded light curves and transit models. |

---

### 2. `DATATEST`
MAST FITS light curve reading, data processing, and polynomial detrending routines.

| Executable | Description |
| :--- | :--- |
| **`detrend51`** | Parallelized out-of-transit running polynomial baseline detrending filter. |
| **`detrend5`** | Detrends light curves while preserving transit signals via input transit models. |
| **`kfitsread`** | Reads MAST Kepler FITS light curve files (CFITSIO). |
| **`sigclip`** | Outlier clipping utility. |

---

### 3. `TRANSITFIND`
Exoplanet transit discovery and search algorithms.

| Executable | Description |
| :--- | :--- |
| **`transitfind5`** | Optimized Box-fitting Least Squares (BLS) transit search algorithm. |
| **`transitfind2`** | Standard BLS transit search. |
| **`pdmsearch`** | Phase Dispersion Minimization (PDM) periodic signal search. |

---

## Citation & References

If you make use of this codebase in your research, please cite:

* **Rowe et al. (2014)**, *Validation of Kepler's Multiple Planet Candidates. III. Light Curve Analysis & Performance of Medea*, ApJ, 784, 45. [[DOI:10.1088/0004-637X/784/1/45](https://doi.org/10.1088/0004-637X/784/1/45)]
* **Rowe et al. (2015)**, *Planetary Candidates Observed by Kepler. V. Planet Sample from Q1–Q16 (47 Months)*, ApJS, 217, 16. [[DOI:10.1088/0067-0049/217/1/16](https://doi.org/10.1088/0067-0049/217/1/16)]
* **Jason Rowe (2016)**, *Kepler: Kepler Transit Model Codebase Release* [Data set]. Zenodo. [[DOI:10.5281/zenodo.60297](http://doi.org/10.5281/zenodo.60297)]
