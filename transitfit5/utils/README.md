# TransitFit5 Python Extensions (f2py)

This directory contains Fortran source files and f2py signature files (`.pyf`) to build Python extension modules for `transitfit5`:

- **`tfit5`** (`tfit5.pyf`): Python bindings for fast transit and occultation models, Keplerian orbit solver, transit duration calculation, TTV interpolation, and quicksort routines.
- **`fittransitmodel`** (`fittransitmodel3.pyf`): Python bindings for Levenberg-Marquardt photometric transit model optimization using MINPACK.

---

## Prerequisites

1. **Python 3** with **NumPy** (`numpy.f2py` or `f2py` CLI).
2. **Fortran Compiler**: `gfortran` (or Intel Fortran `ifx`/`ifort`).
3. **Build Backend**:
   - For modern Python (>= 3.12) and NumPy (>= 1.26), `f2py` uses the **Meson** backend. Ensure `meson` and `ninja` are installed:
     - On Fedora / RHEL: `sudo dnf install python3-meson-python python3-ninja`
     - Via pip: `pip install meson ninja`
4. **OpenMP**: `transitmodel.f` utilizes OpenMP multi-threading (`USE OMP_LIB`). OpenMP runtime support (e.g. `libgomp`) is required.

---

## Compilation Instructions

### Option 1: Compiling from `transitfit5/utils/`

From within this directory (`cd transitfit5/utils`):

1. **Build `tfit5`**:
   ```bash
   python3 -m numpy.f2py -c tfit5.pyf transitmodel.f keplerian.f ttcor.f occultquad.f mandelagol.f rqsort.f transitdur.f --dep openmp
   ```

2. **Build `fittransitmodel`**:
   ```bash
   python3 -m numpy.f2py -c fittransitmodel3.pyf precision.f90 fittermod.f90 fittransitmodel3.f90 getrhosig.f minpack.f transitmodel.f occultquad.f keplerian.f mandelagol.f ttcor.f --dep openmp
   ```

*(Note: `f2py` can be substituted for `python3 -m numpy.f2py` if available in your `$PATH`.)*

---

### Option 2: Compiling from Repository Root

From the root of the Kepler repository:

1. **Build `tfit5`**:
   ```bash
   python3 -m numpy.f2py -c ./transitfit5/utils/tfit5.pyf \
       ./transitfit5/utils/transitmodel.f \
       ./transitfit5/utils/keplerian.f \
       ./transitfit5/utils/ttcor.f \
       ./transitfit5/utils/occultquad.f \
       ./transitfit5/utils/mandelagol.f \
       ./transitfit5/utils/rqsort.f \
       ./transitfit5/utils/transitdur.f \
       --dep openmp
   ```

2. **Build `fittransitmodel`**:
   ```bash
   python3 -m numpy.f2py -c ./transitfit5/utils/fittransitmodel3.pyf \
       ./transitfit5/utils/precision.f90 \
       ./transitfit5/utils/fittermod.f90 \
       ./transitfit5/utils/fittransitmodel3.f90 \
       ./transitfit5/utils/getrhosig.f \
       ./transitfit5/utils/minpack.f \
       ./transitfit5/utils/transitmodel.f \
       ./transitfit5/utils/occultquad.f \
       ./transitfit5/utils/keplerian.f \
       ./transitfit5/utils/mandelagol.f \
       ./transitfit5/utils/ttcor.f \
       --dep openmp
   ```

---

## Compiler Flags & Notes

- **OpenMP Dependency**:
  - Meson backend (NumPy >= 1.26): Use `--dep openmp`.
  - Classic / distutils backend: Replace `--dep openmp` with `--f77flags='-fopenmp' --f90flags='-fopenmp' -lgomp`.
- **Intel Fortran Compiler (optional)**:
  - If compiling with Intel compilers instead of GCC/gfortran:
    ```bash
    f2py -c tfit5.pyf transitmodel.f keplerian.f ttcor.f occultquad.f mandelagol.f rqsort.f transitdur.f --f77flags='-qopenmp' --f90flags='-qopenmp' -liomp5
    ```

---

## Verification

Verify the compiled modules can be imported in Python:

```bash
python3 -c "import tfit5; print('tfit5 loaded successfully')"
python3 -c "import fittransitmodel; print('fittransitmodel loaded successfully')"
```
