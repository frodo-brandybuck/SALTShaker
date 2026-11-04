.. _gettingstarted:

***********************
Getting Started Quickly
***********************

To make it easy to run the latest SALTShaker training, the `trainsalt -g`
command downloads and unpacks the training data/config files into a local directory called
`saltshaker-latest-training`.  You can run the full training (takes a few hours) with the
following commands::

  conda activate saltshaker
  trainsalt -g
  cd saltshaker-latest-training
  trainsalt -c traingradient.conf

Training results will be in the `output/` directory.

Training from a parquet file
============================

Instead of a list of SNANA light-curve files, `snlists` can point to a single
parquet file that holds all of the training SNe::

  trainsalt -c traingradient.conf --snlists SALT3TRAIN_PanPlus.parquet

The parquet file has one row per SN.  Each row needs these metadata columns:

* `SNID`, `SURVEY`, `RA`, `DEC`, `MWEBV`, `SEARCH_PEAKMJD`, `FILTERS`
* `REDSHIFT_HELIO`, `REDSHIFT_CMB` and `VPEC`, each with a matching `_err` column

It also needs two nested columns:

* `lc`, the light curve, with the same columns as the `VARLIST` of a SNANA
  file (`MJD`, `FLT`, `FLUXCAL`, `FLUXCALERR`, ...).
* `spec`, all of the SN's spectra stacked together, with columns `LAMMIN`,
  `LAMMAX`, `FLAM`, `FLAMERR`, `SPECFLAG`, `MJD` and `SPECTRUM_ID`.  Rows
  with the same `SPECTRUM_ID` form one spectrum, so several spectra can share
  an MJD.  If `SPECTRUM_ID` is missing, rows with the same `MJD` are treated
  as one spectrum.  Use null for an SN without spectra.

Values must be the same as in the SNANA files.  In particular, `FLT` must be
the full filter name used in the kcor files (e.g. `CFA42-B/P`), and `SNID` must
match the IDs in the `tmaxlist` and `snparlist` files.  Peak MJDs are matched
by exact SNID, and SNe missing from `tmaxlist` are skipped.

Quick comparison with the SNANA files
-------------------------------------

To check that a parquet file gives the same training as the original files,
train on one survey with each input and compare the `output_*` directories.
First make a copy of `test_gradientspec.conf` named `quick_gradientspec.conf`
and lower `gradientmaxiter` to 50, so that each run takes a few minutes.
Then, with `SALT3TRAIN_PanPlus_CSPDR3.parquet` holding the CSPDR3 SNe::

  cd saltshaker-latest-training
  trainsalt -c traingradient.conf --maxsn None --phaserange -20 50 \
      --trainingconfig quick_gradientspec.conf --outputdir output_dat \
      --snlists SALT3TRAIN_PanPlus/SALT3TRAIN_PanPlus_CSPDR3/SALT3TRAIN_PanPlus_CSPDR3.LIST
  trainsalt -c traingradient.conf --maxsn None --phaserange -20 50 \
      --trainingconfig quick_gradientspec.conf --outputdir output_parquet \
      --snlists SALT3TRAIN_PanPlus_CSPDR3.parquet

Use `--maxsn None` rather than a `maxsn` cut.  With several `.LIST` files,
`maxsn` is split between the lists, but a parquet file counts as one list, so
the two runs would train on different SNe.  `--phaserange -20 50` keeps the
model phase range inside the default SALT2 starting templates, which end at
+50 days.

Both runs read the same data: 13 SNe, with 1498 photometric and 6105
spectroscopic observations after cuts.  A short run like this is not
converged, so repeated runs on the *same* input can end at different
chi-squared values.  Compare the parquet run against that run-to-run spread,
not against an exact match.
