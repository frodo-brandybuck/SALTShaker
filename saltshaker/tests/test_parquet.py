"""Tests for reading SALT3 training data from a nested parquet file.

A small parquet file is built from the CSPDR3 test data, and each SN read
from it is compared with the same SN read from its SNANA .dat file.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from saltshaker.util import readutils, snana

DATADIR = 'testdata/data/SALT3TRAIN_K21_CSPDR3'
SNLIST = f'{DATADIR}/SALT3TRAIN_K21_CSPDR3.LIST'
TMAXLIST = 'testdata/test_pkmjd.LIST'
SNPARLIST = 'testdata/test_pars.LIST'
SPECRECALLIST = 'testdata/test_spcrl.LIST'

LC_COLS = ['MJD', 'FLT', 'FIELD', 'FLUXCAL', 'FLUXCALERR', 'MAG', 'MAGERR']
SPEC_COLS = ['LAMMIN', 'LAMMAX', 'FLAM', 'FLAMERR', 'SPECFLAG']


def _value(s):
    """'340.542  deg' -> 340.542, '0.047 # MW E(B-V)' -> 0.047"""
    return float(str(s).split('#')[0].split()[0])


def _value_err(s):
    """'0.031011 +- 0.000021' -> (0.031011, 0.000021)"""
    val, err = str(s).split('+-')
    return float(val), float(err)


def sn_to_parquet_row(sn):
    """Convert a SuperNova read from a .dat file to one parquet row."""
    row = {
        'SURVEY': sn.SURVEY, 'SNID': str(sn.SNID),
        'RA': _value(sn.RA), 'DEC': _value(sn.DEC), 'MWEBV': _value(sn.MWEBV),
        'SEARCH_PEAKMJD': float(sn.SEARCH_PEAKMJD), 'FILTERS': sn.FILTERS,
    }
    for key in ['REDSHIFT_HELIO', 'REDSHIFT_CMB', 'VPEC']:
        row[key], row[key + '_err'] = _value_err(getattr(sn, key))
    row['lc'] = {c: list(getattr(sn, c)) for c in LC_COLS}

    spectra = list(getattr(sn, 'SPECTRA', {}).values())
    if spectra:
        spec = {c: [] for c in SPEC_COLS + ['MJD', 'SPECTRUM_ID']}
        for specid, s in enumerate(spectra, start=1):
            for c in SPEC_COLS:
                spec[c] += list(np.asarray(s[c], dtype=float))
            spec['MJD'] += [float(s['SPECTRUM_MJD'])] * len(s['LAMMIN'])
            spec['SPECTRUM_ID'] += [specid] * len(s['LAMMIN'])
        row['spec'] = spec
    else:
        row['spec'] = None
    return row


def spectra_in_order(sn):
    return list(getattr(sn, 'SPECTRA', {}).values())


def assert_same_spectra(testcase, p, d):
    """Check that p and d hold the same spectra, one by one, in file order."""
    pspec, dspec = spectra_in_order(p), spectra_in_order(d)
    testcase.assertEqual(len(pspec), len(dspec))
    for ps, ds in zip(pspec, dspec):
        testcase.assertAlmostEqual(float(ps['SPECTRUM_MJD']),
                                   float(ds['SPECTRUM_MJD']))
        for c in ['LAMMIN', 'FLAM']:
            np.testing.assert_allclose(np.asarray(ps[c], float),
                                       np.asarray(ds[c], float))


class test_parquet(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        datfiles = np.atleast_1d(np.genfromtxt(SNLIST, dtype='str'))
        cls.datfiles = {}
        rows = []
        for f in datfiles:
            sn = snana.SuperNova(os.path.join(DATADIR, f))
            cls.datfiles[str(sn.SNID)] = os.path.join(DATADIR, f)
            rows.append(sn_to_parquet_row(sn))

        cls.tmpdir = tempfile.mkdtemp()
        cls.parquetfile = os.path.join(cls.tmpdir, 'test.parquet')
        pd.DataFrame(rows).to_parquet(cls.parquetfile)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)

    def test_snana_parquet_matches_dat(self):
        nwithspec = 0
        for snid, datfile in self.datfiles.items():
            with self.subTest(snid=snid):
                p = snana.SuperNova(snid=snid, parquetfile=self.parquetfile)
                d = snana.SuperNova(datfile)

                np.testing.assert_array_equal(np.asarray(p.MJD), d.MJD)
                np.testing.assert_array_equal(np.asarray(p.FLT), d.FLT)
                np.testing.assert_allclose(np.asarray(p.FLUXCAL, float), d.FLUXCAL)
                np.testing.assert_allclose(np.asarray(p.FLUXCALERR, float), d.FLUXCALERR)
                self.assertAlmostEqual(float(p.SEARCH_PEAKMJD), float(d.SEARCH_PEAKMJD))
                self.assertAlmostEqual(_value_err(p.REDSHIFT_HELIO)[0],
                                       _value_err(d.REDSHIFT_HELIO)[0])

                assert_same_spectra(self, p, d)
                nwithspec += len(spectra_in_order(d)) > 0
        self.assertGreater(nwithspec, 0)

    def test_spectra_on_same_mjd_stay_separate(self):
        # Some SNe (e.g. 1991T) have several spectra on one MJD. Give one SN
        # an extra spectrum on the MJD of its first one and check that the
        # parquet reader keeps them apart.
        snid, datfile = next((k, f) for k, f in self.datfiles.items()
                             if spectra_in_order(snana.SuperNova(f)))
        d = snana.SuperNova(datfile)
        first = spectra_in_order(d)[0]
        extra = dict(first)
        extra['FLAM'] = np.asarray(first['FLAM'], float) * 2
        d.SPECTRA = dict(enumerate([first, extra] + spectra_in_order(d)[1:]))

        row = sn_to_parquet_row(d)
        parquetfile = os.path.join(self.tmpdir, 'samemjd.parquet')
        pd.DataFrame([row]).to_parquet(parquetfile)
        p = snana.SuperNova(snid=snid, parquetfile=parquetfile)
        self.assertEqual(p.NSPECTRA, len(d.SPECTRA))
        assert_same_spectra(self, p, d)

        # Without SPECTRUM_ID the reader can only group by MJD, so the two
        # spectra on the same MJD come back as one.
        del row['spec']['SPECTRUM_ID']
        parquetfile = os.path.join(self.tmpdir, 'samemjd_noid.parquet')
        pd.DataFrame([row]).to_parquet(parquetfile)
        p = snana.SuperNova(snid=snid, parquetfile=parquetfile)
        self.assertEqual(p.NSPECTRA, len(d.SPECTRA) - 1)

    def test_rdAllData_parquet_matches_list(self):
        kwargs = dict(dospec=True, peakmjdlist=TMAXLIST, snparlist=SNPARLIST,
                      specrecallist=SPECRECALLIST)
        fromlist = readutils.rdAllData(SNLIST, False, **kwargs)
        fromparquet = readutils.rdAllData(self.parquetfile, False, **kwargs)

        self.assertGreater(len(fromlist), 0)
        self.assertEqual(sorted(fromlist), sorted(fromparquet))
        for snid in fromlist:
            with self.subTest(snid=snid):
                a, b = fromlist[snid], fromparquet[snid]
                self.assertEqual(a.num_photobs, b.num_photobs)
                self.assertEqual(a.num_spec, b.num_spec)
                self.assertEqual(a.num_specobs, b.num_specobs)
                self.assertAlmostEqual(a.zHelio, b.zHelio)
                self.assertAlmostEqual(a.tpk_guess, b.tpk_guess)
        self.assertTrue(any(sn.num_spec > 0 for sn in fromparquet.values()))


if __name__ == "__main__":
    unittest.main()
