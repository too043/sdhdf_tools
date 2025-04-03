#!/usr/bin/env python
from __future__ import annotations

import argparse

import h5py
import numpy as np
from astropy.table import QTable

__version__ = "2.2"
__author__ = "Lawrence Toomey"


def sdhdf_access_meta(f):
    """
    Demonstrate access of HDF attributes and astropy QTable metadata objects
    :param string f: Path to SDHDF file
    :return: None
    """
    try:
        # with the SDHDF file open read-only
        with h5py.File(f, "r") as h5:
            # the metadata can be accessed directly from the HDF attributes...e.g.:



            # ...or from an astropy.QTable object...e.g.:
            QTable.read(h5, path="metadata/primary_header")
            QTable.read(h5, path="beam_0/metadata/band_params")



            # the QTable object dtypes can be accessed with e.g.:

    except Exception:
        pass


def sdhdf_access_data(f):
    """
    Demonstrate access of HDF attributes and astropy QTable data objects
    :param string f: Path to SDHDF file
    :return: None
    """
    try:
        # with the SDHDF file open read-only
        with h5py.File(f, "r") as h5:
            # list the HDF attributes

            # describe the dimensions of the data

            # describe the frequency axis unit

            # load the data for a particular sub-band into a numpy array e.g.:
            np.array(h5["beam_0/band_SB0/astronomy_data/data"])

            # load the frequency axis for sub-band 5 into an numpy array e.g.:
            np.array(h5["beam_0/band_SB0/astronomy_data/frequency"])
    except Exception:
        pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--filename", help="Path to SDHDF file to read", required=True)
    args = ap.parse_args()

    sdhdf_access_meta(args.filename)
    sdhdf_access_data(args.filename)
