#!/usr/bin/env python
from __future__ import annotations

import argparse

import h5py
import numpy as np
import pandas as pd
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


def sdhdf_access_data(f, dset, out_pth):
    """
    Demonstrate access of HDF attributes and astropy QTable data objects
    :param string f: Path to SDHDF file
    :return: None
    """
    try:
        # with the SDHDF file open read-only
        with h5py.File(f, "r") as h5:
            # list the HDF attributes

            # load the data for a particular sub-band into a numpy array e.g.:
            data = np.array(h5[dset])

            # print out data
            df = pd.DataFrame(data)
            for col, dtype in df.dtypes.items():
                if dtype == object:
                    # Only process byte object columns.
                    df[col] = df[col].apply(lambda x: x.decode("utf-8"))
            df.to_csv(out_pth, sep=" ", index=False, header=True)

    except Exception:
        pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--filename", help="Path to SDHDF file to read", required=True)
    ap.add_argument("--dataset", help="Path to dataset to read", required=True)
    ap.add_argument("--outfile", help="CSV file to write to", required=True)
    args = ap.parse_args()

    sdhdf_access_meta(args.filename)
    sdhdf_access_data(args.filename, args.dataset, args.outfile)
