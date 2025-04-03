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
            print(
                "\nPrimary header attribute keys: \n%s"
                % h5["metadata/primary_header"].attrs.keys()
            )
            print(
                "\nBand params attribute keys: \n%s"
                % h5["beam_0/metadata/band_params"].attrs.keys()
            )

            print("\nPID attributes: \n%s" % h5["metadata/primary_header"].attrs["PID"])
            print("\nPID data: \n%s" % h5["metadata/primary_header"]["PID"])
            print("\nPID data type: \n%s" % h5["metadata/primary_header"]["PID"].dtype)

            print(
                "\nN_CHANS attributes: \n%s"
                % h5["beam_0/metadata/band_params"].attrs["N_CHANS"]
            )
            print("\nN_CHANS data: \n%s" % h5["beam_0/metadata/band_params"]["N_CHANS"])

            # ...or from an astropy.QTable object...e.g.:
            tb = QTable.read(h5, path="metadata/primary_header")
            op = QTable.read(h5, path="beam_0/metadata/band_params")

            print("\nPID attributes from astropy.QTable: \n%s" % tb.meta["PID"])
            print("\nN_CHANS attributes from astropy.QTable: \n%s" % op.meta["N_CHANS"])

            print("\nPID data from astropy.QTable: \n%s" % tb["PID"])
            print("\nN_CHANS data from astropy.QTable: \n%s" % op["N_CHANS"])

            # the QTable object dtypes can be accessed with e.g.:
            print("\nPrimary header data types: \n%s" % tb.info)
            print("\nBand params data types: \n%s" % op.info)

    except Exception as e:
        print("ERROR: failed to read file %s" % f, e)


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
            print(
                "\nData attribute keys: \n%s"
                % h5["beam_0/band_SB0/astronomy_data/data"].attrs.keys()
            )
            print(
                "\nFrequency attribute keys: \n%s"
                % h5["beam_0/band_SB0/astronomy_data/frequency"].attrs.keys()
            )

            # describe the dimensions of the data
            print(
                "\nData dimension labels: \n%s"
                % h5["beam_0/band_SB0/astronomy_data/data"].attrs["DIMENSION_LABELS"]
            )

            # describe the frequency axis unit
            print(
                "\nFrequency unit: \n%s"
                % h5["beam_0/band_SB0/astronomy_data/frequency"].attrs["UNIT"]
            )

            # load the data for a particular sub-band into a numpy array e.g.:
            data = np.array(h5["beam_0/band_SB0/astronomy_data/data"])
            print("\nData array shape: \n%s" % str(data.shape))

            # load the frequency axis for sub-band 5 into an Numpy array e.g.:
            freq_axis = np.array(h5["beam_0/band_SB0/astronomy_data/frequency"])
            print("\nFrequency array shape: \n%s" % str(freq_axis.shape))
    except Exception as e:
        print("ERROR: failed to read file %s" % f, e)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--filename", help="Path to SDHDF file to read", required=True)
    args = ap.parse_args()

    sdhdf_access_meta(args.filename)
    sdhdf_access_data(args.filename)
