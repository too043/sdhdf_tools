#!/usr/bin/env python
from __future__ import annotations

import argparse
import glob
import os

import h5py
from astropy.table import QTable
from tabulate import tabulate

__version__ = "2.2"
__author__ = "Lawrence Toomey"


def get_metadata(pth):
    """
    Retrieve metadata from SDHDF dataset

    :param string pth: Path to SDHDF file to read
    :return: list row_d: Metadata as a list
    """
    row_d = []

    try:
        with h5py.File(pth, "r") as h5:
            ph = QTable.read(h5, path="/metadata/primary_header")
            bp = QTable.read(h5, path="/metadata/beam_params")
            f = os.path.basename(pth)
            nbeams = len(bp)
            pid = ph["PID"][0]
            utc_start = ph["UTC_START"][0]
            obs_type = ph["OBS_TYPE"][0]
            hdr_defn_ver = ph["HDR_DEFN_VERSION"][0]
            sched_blk_id = ph["SCHED_BLOCK_ID"][0]

            for beam in range(nbeams):
                beam_id = bp["LABEL"][beam]
                source = bp["SOURCE"][beam]
                nbands = bp["N_BANDS"]
                row = (
                    f,
                    hdr_defn_ver,
                    sched_blk_id,
                    pid,
                    beam_id,
                    source,
                    obs_type,
                    utc_start,
                    nbands,
                )
                row_d.append(row)

    except Exception as e:
        print("ERROR: failed to read file %s" % f, e)

    return row_d


def quick_list(pth):
    """
    Format metadata from SDHDF file
    :param string pth: Path to SDHDF file to read
    :return: None
    """
    row_data = []
    hdr = (
        "File",
        "SDHDF Version",
        "Sched Block ID",
        "Project ID",
        "Beam",
        "Source",
        "Obs Type",
        "UTC start",
        "No. bands",
    )

    if os.path.isdir(pth):
        for fpth in glob.glob(pth + "/*.hdf"):
            row = get_metadata(fpth)
            row_data += row

    else:
        row_data = get_metadata(pth)

    print(tabulate(row_data, headers=hdr))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", help="Path to directory containing SDHDF files")
    ap.add_argument("--filename", help="Path to SDHDF file")
    args = ap.parse_args()

    if args.dir:
        quick_list(args.dir)
    else:
        quick_list(args.filename)
