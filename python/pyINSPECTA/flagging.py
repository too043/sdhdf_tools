#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SDHDF flagging utilities"""

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import h5py
import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pkg_resources
from astropy.stats import mad_std, sigma_clip
from astropy.table import Table
from tqdm.auto import tqdm
from xarray import DataArray, Dataset, Variable


class AutoFlagError(Exception):
    def __init__(self, msg):
        self.msg = msg


def get_persistent_rfi(telescope: str = "Parkes") -> pd.DataFrame:
    """Read persistent RFI file
    Returns:
        pd.Dataframe: Persistent RFI data.
    """
    rfi_file = pkg_resources.resource_filename(
        "pyINSPECTA", f"{telescope.lower()}_rfi.csv"
    )
    if not Path(rfi_file).exists():
        raise NotImplementedError(
            f"Persistent RFI file for {telescope} not found at '{Path(rfi_file).absolute()}'."
        )
    rfi_df = pd.read_csv(
        rfi_file,
        sep=",",
        # skip_blank_lines=True,
        comment="#",
        names=[
            "type",
            "observatory label",
            "receiver label",
            "freq0 MHz",
            "freq1 MHz",
            "MJD0",
            "MJD1",
            "text string for label",
        ],
    )
    return rfi_df


def box_filter(spectrum: np.ndarray | DataArray, sigma=3, n_windows=100):
    """
    Filter a spectrum using a box filter.
    """
    # Divide spectrum into windows
    window_size = len(spectrum) // n_windows
    dat_filt = np.zeros_like(spectrum).astype(bool)
    # Iterate through windows
    for i in range(n_windows):
        _dat = spectrum[i * window_size : window_size + i * window_size]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Use sigma clipping to remove outliers
            _dat_filt = sigma_clip(
                _dat, sigma=sigma, maxiters=None, stdfunc=mad_std, masked=True
            )
        dat_filt[i * window_size : window_size + i * window_size] = _dat_filt.mask
    if isinstance(spectrum, DataArray):
        return DataArray(dat_filt, dims=spectrum.dims, coords=spectrum.coords)
    
    return dat_filt
    

