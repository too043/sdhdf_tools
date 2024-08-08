#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SDHDF history utilities"""

import datetime
import inspect
import socket
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pkg_resources
from astropy.stats import mad_std, sigma_clip
from astropy.table import Table
from IPython import embed
from tqdm.auto import tqdm
from xarray import DataArray, Dataset, Variable
from .tables import SDHDFTable


def generate_history_row() -> pd.DataFrame:
    """Generate a history row.
    Returns:
        pd.DataFrame: History row.
    """

    # Get the calling function from inspect
    process_name = inspect.stack()[1][3]
    # Get the calling function's docstring
    process_description = inspect.stack()[1][0].f_locals["self"].__doc__
    # Get the calling function's arguments
    process_arguments = str(inspect.stack()[1][0].f_locals["self"].__dict__)

    history_row = pd.DataFrame(
        {
            "DATE": datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S"),
            "PROC": process_name,
            "PROC_DESCR": process_description,
            "PROC_ARGS": process_arguments,
            "PROC_HOST": socket.getfqdn(),
        },
        index=[0],
    )
    return history_row
