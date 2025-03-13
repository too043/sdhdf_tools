#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read and write SDHDF tables"""

__author__ = ["Danny Price", "Alec Thomson"]

import pandas as pd
import h5py
from astropy.table import QTable

# Ignore astropy warnings
import warnings
warnings.filterwarnings('ignore', category=Warning, append=True)


class SDHDFAttribute:

    def __init__(self, sdhdf_dataset: h5py.Dataset, *args, **kwargs):
        """Read SDHDF attributes

        Args:
            sdhdf_dataset (h5py.Dataset): SDHDF table dataset
        """
        self.attrs = dict(sdhdf_dataset.attrs)
        self.table = pd.DataFrame((self.attrs.keys(), self.attrs.values()))

    def __repr__(self):
        return self.table.__repr__()

    def __str__(self):
        return self.table.__str__()

    def _repr_html_(self):
        return self.table._repr_html_()

    def __getitem__(self, key):
        return self.table.__getitem__(key)

    def __setitem__(self, key, value):
        return self.table.__setitem__(key, value)

    def __len__(self):
        return self.table.__len__()

    def __iter__(self):
        return self.table.__iter__()

    def __contains__(self, key):
        return self.table.__contains__(key)

    @staticmethod
    def _decode_df(df: pd.DataFrame) -> pd.DataFrame:
        """Decode a pandas dataframe to a string"""
        str_df = df.select_dtypes([object])
        str_df = str_df.stack().str.decode("utf-8").unstack()
        for col in str_df:
            df[col] = str_df[col]
        return df