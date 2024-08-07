#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core SDHDF module
"""
from __future__ import annotations

import json
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pkg_resources
import xarray as xr
from astropy.table import Table
from dask.diagnostics import ProgressBar
from dask.distributed import Client, get_client, get_task_stream, progress
from tqdm.auto import tqdm
from xarray import DataArray, Dataset, Variable

from pyINSPECTA import flagging, history
from pyINSPECTA.logger import logger
from pyINSPECTA.tables import SDHDFTable
from pyINSPECTA.attributes import SDHDFAttribute


def _get_sdhdf_version(filename: Path) -> tuple[float, Path]:
    """Get the SDHDF version of a file and return the path
       to the definition template

    Args:
        filename (Path): Path to the SDHDF file

    Returns:
        str: SDHDF version
        str: Path to definition json template
    """
    with h5py.File(filename, "r") as f:
        hdr_keys = ["HDR_DEFN_VERSION", "HEADER_DEFINITION_VERSION"]
        primary_header = SDHDFTable(f["metadata/primary_header"])
        version = None
        for hdr_key in hdr_keys:
            if hdr_key in primary_header:
                version = float(primary_header[hdr_key])

        if version is None:
            raise ValueError(f"SDHDF version not found in file '{filename}'")
        elif version <= 2.0:
            version = 2.0
        elif (version > 2.0) and (version <= 2.1):
            version = 2.1
        elif (version > 2.2) and (version <= 2.9):
            version = 2.9

        try:
            definition_file = Path(pkg_resources.resource_filename(
                "pyINSPECTA", f"definitions/sdhdf_def_v{version}.json"
            ))
        except ValueError as e:
            raise ValueError(f"SDHDF definition template {definition_file} not found.") from e

    return version, definition_file


@dataclass
class MetaData:
    """An SDHDF metadata object

    Args:
        filename (Path): Path to the SDHDF file

    Attributes:
        beam_params (SDHDFTable): The beam parameters
        history (SDHDFTable): File history
        primary_header (SDHDFTable): Primary header
        backend_config (SDHDFTable): Backend configuration
        cal_backend_config (SDHDFTable): Calibration backend configuration
        software (SDHDFTable): Software versions used in creation of the file
        schedule (SDHDFTable): Observation schedule metadata (if available)

    Methods:
        print_obs_metadata: Quickly list the observation metadata
        print_obs_config: Quickly list the observation configuration
        write: Write metadata to file [NOT YET IMPLEMENTED]

    """

    filename: Path

    def __post_init__(self):
        version, definition_file_path = _get_sdhdf_version(self.filename)
        logger.info(f"SDHDF version: {version}")
        logger.debug(f"Loading SDHDF definition template: {definition_file_path}")

        # load the definition
        with open(definition_file_path, "r") as definition_file:
            self.definition = json.load(definition_file)

        # load the metadata
        with h5py.File(self.filename, "r") as file:
            logger.debug(json.dumps(self.definition, indent=4))
            for key in file.keys():
                logger.info(f"Found SDHDF group '{key}'...")
                group_key = "beam" if "beam" in key else key
                if group_key in self.definition:
                    logger.info(f"Loading metadata for '{group_key}'...")
                    self.load_metadata(file, self.definition[group_key], key)
                else:
                    logger.warning(
                        f"Key '{group_key}' not found in definition file '{definition_file_path}'. Ignoring..."
                    )

    def load_metadata(self, f, definition, path=""):
        for k, v in definition.items():
            current_path = f"{path}/{k}" if path else k
            if k == "attributes":
                attr = SDHDFAttribute(f[path])
                setattr(self, k, attr)
                logger.debug(getattr(self, k))
            elif k == "band":
                bands = list(filter(lambda element: "band" in element, f[path].keys()))
                for band in bands:
                    logger.info(f"Loading metadata for beam '{path}' band '{band}'...")
                    self.load_metadata(f, v, f"{path}/{band}")
            elif isinstance(v, dict):
                self.load_metadata(f, v, current_path)
            elif v in f:
                tab = SDHDFTable(f[v])
                logger.debug(tab)
                setattr(self, k, tab)
            else:
                logger.warning(f"No object '{current_path}' found in file!")

    def print_obs_metadata(self, format: str = "grid") -> None:
        """Print observation metadata to the terminal"""
        for key in self.definition["metadata"].keys():
            if key in self.__dict__:
                df = self.__dict__[key]
                logger.info(f"{key}:")
                print(df.table.to_markdown(tablefmt=format, headers=[]))
            else:
                logger.warning(
                    f"No metadata found for key '{key}'. Ignoring..."
                )

    def print_obs_config(self, format: str = "grid") -> None:
        """Print the observation configuration to the terminal"""
        for key in self.definition["config"].keys():
            if key in self.__dict__:
                df = self.__dict__[key]
                logger.info(f"{key}:")
                print(df.table.to_markdown(tablefmt=format, headers=[]))
            else:
                logger.warning(
                    f"No metadata found for key '{key}'. Ignoring..."
                )

    def write(self, filename: str | Path, overwrite:bool=False) -> pd.DataFrame:
        """Write the metadata to a file

        Args:
            filename (str | Path): Path to the output file
        """
        if isinstance(filename, str):
            filename = Path(filename)

        if filename.exists() and not overwrite:
            raise FileExistsError(f"File '{filename}' already exists")
        for name in ("metadata", "config"):
            for key, val in tqdm(self.definition[name].items(), desc=f"Writing {name}"):
                df = self.__dict__[key]
                df.to_hdf(filename, key=f"{val}", mode="a", data_columns=True)

        return history.generate_history_row()


@dataclass
class SubBand:
    """An SDHDF sub-band data object

    Args:
        label (str): Sub-band label
        filename (Path): Path to the SDHDF file
        definition (dict): SDHDF definition
        beam_label (str): Beam label
        in_memory (bool, optional): Load the data into memory. Defaults to False.
        client (Client, optional): Dask client. Defaults to None.

    Attributes:
        data (DataArray): The sub-band data as an xarray DataArray
        flag (DataArray): The sub-band flag as an xarray DataArray
        meta (DataFrame): The sub-band metadata as a pandas DataFrame

    Methods:
        plot_waterfall: Plot the sub-band data as a waterfall plot
        plot_spectrum: Plot a single spectrum from the sub-band data

    """

    label: str
    filename: Path
    definition: dict
    beam_label: str
    in_memory: bool = False
    client: Client | None = None

    def __post_init__(self):
        # Get the astronomy data
        self.astronomy_dataset = self._get_data()
        # Now get the calibrator data
        self.calibrator_dataset = self._get_cal()
        # TODO: Get the calibrator data

    def _get_cal(self):
        return

    @staticmethod
    def _get_data_dimensions(
        dim_labels: str,
        data_shape: tuple[int],
        meta: SDHDFTable,
    ) -> list[str]:
        if dim_labels != "NOT SET":
            return dim_labels.split(",")

        logger.warning("No dimension labels found in file! Using default labels.")
        usual_dims = ["time", "polarization", "channel", "beam"]
        dims = []
        for i, shape in enumerate(data_shape):
            if shape == len(meta):
                dims.append("time")
            else:
                # dims.append(f"dim_{i}")
                dims.append(usual_dims[i])
        return dims

    @staticmethod
    def _get_freq_dimensions(
        dims: list[str],
        freq_dim_labels: str,
        data_shape: tuple[int],
        freq_shape: tuple[int],
    ) -> list[str]:

        if freq_dim_labels != "NOT SET":
            return freq_dim_labels.split(",")

        logger.warning(
            "No frequency dimension labels found in file! Attempting to match dimensions to data shape..."
        )
        # Attempt to match the dimensions to the data shape
        freq_dims = []
        for i_shape in freq_shape:
            for i in range(len(data_shape)):
                if i_shape == data_shape[i]:
                    freq_dims.append(dims[i])
                    break

        if len(freq_dims) == len(freq_shape):
            return freq_dims

        logger.warning(
            "Could not automatically match frequency dimensions to data dimensions! Using default labels."
        )
        return [f"freq_{i}" for i in range(len(freq_shape))]

    def _get_data(self):
        """Get the astronomy sub-band data"""
        astro_def = self.definition["beam"]["band"]["astronomy"]
        meta_def = self.definition["beam"]["band"]["metadata"]
        sb_path = f"{self.beam_label}/{self.label}"

        with h5py.File(self.filename, "r") as h5:
            data_path = f"{sb_path}/{astro_def['data']}"
            freq_path = f"{sb_path}/{astro_def['frequency']}"
            meta_path = f"{sb_path}/{meta_def['obs_params']}"

            data = h5[data_path]
            freqs = h5[freq_path]
            meta = SDHDFTable(h5[meta_path])
            self.metadata = meta

            # Get the flags (if they exist)
            flag_path = f"{sb_path}/{astro_def['flags']}"
            if "flags" in astro_def.keys() and flag_path is True:
                flags = h5[flag_path][:]
                # Ensure flag has same shape as data
                flag_reshape = flags[:].copy()
                for i, s in enumerate(data.shape):
                    if i > len(flag_reshape.shape) - 1:
                        flag_reshape = np.expand_dims(flag_reshape, axis=-1)
                    else:
                        if flag_reshape.shape[i] == s:
                            continue
                        else:
                            flag_reshape = np.expand_dims(flag_reshape, axis=i)
                flags = flag_reshape
            else:
                logger.warning(f"No flags found for sub-band '{self.label}' in file '{self.filename}'!")
                logger.warning(f"Band '{self.label}' flags will be set to all zeros.")
                flags = np.zeros_like(data)

            # Load into memory if requested
            if self.in_memory:
                logger.info(f"Loading {self.label} into memory...")
                data = np.array(data)
                freqs = np.array(freqs)
                flags = np.array(flags)

            # Process into xarray
            coords = {col: ("time", meta[col].values) for col in meta.table.columns}

            dims = self._get_data_dimensions(
                dim_labels=h5[data_path].attrs["DIMENSION_LABELS"].decode(),
                data_shape=data.shape,
                meta=meta,
            )

            freq_dims = self._get_freq_dimensions(
                dims=dims,
                freq_dim_labels=h5[freq_path].attrs["DIMENSION_LABELS"].decode(),
                data_shape=data.shape,
                freq_shape=freqs.shape,
            )
            coords["frequency"] = Variable(
                dims=freq_dims,
                data=freqs,
                attrs={"units": h5[freq_path].attrs["UNIT"].decode()},
            )

            attrs = dict(h5[data_path].attrs)
            for key, val in attrs.items():
                if isinstance(val, bytes):
                    attrs[key] = val.decode()

            data_xr = DataArray(
                data,
                dims=dims,
                coords=coords,
                name=f"{self.label}_data",
                attrs=attrs,
            )
            # Check if data has beam dimension
            if "beam" in data_xr.dims:
                data_xr = data_xr.isel(beam=0)
            data_xr.attrs["units"] = data_xr.UNIT
            self.attrs = attrs

            flag_xr = DataArray(
                flags,
                dims=dims,
                coords=coords,
                name=f"{self.label}_flag",
            )
            # Same as above
            if "beam" in flag_xr.dims:
                flag_xr = flag_xr.isel(beam=0)

            astronomy_dataset = Dataset(
                {
                    "data": data_xr,
                    "flag": flag_xr,
                    "metadata": xr.DataArray(meta.table, dims=["time", "meta"])
                }
            )

            return astronomy_dataset

    def plot_waterfall(
        self,
        polarization: int = 0,
        # bin: int = 0,
        flag: bool = False,
        **plot_kwargs,
    ):
        """Waterfall plot of the data

        Args:
            polarization (int, optional): Polarization to select. Defaults to 0.
            bin (int, optional): Bin to select. Defaults to 0.
            flag (bool, optional): Blank flagged data. Defaults to False.
        """
        sub_data = self.astronomy_dataset.isel(polarization=polarization,)
        if flag:
            sub_data = sub_data.where(sub_data.flag == 0)
        sub_data.data.plot(x="frequency", **plot_kwargs)
        ax = plt.gca()
        return ax

    def plot_spectrum(
        self,
        time,
        polarization: int = 0,
        bin: int = 0,
        flag: bool = False,
        **plot_kwargs,
    ):
        sub_data = self.astronomy_dataset.isel(
            time=time, polarization=polarization,
        )
        if flag:
            sub_data = sub_data.where(sub_data.flag == 0)
        sub_data.data.plot(**plot_kwargs)
        ax = plt.gca()
        return ax

    def autoflag(self, sigma=3, n_windows=100):
        """Automatic flagging using rolling sigma clipping"""
        data_xr_flg = self.astronomy_dataset.data.where(
            ~self.astronomy_dataset.flag.astype(bool)
        )
        # Set chunks for parallel processing
        chunks = {d: 1 for d in data_xr_flg.dims}
        chunks["channel"] = len(self.astronomy_dataset.data.channel)
        data_xr_flg = data_xr_flg.chunk(chunks)
        mask = xr.apply_ufunc(
            flagging.box_filter,
            data_xr_flg,
            input_core_dims=[["channel"]],
            output_core_dims=[["channel"]],
            kwargs={"sigma": sigma, "n_windows": n_windows},
            dask="parallelized",
            vectorize=True,
            output_dtypes=(bool),
        )
        self.astronomy_dataset["flag"] = mask.astype(int).compute()
        hist = history.generate_history_row()
        return hist

    def decimate(
        self, bins: float | int, axis: str = "frequency", use_median: bool = False
    ) -> pd.DataFrame:
        """Average the data along the an axis

        Args:
            bins (float | int): If int, the number of channels to bin in an average.
                If float, the desired width of a channel after averaging.
            axis (str, optional): The axis to decimate along. Defaults to "frequency".
            use_median (bool, optional): Use the median instead of the mean. Defaults to False.

        Returns:
            pd.DataFrame: The history row

        Raises:
            NotImplementedError: Decimation along the time axis is not yet implemented

        """

        if axis == "time":
            # TODO: Figure out how to decimate along the time axis - includeing the metadata / coords
            raise NotImplementedError(
                "Decimation along the time axis is not yet implemented"
            )
        dataset = self.astronomy_dataset
        if isinstance(bins, float):
            # Convert to integer number of bins
            try:
                unit = dataset.data[axis].units
            except AttributeError:
                unit = "units"
            logger.info(f"Asked for a bin width of {bins} {unit}")
            logger.info(
                f"Dimension {axis} has range {dataset[axis].min()} to {dataset[axis].max()}: {dataset[axis].max() - dataset[axis].min()} {unit}"
            )
            bins = int((dataset[axis].max() - dataset[axis].min()) / bins)

        logger.info(f"Using {bins} channels per bin")

        # Apply CASA-style decimation
        flagged = dataset.where(dataset.flag == 0)
        unflagged = dataset

        if use_median:
            unflagged_dec = (
                unflagged.coarsen(**{axis: bins}, boundary="trim")
                .construct(**{axis: ("decimated", "original")})
                .median(dim="original", skipna=True)
                .rename({"decimated": axis})
            )

            flagged_dec = (
                flagged.coarsen(**{axis: bins}, boundary="trim")
                .construct(**{axis: ("decimated", "original")})
                .median(dim="original", skipna=True)
                .rename({"decimated": axis})
            )
            axis_dec = unflagged[axis].coarsen(**{axis: bins}, boundary="trim").median()

        else:
            unflagged_dec = (
                unflagged.coarsen(**{axis: bins}, boundary="trim")
                .construct(**{axis: ("decimated", "original")})
                .mean(dim="original", skipna=True)
                .rename({"decimated": axis})
            )

            flagged_dec = (
                flagged.coarsen(**{axis: bins}, boundary="trim")
                .construct(**{axis: ("decimated", "original")})
                .mean(dim="original", skipna=True)
                .rename({"decimated": axis})
            )
            axis_dec = unflagged[axis].coarsen(**{axis: bins}, boundary="trim").mean()

        unflagged_dec[axis] = axis_dec
        flagged_dec[axis] = axis_dec
        new_flag = flagged_dec.flag.fillna(1)
        new_data = flagged_dec.data
        new_data = new_data.fillna(unflagged_dec.data)
        dataset_dec = xr.Dataset(
            {
                "data": new_data,
                "flag": new_flag,
                "metadata": dataset.metadata,
            },
            attrs=dataset.attrs,
        )
        self.astronomy_dataset = dataset_dec
        hist = history.generate_history_row()
        return hist

    def _write_astronomy_dataset(self, filename: Path) -> pd.DataFrame:
        astro_def = self.definition["subband"]["astronomy"]
        sb_path = f"{self.beam_label}/{self.label}"
        with h5py.File(filename, "w") as f:
            f[f"{sb_path}/{astro_def['data']}"] = self.astronomy_dataset.data.values
            f[f"{sb_path}/{astro_def['frequency']}"] = self.astronomy_dataset.frequency.values
            if "flags" in self.definition["subband"]["astronomy"]:
                f[f"{sb_path}/{astro_def['flags']}"] = self.astronomy_dataset.flag.values
            else:
                logger.warning("No flags in definition")
                logger.info("Saving flags to /astronomy_data/flags")
                f[f"{sb_path}/astronomy_data/flags"] = self.astronomy_dataset.flag.values

        self.astronomy_dataset.metadata.to_dataframe().to_hdf(
            filename,
            f"{sb_path}/{astro_def['metadata']}",
            mode="a",
        )
        return history.generate_history_row()

    def _write_cal_dataset(self, filename: Path):
        # TODO: Write the cal dataset
        return history.generate_history_row()

    def write(self, filename: str | Path, overwrite: bool = False) -> list[pd.DataFrame]:
        """Write the dataset to a file

        Args:
            filename (str | Path): The filename to write to
            overwrite (bool, optional): Overwrite the file if it exists. Defaults to False.

        Raises:
            FileExistsError: The file exists and overwrite is False

        """
        if isinstance(filename, str):
            filename = Path(filename)
        if filename.exists() and not overwrite:
            raise FileExistsError(f"{filename} already exists")

        astro_hist = self._write_astronomy_dataset(filename)
        cal_hist = self._write_cal_dataset(filename)

        return [astro_hist, cal_hist, history.generate_history_row()]


@dataclass
class Beam:
    """An SDHDF beam data object

    Args:
        label (str): The beam label
        filename (Path): The SDHDF file
        definition (dict): The SDHDF definition
        in_memory (bool, optional): Load data into memory. Defaults to False.
        client (Client, optional): Dask client. Defaults to None.

    Attributes:
        subbands (list[SubBand]): A list of subbands

    Methods:
        plot_waterfall: Plot a waterfall plot of the data
        plot_spectrum: Plot a spectrum of the data
        plot_wide: Plot spectra from all subbands

    """

    label: str
    filename: Path
    definition: dict
    in_memory: bool = False
    client: Client | None = None

    def __post_init__(self):
        meta_def = self.definition["beam"]["metadata"]
        with h5py.File(self.filename, "r") as f:
            sb_avail = Table.read(f, path=self.label + f"/{meta_def['band_params']}")
            self.subbands = [
                SubBand(
                    label=sb,
                    filename=self.filename,
                    definition=self.definition,
                    beam_label=self.label,
                    in_memory=self.in_memory,
                    client=self.client,
                )
                for sb in sb_avail["LABEL"]
            ]
            for sb in self.subbands:
                self.__dict__[sb.label] = sb

    def plot_waterfall(
        self,
        subband: int | str,
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(subband, int):
            subband = self.subbands[subband]
        elif isinstance(subband, str):
            subband = self.__dict__[subband]

        ax = subband.plot_waterfall(
            polarization=polarization,
            # bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_spectrum(
        self,
        subband: int | str,
        time: int = 0,
        polarization: int = 0,
        # bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(subband, int):
            subband = self.subbands[subband]
        elif isinstance(subband, str):
            subband = self.__dict__[subband]

        ax = subband.plot_spectrum(
            time=time,
            polarization=polarization,
            # bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_wide(
        self,
        time: int = 0,
        polarization: int = 0,
        # bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        fig, ax = plt.subplots()
        for i, sb in enumerate(self.subbands):
            sb.plot_spectrum(
                time=time,
                polarization=polarization,
                # bin=bin,
                flag=flag,
                ax=ax,
                label=sb.label,
                **plot_kwargs,
            )
        ax.legend()
        return ax

    def autoflag(self, sigma=3, n_windows=100) -> list[pd.DataFrame]:
        """Automatic flagging using rolling sigma clipping"""
        hists = []
        for sb in tqdm(self.subbands, desc="Flagging subbands"):
            hist = sb.autoflag(
                sigma=sigma,
                n_windows=n_windows,
            )
            hists.append(hist)
        return hists

    def decimate(
        self, bins: float | int, axis: str = "frequency", use_median: bool = False
    ) -> list[pd.DataFrame]:
        """Decimate the data

        Args:
            bins (float | int): If int, the number of channels to bin in an average.
                If float, the desired width of a channel after averaging.
            axis (str, optional): The axis to decimate along. Defaults to "frequency".
            use_median (bool, optional): Use the median instead of the mean. Defaults to False.

        Returns:
            list[pd.DataFrame]: list of history rows
        """
        hists = []
        for sb in tqdm(self.subbands, desc="Decimating subbands"):
            hist = sb.decimate(
                bins=bins,
                axis=axis,
                use_median=use_median,
            )
            hists.append(hist)
        return hists

    def write(self, filename: str | Path, overwrite: bool = False) -> list[pd.DataFrame]:
        """Write the data to a new file

        Args:
            filename (str | Path): The filename to write to
            overwrite (bool, optional): Overwrite the file if it exists. Defaults to False.

        Returns:
            list[pd.DataFrame]: list of history rows
        """
        hists = []
        for sb in tqdm(self.subbands, "Writing subbands"):
            hists.extend(sb.write(filename, overwrite=overwrite))

        return hists + [history.generate_history_row()]


@dataclass
class SDHDF:
    """An SDHDF data object

    Args:
        filename (Path): Path to the SDHDF file
        in_memory (bool, optional): Load data into memory. Defaults to False.
        parallel (bool, optional): Use dask for parallel processing. Defaults to False.

    Attributes:
        metadata (MetaData): Observation metadata
        beams (list[Beam]): list of beams

    Methods:
        plot_waterfall: Waterfall plot of the data
        plot_spectrum: Spectrum plot of the data
        plot_wide: Plot spectra from all subbands
        print_obs_metadata: list the observation metadata in the file
        write: Write the data to a new file

    """

    filename: Path
    in_memory: bool = False
    parallel: bool = False

    def __post_init__(self):
        self.client = Client() if self.parallel else None
        if self.parallel:
            logger.info(f"Dask dashboard at: {self.client.dashboard_link}")
        self.metadata = MetaData(self.filename)
        self.definition = self.metadata.definition
        with h5py.File(self.filename, "r") as f:
            keys = list(f.keys())
            #self.attrs = list(f.attrs) # TODO FIX THIS
            self.beams = [
                Beam(
                    label=key,
                    filename=self.filename,
                    in_memory=self.in_memory,
                    definition=self.definition,
                    client=self.client,
                )
                for key in keys
                if "beam_" in key
            ]
            for beam in self.beams:
                self.__dict__[beam.label] = beam

    def plot_waterfall(
        self,
        beam: int | str,
        subband: int | str,
        polarization: int = 0,
        # bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        """Waterfall plot of the data

        Args:
            beam (int | str): Beam to select.
            subband (int | str): Subband to select.
            polarization (int, optional): Polarization to select. Defaults to 0.
            bin (int, optional): Bin to select. Defaults to 0.
            flag (bool, optional): Blank flagged data. Defaults to False.
        """
        if isinstance(beam, int):
            beam = self.beams[beam]
        elif isinstance(beam, str):
            beam = self.__dict__[beam]
        ax = beam.plot_waterfall(
            subband=subband,
            polarization=polarization,
            # bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_spectrum(
        self,
        beam: int | str,
        subband: int | str,
        time: int = 0,
        polarization: int = 0,
        # bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(beam, int):
            beam = self.beams[beam]
        elif isinstance(beam, str):
            beam = self.__dict__[beam]

        ax = beam.plot_spectrum(
            subband=subband,
            time=time,
            polarization=polarization,
            # bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_wide(
        self,
        beam: int | str,
        time: int = 0,
        polarization: int = 0,
        # bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(beam, int):
            beam = self.beams[beam]
        elif isinstance(beam, str):
            beam = self.__dict__[beam]

        ax = beam.plot_wide(
            time=time,
            polarization=polarization,
            # bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def print_obs_metadata(self, format: str = "grid"):
        self.metadata.print_obs_metadata(format=format)

    def print_obs_config(self, format: str = "grid"):
        self.metadata.print_obs_config(format=format)

    #def print_attributes(self):
    #    self.metadata.print_attributes()
    #    #self.attributes.print_attributes()

    def flag_persistent_rfi(self):
        """Flag persistent RFI in all subbands."""
        telescope = self.metadata.primary_header["TELESCOPE"][0]
        rfi = flagging.get_persistent_rfi(telescope=telescope)
        for i, x in tqdm(
            rfi.iterrows(), desc="Flagging persistent RFI", total=len(rfi)
        ):
            for beam in self.beams:
                for sb in beam.subbands:
                    sb.astronomy_dataset.flag.loc[
                        dict(frequency=slice(x["freq0 MHz"], x["freq1 MHz"]))
                    ] = 1
        row = history.generate_history_row()
        self.metadata.history = pd.concat([self.metadata.history, row])

    def auto_flag_rfi(
        self,
        sigma=3,
        n_windows=100,
    ):
        """Automatic flagging using rolling sigma clipping"""
        self.flag_persistent_rfi()
        hists = []
        for beam in tqdm(self.beams, desc="Flagging beams"):
            hists.extend(beam.autoflag(sigma=sigma, n_windows=n_windows))

        self.metadata.history = pd.concat([self.metadata.history] + hists)

    def decimate(
        self, bins: float | int, axis: str = "frequency", use_median=False
    ):
        """Decimate the data in all subbands.

        Args:
            bins (float | int): If int, the number of channels to bin in an average.
                If float, the desired width of a channel after averaging.
            axis (str, optional): Axis to decimate along. Defaults to 'frequency'.
            use_median (bool, optional): Use median instead of mean. Defaults to False.

        """
        hists = []
        for beam in tqdm(self.beams, desc="Decimating beams"):
            hists.extend(beam.decimate(bins=bins, axis=axis, use_median=use_median))

        self.metadata.history = pd.concat([self.metadata.history] + hists)

    def write(self, filename: str | Path, overwrite: bool = False):
        """Write the SDHDF object to a file.

        Args:
            filename (Path): Filename to write to.
        """
        hists = []
        for beam in tqdm(self.beams, desc="Writing beams"):
            hists.extend(beam.write(filename, overwrite=overwrite))

        self.metadata.history = pd.concat([self.metadata.history] + hists + [history.generate_history_row()])
        self.metadata.write(filename, overwrite=overwrite)
