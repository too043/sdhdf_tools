#!/usr/bin/env python
from __future__ import annotations

import h5py
import matplotlib as mpl
import numpy as np
import pandas as pd
import pylab as plt
from astropy.table import QTable

__version__ = "2.2"
__author__ = "Lawrence Toomey"

# increase matplotlib chunk size above default -
# this improves speed and prevents Agg rendering failure
# when plotting large data sets
mpl.rcParams["agg.path.chunksize"] = 20000


def print_hdr(tb):
    """
    Format SDHDF header metadata output

    :param astropy.QTable tb: astropy.QTable metadata object
    :return: None
    """
    params = []
    for col in tb.colnames:
        params.append([col, tb[col][0]])
    df = pd.DataFrame(params, columns=("-- Key --", "-- Value --"))
    print(df)


def read_sdhdf_header(f_pth, dset_pth):
    """
    Read SDHDF header metadata

    :param string f_pth: Path to SDHDF file
    :param HDF dataset dset_pth: Path to HDF metadata dataset
    :return astropy.QTable tb: astropy.QTable metadata object
    """
    try:
        with h5py.File(f_pth, "r") as h5:
            tb = QTable.read(h5, path=dset_pth)
    except Exception as e:
        print("ERROR: failed to read file %s" % f_pth, e)

    return tb


def get_available_subbands(sb_dict):
    """
    Retrieve a list of available sub-band groups
    from SDHDF format file

    :param dict sb_dict: Dictionary of sub-bands to check
    :return dict sb_avail_dict: Dictionary of sub-bands and metadata
    """
    sb_avail_dict = {}
    for k, v in sb_dict.items():
        sb_avail_dict[v[2]] = [k, v[0], v[1]]

    return sb_avail_dict


def get_channel_range(freq_arr, c_freq, z_width):
    """
    Retrieve a range of frequency channels given
    a specific centre frequency and zoom window width

    :param numpy.ndarray freq_arr: Array containing the frequency axis
    :param int c_freq: User defined centre frequency of zoom window (MHz)
    :param float z_width: User defined half width of zoom window (MHz)
    :return int z_min: Zoom window frequency channel minimum
    :return int z_max: Zoom window frequency channel maximum
    """
    n_chan = freq_arr.shape[0]
    freq_min = int(freq_arr[0])
    freq_max = int(freq_arr[-1])
    freq_range = c_freq - freq_min
    bw = freq_max - freq_min
    ch_bw = bw / n_chan
    z_min = freq_range / ch_bw - (z_width / ch_bw)
    z_max = freq_range / ch_bw + (z_width / ch_bw)

    return int(z_min), int(z_max)


def plot_sdhdf(f, sp, wf):
    """
    Plot the spectra (uncalibrated flux vs frequency),
    and waterfall (time vs frequency) for each polarisation product,
    for a specified sub-band from an SDHDF format file,
    and zoom if specified by the user (Default is no zoom)

    :param string f: Name of SDHDF file to read
    :param bool sp: Plot spectra [True|False]
    :param bool wf: Plot waterfall [True|False]
    :return None
    """
    h5 = h5py.File(f, "r")

    bp = QTable.read(h5, path="/metadata/beam_params")
    ph = QTable.read(h5, path="/metadata/primary_header")

    # display primary header astropy.QTable object
    print("----------------------------------------------------")
    hdr = read_sdhdf_header(f, "/metadata/primary_header")
    print_hdr(hdr)

    # get header version
    hdr_ver = float(hdr["HDR_DEFN_VERSION"][0])

    # display available beam groups in file
    print("\n----------------------------------------------------")
    print("Available beams to plot are:\n")
    print(bp)

    # user selects beam to plot (need single quotes around input for python 2.7)
    beam = input("\nWhich beam do you wish to plot (e.g. 0) ?\n")
    beam_label = "beam_" + beam

    # display available sub-band groups in file
    print("\n----------------------------------------------------")
    print("Available sub-bands to plot are:\n")
    sb_avail = QTable.read(h5, path=beam_label + "/metadata/band_params")
    print(sb_avail)

    # user selects sub-band to plot (need single quotes around input for python 2.7)
    sb = input("\nWhich sub-band do you wish to plot (e.g. 0) ?\n")
    sb_label = "band_SB" + sb

    if sb_label in sb_avail["LABEL"]:
        # set paths to data
        sb_data = beam_label + "/" + sb_label + "/astronomy_data/data"
        sb_freq = beam_label + "/" + sb_label + "/astronomy_data/frequency"
        print("----------------------------------------------------")
        print("Processing sub-band: %r" % sb)
        # print(h5[sb_freq[:][0]])
        sb_min = h5[sb_freq][0]
        sb_max = h5[sb_freq][-1]
        # print(sb_min, sb_max)
        # exit()
        print("Freq min: %s Freq max: %s" % (str(sb_min), str(sb_max)))
        # exit()
        zoom = input("Zoom? (y/n)  ")
        if zoom == "y":
            zoom_centre = input("Enter zoom band centre frequency (integer MHz)  ")
            zoom_width = input("Enter zoom band window half width (float MHz)  ")
            z_window_lo = int(zoom_centre) - float(zoom_width)
            z_window_hi = int(zoom_centre) + float(zoom_width)

            for row in range(len(sb_avail)):
                if sb_avail["LABEL"][row] == sb_label:
                    # sb_min = sb_avail[row]['LOW_FREQ']
                    # sb_max = sb_avail[row]['HIGH_FREQ']
                    sb_min = h5[sb_freq][0]
                    sb_max = h5[sb_freq][-1]

            # check that the user values are within the specified sub-band
            if z_window_lo >= sb_min and z_window_hi <= sb_max:
                zoom = True
            else:
                raise ValueError(
                    "ERROR: input values are not within range of sub-band %s" % sb
                )
        else:
            zoom = False
        print("----------------------------------------------------")
    else:
        raise ValueError("ERROR: sub-band %s does not exist in data file" % sb)

    # sb_data = beam_label + '/' + sb_label + '/astronomy_data/data'
    # sb_freq = beam_label + '/' + sb_label + '/astronomy_data/frequency'

    # print(h5[sb_freq][0])
    # exit()

    if hdr_ver >= 2.1:
        n_pol = h5[sb_data].shape[1]
    else:
        n_pol = h5[sb_data].shape[2]
    op = QTable.read(h5, path=beam_label + "/" + sb_label + "/metadata/obs_params")

    if ph["CAL_MODE"] == "ON":
        sb_cal_data_on = beam_label + "/" + sb_label + "/calibrator_data/cal_data_on"
        sb_cal_data_off = beam_label + "/" + sb_label + "/calibrator_data/cal_data_off"
        sb_cal_freq = beam_label + "/" + sb_label + "/calibrator_data/cal_frequency"
        if hdr_ver >= 2.1:
            cal_n_pol = h5[sb_cal_data_on].shape[1]
        else:
            cal_n_pol = h5[sb_cal_data_on].shape[2]

    print("Data array shape: %r" % str(h5[sb_data].shape))
    print("Frequency array shape: %r" % str(h5[sb_freq].shape))
    print("Number of polarisations: %r" % n_pol)

    if ph["CAL_MODE"] == "ON":
        print(
            "Calibration data array shape (cal on): %r" % str(h5[sb_cal_data_on].shape)
        )
        print(
            "Calibration data array shape (cal off): %r"
            % str(h5[sb_cal_data_off].shape)
        )
        print("Calibration frequency array shape: %r" % str(h5[sb_cal_freq].shape))
        print("Number of polarisations (calibration data): %r" % cal_n_pol)
    print("----------------------------------------------------")

    # get range of channels for zoom window
    if zoom is True:
        z_min, z_max = get_channel_range(
            h5[sb_freq], int(zoom_centre), float(zoom_width)
        )
        print("Zoom enabled")
        print("Zoom window centred at %s MHz" % int(zoom_centre))
        print("Zoom window width: %s to %s MHz" % (z_window_lo, z_window_hi))
        print("Zoom channel min: %i" % z_min)
        print("Zoom channel max: %i" % z_max)
        print("----------------------------------------------------")

    # plot spectra
    if n_pol == 2 or n_pol == 4:
        plt.figure(figsize=(8, 8))
        color = ("b", "g", "r", "c")
        for ii in range(n_pol):
            print("Processing polarisation %r..." % ii)
            lc = color[ii]
            plt.subplot(n_pol, 1, ii + 1)
            print("Averaging over integrations\n")
            np.mean(h5[sb_data], axis=0)
            if ii <= 1:
                plt.yscale("log")
                plt.ylabel("Flux [Log counts]")
            else:
                plt.ylabel("Flux [counts]")
            if zoom is True:
                if hdr_ver >= 2.1:
                    data_arr = h5[sb_data][0, ii, z_min:z_max]
                else:
                    data_arr = h5[sb_data][0, 0, ii, z_min:z_max]
                plt.plot(
                    h5[sb_freq][z_min:z_max],
                    data_arr,
                    linewidth=1,
                    color=lc,
                    label="Pol" + str(ii),
                )
            else:
                if hdr_ver >= 2.1:
                    data_arr = h5[sb_data][0, ii, :]
                else:
                    data_arr = h5[sb_data][0, 0, ii, :]
                plt.plot(
                    h5[sb_freq][:],
                    data_arr,
                    linewidth=1,
                    color=lc,
                    label="Pol" + str(ii),
                )
            plt.legend(loc="upper right")
        plt.xlabel("Frequency [MHz]")

    elif n_pol == 1:
        plt.figure(figsize=(8, 8))
        if zoom is True:
            if hdr_ver >= 2.1:
                data_arr = h5[sb_data][0, :, z_min:z_max]
            else:
                data_arr = h5[sb_data][0, 0, :, z_min:z_max]
            plt.plot(
                h5[sb_freq][z_min:z_max],
                data_arr,
                linewidth=1,
                label="Pol" + str(n_pol),
            )
        else:
            if hdr_ver >= 2.1:
                data_arr = h5[sb_data][0, :, :]
            else:
                data_arr = h5[sb_data][0, 0, :, :]
            plt.plot(h5[sb_freq][:], data_arr, linewidth=1, label="Pol" + str(n_pol))
        plt.xlabel("Frequency [MHz]")
        plt.ylabel("Flux [counts]")

    plt.suptitle("Mean Flux - Sub-band %r" % sb, fontsize=14)
    plt.subplots_adjust(top=0.9)
    if wf is True:
        plt.show(block=False)
    else:
        plt.show()

    # plot waterfall data (time/frequency)
    if wf is True:
        plt.figure(figsize=(10, 8))
        plt.title("Waterfall - Sub-band %r" % sb, fontsize=14)
        if zoom is True:
            if hdr_ver >= 2.1:
                av = np.mean(h5[sb_data], axis=3)
                av_arr = av[:, 0, z_min:z_max]
            else:
                av = np.mean(h5[sb_data], axis=4)
                av_arr = av[:, 0, 0, z_min:z_max]

            if len(op["MJD"]) > 1:
                plt.imshow(
                    av_arr,
                    aspect="auto",
                    extent=(z_window_lo, z_window_hi, op["MJD"][0], op["MJD"][-1]),
                    interpolation="nearest",
                )
            else:
                plt.imshow(
                    av_arr,
                    aspect="auto",
                    extent=(z_window_lo, z_window_hi, 0, len(av_arr)),
                    interpolation="nearest",
                )

        else:
            if hdr_ver >= 2.1:
                av = np.mean(h5[sb_data], axis=3)
                av_arr = av[:, 0, :]
            else:
                av = np.mean(h5[sb_data], axis=4)
                av_arr = av[:, 0, 0, :]

            if len(op["MJD"]) > 1:
                plt.imshow(
                    av_arr,
                    aspect="auto",
                    extent=(
                        h5[sb_freq][0],
                        h5[sb_freq][-1],
                        op["MJD"][0],
                        op["MJD"][-1],
                    ),
                    interpolation="nearest",
                )
            else:
                plt.imshow(
                    av_arr,
                    aspect="auto",
                    extent=(h5[sb_freq][0], h5[sb_freq][-1], 0, len(av_arr)),
                    interpolation="nearest",
                )

        plt.xlabel("Frequency [MHz]")
        plt.ylabel("Time [MJD]")
        cbar = plt.colorbar()
        cbar.set_label("counts")
        plt.show()


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--filename", help="Path to SDHDF file to read", required=True)
    ap.add_argument(
        "--plot_spectra",
        help="Plot spectra (flux vs. frequency)",
        required=True,
        default=True,
    )
    ap.add_argument(
        "--plot_waterfall",
        help="Plot waterfall (time vs. frequency)",
        required=True,
        default=True,
    )
    args = ap.parse_args()

    plot_sdhdf(args.filename, bool(args.plot_spectra), bool(args.plot_waterfall))
