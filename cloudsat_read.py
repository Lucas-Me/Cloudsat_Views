# -*- coding: utf-8 -*-
"""Functions to read and visualise CloudSat data."""
from datetime import datetime, timedelta
import h5py
import numpy as np


VARNAMES = ["Longitude", "Latitude", "Height", "Profile_time", "DEM_elevation"]


def get_geodata(h5name, varnames=VARNAMES, proftime2datetime=True, return_list=False):
    """Get geographic data from the CloudSat data."""
    with h5py.File(h5name, "r") as f:
        for ikey in f.keys():
            if "Data Fields" in f[ikey]:
                topkey = ikey
                break

        var_dict = {}
        for ivar in varnames:
            fld_dtype = f[topkey]["Geolocation Fields"][ivar][:][0][0].dtype
            var_dict[ivar] = f[topkey]["Geolocation Fields"][ivar][:].astype(fld_dtype)

        sw_attr = f[topkey]["Swath Attributes"]
        if "Profile_time" in var_dict:
            start_time = sw_attr["start_time"][0][0].decode("UTF-8")
            start_time = datetime.strptime(start_time, "%Y%m%d%H%M%S")
            if proftime2datetime:
                var_dict["Profile_time"] = np.array(
                    [
                        start_time + timedelta(seconds=float(i))
                        for i in var_dict["Profile_time"]
                    ]
                )  # noqa

        if "DEM_elevation" in var_dict:
            var_dict["DEM_elevation"][var_dict["DEM_elevation"] < 0] = 0.0

        if return_list:
            return [var_dict[i] for i in varnames]
        else:
            return var_dict


def read_data(h5name, data_field="Radar_Reflectivity", limits=None, fillmask=False):
    """Read CloudSat data from a HDF5 file and return a masked numpy array."""
    # TODO: return units f['2B-CWC-RO']['Swath Attributes']
    #                         ['RO_liq_water_content.units'][0]
    with h5py.File(h5name, "r") as f:
        for ikey in f.keys():
            if "Data Fields" in f[ikey]:
                topkey = ikey
                break

        data = f[topkey]["Data Fields"][data_field][:].astype(np.float32)

        sw_attr = f[topkey]["Swath Attributes"]

        missval = sw_attr[data_field + ".missing"][:][0][0]
        fillval = sw_attr["_FV_" + data_field][:][0][0]

        for imask in [missval, fillval]:
            data = np.ma.masked_equal(data, imask)

        if limits is not None:
            lim_mask = np.logical_or(data < limits[0], data > limits[1])
            data = np.ma.masked_equal(data, lim_mask)

        if fillmask:
            data = data.filled(np.nan)

        data_factor = sw_attr[data_field + ".factor"][0][0]
        data_offset = sw_attr[data_field + ".offset"][0][0]
        data = (data - data_offset) / data_factor

    return data


def read_cldclass(h5name):
    """
    Convert cloud scenario codes to one of the 8 classes

    disc.sci.gsfc.nasa.gov/measures/documentation/README.AIRS_CloudSat.pdf
    """
    with h5py.File(h5name, "r") as f:
        cldclass = f["2B-CLDCLASS/Data Fields/cloud_scenario"][:]
        cldclass = np.right_shift(np.array(cldclass), 1) & int("1111", 2)
        cloudcodes = {
            "0": "clear",
            "1": "Ci",
            "2": "As",
            "3": "Ac",
            "4": "St",
            "5": "Sc",
            "6": "Cu",
            "7": "Ns",
            "8": "DC",
        }

    return cldclass, cloudcodes