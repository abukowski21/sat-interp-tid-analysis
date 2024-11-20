import numpy as np
import xarray as xr
import sys
utils = '/Users/prasoonv/repo/sat-interp-tid-analysis/'
sys.path.append(f'{utils}Models/new/GITM/srcPython')
import gitm


def gitm_to_xarray(gitm_bin):
    # Extract attributes from the GITM binary file
    attrs = gitm_bin.attrs
    nLon = attrs['nLon']  # Number of longitude grid points
    nLat = attrs['nLat']  # Number of latitude grid points
    nAlt = attrs['nAlt']  # Number of altitude grid points
    nVars = attrs['nVars']  # Number of variables

    # Extract data and reshape it
    data_vars = {}
    for key in list(gitm_bin.keys())[1:]:  # Skip the first key if it's not a data variable
        data = np.array(gitm_bin[key])  # Convert data to a NumPy array
        data = data.reshape((nLon, nLat, nAlt))  # Reshape data to match the grid dimensions
        data_vars[key] = (['lon', 'lat', 'alt'], data)  # Store data with corresponding coordinates

    # Define coordinates for the xarray.Dataset
    coords = {
        'lon': gitm_bin['Longitude'][:, 0, 0],  # Longitude values
        'lat': gitm_bin['Latitude'][0, :, 0],  # Latitude values
        'alt': gitm_bin['Altitude'][0, 0, :],  # Altitude values
    }

    # Create an xarray.Dataset with the data variables and coordinates
    ds = xr.Dataset(data_vars, coords=coords, attrs=attrs)
    return ds