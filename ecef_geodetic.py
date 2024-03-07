
import math
import numpy as np

def ecef_to_lla(x, y, z):
    # WGS84 parameters
    a = 6378137.0  # semi-major axis in meters
    f = 1 / 298.257223563  # flattening
    b = (1 - f) * a  # semi-minor axis

    e = math.sqrt(1 - (b**2) / (a**2))  # eccentricity

    lon = math.atan2(y, x)
    p = math.sqrt(x**2 + y**2)

    # Iterative solution for latitude
    lat = math.atan2(z, p * (1 - e**2))

    # Iterative solution for altitude
    N = a / math.sqrt(1 - e**2 * math.sin(lat)**2)
    alt = p / math.cos(lat) - N

    # Iterate to refine latitude and altitude
    for _ in range(5):
        lat_new = math.atan2(z, p * (1 - e**2 * N / (N + alt)))
        N = a / math.sqrt(1 - e**2 * math.sin(lat_new)**2)
        alt = p / math.cos(lat_new) - N
        if abs(lat - lat_new) < 1e-10:
            break
        lat = lat_new

    # Convert radians to degrees
    lat = math.degrees(lat)
    lon = math.degrees(lon)

    return lat, lon, alt


def lla_to_ecef(lat, lon, alt):
    # WGS84 parameters
    a = 6378137.0  # semi-major axis in meters
    f = 1 / 298.257223563  # flattening
    b = (1 - f) * a  # semi-minor axis

    e = math.sqrt(1 - (b**2) / (a**2))  # eccentricity

    # Calculate N, the radius of curvature in the prime vertical
    N = a / math.sqrt(1 - e**2 * math.sin(math.radians(lat))**2)

    # Convert geodetic coordinates to ECEF
    x = (N + alt) * math.cos(math.radians(lat)) * math.cos(math.radians(lon))
    y = (N + alt) * math.cos(math.radians(lat)) * math.sin(math.radians(lon))
    z = (N * (1 - e**2) + alt) * math.sin(math.radians(lat))

    return x, y, z

def ut_to_lt(time_array, glon):
    """
    Compute local time from date and longitude.

    Parameters
    ----------
    time_array : array-like
        Array-like of datetime objects in universal time
    glon : array-like or float
        Float or array-like of floats containing geographic longitude
        in degrees. If single value or array of a different shape, all
        longitudes are applied to all times. If the shape is the same as
        `time_array`, the values are paired in the SLT calculation.

    Returns
    -------
    array of floats
        List of local times in hours
    """
    time_array = np.asarray(time_array)
    glon = np.asarray(glon)

    # Get UT seconds of day
    try:  # if numpy timestamps
        utsec = [(ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                  + ut.microsecond * 1.0e-6) / 3600.0 for ut in time_array]
    except BaseException:
        utsec = []
        for ut in time_array:
            ut = pd.Timestamp(ut)
            utsec.append((ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                          + ut.microsecond * 1.0e-6) / 3600.0)
    
    # Determine if the calculation is paired or broadcasted
    if glon.shape == time_array.shape:
        lt = np.array([utime + glon[i] / 15.0 for i,
                      utime in enumerate(utsec)])
    else:
        lt = np.array([utime + glon / 15.0 for utime in utsec])
    # Adjust to ensure that 0.0 <= lt < 24.0
    while np.any(lt < 0.0):
        lt[lt < 0.0] += 24.0
    
    while np.any(lt >= 24.0):
        lt[lt >= 24.0] -= 24.0
    
    return lt



def spher_cart(lat, lon, rad):
    z = rad*math.sin(math.radians(lat))
    x = rad*math.cos(math.radians(lat))*math.cos(math.radians(lon))
    y = rad*math.cos(math.radians(lat))*math.sin(math.radians(lon))
    return x, y, z

def cart_spher(x,y,z):
    r = math.sqrt(x*x + y*y + z*z)
    lon = math.degrees(math.atan(y/x))
    lat = math.degrees(math.atan(z/math.sqrt(y*y + x*x)))

    return lat, lon, r
