from geopy.distance import geodesic
import math

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
