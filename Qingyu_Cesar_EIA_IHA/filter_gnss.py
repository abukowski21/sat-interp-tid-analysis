import pandas as pd
import argparse
import apexpy
from datetime import datetime

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--mlat_lim", type=float, default=30)
parser.add_argument("--mlon_min", type=float, default=350)  # change the filter in lines 36,37
parser.add_argument("--mlon_max", type=float, default=10)
args = parser.parse_args()

# Read ASCII file (skip headers automatically if Madrigal format)
df = pd.read_csv(args.filename, delim_whitespace=True, comment='#')

# Keep only relevant columns
df = df[["time", "gdlat", "glon", "tec"]]

# Convert to magnetic coordinates
apex = apexpy.Apex(datetime.utcnow().year)  # use current epoch
mlats, mlons = [], []
for lat, lon in zip(df["gdlat"], df["glon"]):
    mlat, mlon = apex.convert(lat, lon, "geo", "apex", height=350)
    mlats.append(mlat)
    mlons.append(mlon)

df["mlat"] = mlats
df["mlon"] = mlons

# Apply filters
df = df[
    (df["mlat"] >= -args.mlat_lim) &
    (df["mlat"] <= args.mlat_lim) &
    ((df["mlon"] >= args.mlon_min) |
    (df["mlon"] <= args.mlon_max))
]

# Save filtered file
outname = args.filename.replace(".txt", "_filtered.csv")
df.to_csv(outname, index=False)
