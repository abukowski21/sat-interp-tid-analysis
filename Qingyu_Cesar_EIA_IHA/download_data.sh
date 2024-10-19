#!/bin/bash
 
# Loop over years from 2000 to 2023
for year in {2000..2024}; do
    for month in march sept june dec; do
        
        # Set dates and phase based on the month
        if [ "$month" == "march" ]; then
            start_date="03/01/${year}"
            end_date="04/11/${year}"
            phase="equinox"
        elif [ "$month" == "june" ]; then
            start_date="05/30/${year}"
            end_date="07/11/${year}"
            phase="solstice"
        elif [ "$month" == "sept" ]; then
            start_date="09/02/${year}"
            end_date="10/13/${year}"
            phase="equinox"
        elif [ "$month" == "dec" ]; then
            start_date="12/01/${year}"
            end_date="01/31/$((year + 1))"
            phase="solstice"
        fi

        
        # Create the directory structure
        mkdir -p "/home/pxv220016/scratch/Qingyu_Cesar_EIA/${month}_data"
        mkdir -p "/home/pxv220016/scratch/Qingyu_Cesar_EIA/${month}_data/${year}_${month}_${phase}"

        # Run the globalDownload.py script with the specified dates
        globalDownload.py --verbose --url="http://cedar.openmadrigal.org" \
                          --outputDir="/home/pxv220016/scratch/Qingyu_Cesar_EIA/${month}_data/${year}_${month}_${phase}" \
                          --user_fullname="Prasoon" \
                          --user_email="pxv220016@utdallas.edu" \
                          --user_affiliation="None" \
                          --format="ascii" \
                          --startDate="$start_date" \
                          --endDate="$end_date" \
                          --inst=8000 \
                          --kindat=3500 

        # Unzip downloaded files in the output directory
        gunzip"/home/pxv220016/scratch/Qingyu_Cesar_EIA/${month}_data/${year}_${month}_${phase}/gps*gz"
        
    done
done