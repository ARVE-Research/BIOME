#!/usr/bin/env bash
# this script plots soil water and evapotranspiration variables using the Generic Mapping Tools (Classic Mode)

gmt gmtset PS_MEDIA a3
gmt gmtset FORMAT_DATE_MAP o
gmt gmtset FORMAT_TIME_PRIMARY_MAP Character
gmt gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black
gmt gmtset FONT_LABEL 11p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 13p,Helvetica-Bold,black

infile=${1}
output=${place_clean}_soils.ps

# Get lat/lon from .nc file
lon=$(gmt convert singlepixel.nc?lon | tail -1)
lat=$(gmt convert singlepixel.nc?lat | tail -1)

# Reverse geocode to get place name
place=$(curl -s "https://nominatim.openstreetmap.org/reverse?format=json&lat=${lat}&lon=${lon}&zoom=10" | \
        grep -o '"display_name":"[^"]*' | \
        cut -d'"' -f4 | \
        cut -d',' -f1-2)

if [ -z "$place" ]; then
    place="Unknown location"
fi

echo "Plotting data for: ${place} (Lon: ${lon}, Lat: ${lat})"

# get bounds - need columns 15 (soilw) and 17 (pet)
bounds=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg\n",$1,$2,$15,$17}' $infile | gmt gmtinfo -I1/10/1 -C) )

t0=${bounds[0]}
t1=${bounds[1]}
soilw_max=${bounds[3]}
pet_max=${bounds[5]}

# round up soilw_max to nearest 50
soilw_max=$(echo "scale=0; ((${soilw_max%.*}/50)+1)*50" | bc)

# round up aet_max to nearest 5
pet_max=$(echo "scale=0; ((${pet_max%.*}/5)+1)*5" | bc)

echo "Time range: $t0 to $t1"
echo "Soil water max: $soilw_max mm"
echo "PET max: $pet_max mm"

# ---
# GRAPH 1: Soil Water Content (mm)
# ---

gmt psbasemap -R$t0/$t1/0/$soilw_max -JX19/7 -Bpxa1O -Bpya50f10+l"Soil Water (mm)" -BWSen+t"Soil Water Content - ${place} (${lon}, ${lat})" -X4 -Y26 -P -K > $output

# soil water content bars (column 15)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$15}' $infile | gmt psxy -R -J -Sb0.02 -Gbisque3@30 -Wthin,bisque3 -O -P -K >> $output

# ---
# GRAPH 2: Actual Evapotranspiration (mm)
# ---

gmt psbasemap -R$t0/$t1/0/$pet_max -JX19/7 -Bpxa1O -Bpya+l"AET (mm)" -BWSen+t"Actual Evapotranspiration" -Y-10 -P -O -K >> $output

# AET bars (column 14)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$14}' $infile | gmt psxy -R -J -Sb0.02 -Ggold2@30 -Wthin,lightblue -O -P -K >> $output

# PET line (column 17)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$17}' $infile | gmt psxy -R -J -Wthick,violetred1 -O -P -K >> $output

# ---
# GRAPH 3: Relative Saturation and Alpha (both 0-1)
# ---

gmt psbasemap -R$t0/$t1/0/1.1 -JX19/7 -Bpxa1O -Bpya0.2f0.1+l"Fraction" -BWSen+t"Relative Saturation and AET/PET Ratio" -Y-10 -P -O -K >> $output

# relative saturation line (column 16)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$16}' $infile | gmt psxy -R -J -Wthick,springgreen3 -O -P -K >> $output

# alpha line (column 12)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$12}' $infile | gmt psxy -R -J -Wthick,gold2 -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w4c+o0.2c -F+p0.5p+gwhite -O -P << EOF >> $output
S 0.2c - 0.5c - thick,springgreen3 0.6c Relative Sat. (w/whc)
S 0.2c - 0.5c - thick,gold2 0.6c Alpha (AET/PET)
EOF

# ---
# Convert to PDF
# ---

gmt psconvert -A+m0.5c -Tf -Z $output

echo "Plotting complete. PDF saved as ${place_clean}_soils.pdf"