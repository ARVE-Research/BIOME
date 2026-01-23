#!/bin/bash
# Plot biomes from BIOME1 model output using native projected coordinates
# Usage: ./plot_biomes_proj.sh <netcdf_file>
#
# Projection info from climate file input
#   Lambert Azimuthal Equal Area
#   latitude_of_projection_origin = 50
#   longitude_of_projection_origin = -100
#   WGS84 ellipsoid

# Check if input file provided
if [ -z "$1" ]; then
    echo "Usage: $0 <netcdf_file>"
    exit 1
fi

penwid=thinner

# GMT settings
gmt gmtset GMT_VERBOSE normal
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_FRAME_PEN $penwid,black
gmt gmtset PS_MEDIA a2
gmt gmtset FONT_ANNOT_PRIMARY 8p,Helvetica,black
gmt gmtset FONT_LABEL 8p,Helvetica,black
gmt gmtset FONT_TITLE 9p,Helvetica,black
gmt gmtset MAP_TICK_LENGTH_PRIMARY -5p
gmt gmtset MAP_TICK_LENGTH_SECONDARY -3p
gmt gmtset MAP_ANNOT_OFFSET -2p
gmt gmtset MAP_TICK_PEN $penwid,black

# -----
infile=${1}
nopath=${infile##*/}
output=${nopath%%.*}_biomes_proj.ps


echo "Plotting $output from $infile"

title="BIOME1 Output"

# --------
# Get the bounds directly from the biome variable (in meters)
# --------
bounds=$(gmt grdinfo -Ir ${infile}?biome)
echo "Projected bounds (meters): $bounds"

# --------
# Create the map using linear projection (data already projected)
# --------
echo "Creating map..."

# Extract xmin, xmax, ymin, ymax from bounds string
xmin=$(echo $bounds | sed 's/-R//' | cut -d'/' -f1)
xmax=$(echo $bounds | sed 's/-R//' | cut -d'/' -f2)
ymin=$(echo $bounds | sed 's/-R//' | cut -d'/' -f3)
ymax=$(echo $bounds | sed 's/-R//' | cut -d'/' -f4)

# Calculate aspect ratio
width_km=$(echo "scale=2; ($xmax - $xmin) / 1000" | bc)
height_km=$(echo "scale=2; ($ymax - $ymin) / 1000" | bc)
aspect=$(echo "scale=4; $height_km / $width_km" | bc)

map_width=6
map_height=$(echo "scale=2; $map_width * $aspect" | bc)

echo "Map dimensions: ${map_width}i x ${map_height}i"

proj="-JX${map_width}i/${map_height}i"

# Initialize PostScript
gmt psbasemap $bounds $proj -B0 -P -K > $output

# Extract biome data and plot with psxy (grdimage is crashing)
echo "Extracting data for plotting..."
gmt grd2xyz ${infile}?biome | awk '$3 != "NaN" && $3 > 0 && $3 <= 17 {print $1, $2, $3}' > biome_proj.xyz

# Plot biomes
gmt psxy biome_proj.xyz $bounds $proj -Ss0.050i -Cbiome17.cpt -O -P -K >> $output

# Add frame
gmt psbasemap $bounds $proj -Bxaf+l"Easting (m)" -Byaf+l"Northing (m)" -BWSne+t"$title" -O -P -K >> $output

# Add legend
gmt gmtset FONT_ANNOT_PRIMARY 6p,Helvetica,black
sed "s/TITLE/$title/g" biome17.legend | gmt pslegend -Dx0.2i/0.2i+jBL+l1.5+w2.5i -F+pthinnest,black+gwhite -O -P >> $output

# Convert to PDF
gmt psconvert -A -Tf -Z $output

# Clean
rm -f gmt.conf gmt.history biome_proj.xyz

echo "Done! Output: ${nopath%%.*}_biomes_proj.pdf"
