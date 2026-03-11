#!/usr/bin/env bash

# Map June PET from BIOME output (california.nc)
# Plots mpet for month 6 (index 5, zero-based)

# NB THIS SCRIPT WILL NOT RUN ON THE ARC HEAD NODE. SALLOC TO A COMPUTE NODE TO RUN

# module load GMT

penwid=thinner

gmt gmtset GMT_VERBOSE normal
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_FRAME_PEN $penwid,black
gmt gmtset PS_MEDIA a2
gmt gmtset FONT_ANNOT_PRIMARY 8p,Helvetica,black
gmt gmtset FONT_LABEL 8p,Helvetica,black
gmt gmtset FONT_TITLE 9p,Helvetica,black
gmt gmtset FORMAT_GEO_MAP ddd:mmG

gmt gmtset MAP_TICK_LENGTH_PRIMARY 5p
gmt gmtset MAP_TICK_LENGTH_SECONDARY 3p
gmt gmtset MAP_ANNOT_OFFSET 2p
gmt gmtset MAP_TICK_PEN $penwid,black

# --- paths ---

# use this for local machine
ne=/Users/maycolgan/Desktop/Calgary/datasets

# use this for cluster
# ne=/work/kaplan_lab/datasets/naturalearth

ocean=$ne/ne_10m_ocean_blocks.gmt
rivers=$ne/ne_10m_rivers_lake_centerlines.gmt
lakes=$ne/ne_10m_lakes.gmt

infile=${1:-california.nc}

# --- make color palette (equal intervals, roma) ---
# adjust -T range to suit your data: min/max/step

gmt makecpt -Cturbo -T100/400/20 -Z > junepet.cpt

# --- map setup ---

scale=1:4e6

boundsp=$(gmt grdinfo -Ir $infile?mpet)
boundsu=$(gmt grdinfo -Io $infile?mpet)

output=california_junepet.ps

gmt psbasemap $boundsp -Jx$scale -B0 -P -K > $output

gmt grdimage -R -J "$infile?mpet[5]" -Cjunepet.cpt -nn -O -P -K >> $output

gmt psbasemap $boundsu+ue -Ja-115/45/$scale -B0 -O -P -K >> $output

gmt psxy $ocean -R -J -Gslategray1 -O -P -K >> $output
gmt psxy $rivers -R -J -Wthin,slategray1 -O -P -K >> $output
gmt psxy $lakes -R -J -Gslategray1 -O -P -K >> $output

gmt psbasemap -R -J -Ba2 -O -P -K >> $output

# --- color bar ---

gmt gmtset FONT_ANNOT_PRIMARY 12p,Helvetica,black
gmt gmtset FONT_LABEL 13p,Helvetica,black
    
gmt psscale -Cjunepet.cpt -Dx0.25i/0.6i+jBL+w4i/0.2i+h \
    -Baf+l"Modeled June PET (mm/day)" \
    -F+pthinnest,black+gwhite+c0.1i -O -P >> $output
# --- convert to PDF ---

gmt psconvert -A -Tf -Z $output

# --- clean up ---

rm -f junepet.cpt gmt.conf gmt.history