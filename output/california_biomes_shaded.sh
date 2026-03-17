#!/usr/bin/env bash

# NB THIS SCRIPT WILL NOT RUN ON THE ARC HEAD NODE. SALLOC TO A COMPUTE NODE TO RUN

# module load GMT

penwid=thinner

gmt gmtset GMT_VERBOSE normal
gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_FRAME_PEN $penwid,black
gmt gmtset PS_MEDIA a2
gmt gmtset FONT_ANNOT_PRIMARY 10p,Helvetica,black
gmt gmtset FONT_LABEL 10p,Helvetica,black
gmt gmtset FONT_TITLE 10p,Helvetica,black
gmt gmtset FORMAT_GEO_MAP ddd:mmG

gmt gmtset MAP_TICK_LENGTH_PRIMARY 5p
gmt gmtset MAP_TICK_LENGTH_SECONDARY 3p
gmt gmtset MAP_ANNOT_OFFSET 2p
gmt gmtset MAP_TICK_PEN $penwid,black

# use this for local machine
ne=/Users/maycolgan/Desktop/Calgary/datasets

ocean=$ne/ne_10m_ocean_blocks.gmt
rivers=$ne/ne_10m_rivers_lake_centerlines.gmt
lakes=$ne/ne_10m_lakes.gmt
minorlakes=$ne/ne_10m_lakes_north_america.gmt
states=$ne/ne_10m_admin_1_states_provinces_lines.gmt
countries=$ne/ne_10m_admin_0_boundary_lines_land.gmt
relief=/Users/maycolgan/Desktop/Calgary/BIOME/output/relief_ca_light.nc

cpt=/Users/maycolgan/Desktop/Calgary/BIOME/output/biome17.cpt

# use this for cluster

# ne=/work/kaplan_lab/datasets/naturalearth
# 
# ocean=$ne/ne_10m_ocean_blocks.gmt
# rivers=$ne/ne_10m_rivers_lake_centerlines.gmt
# lakes=$ne/ne_10m_lakes.gmt
# 
# cpt=/work/kaplan_lab/projects/may/BIOME/output/biome17.cpt

# ---------------------------

scale=1:4e6

infile=${1}

tmp=${infile##*/}

output=${tmp%%.*}_biomes.ps

boundsp=$(gmt grdinfo -Ir $infile?biome)
boundsu=$(gmt grdinfo -Io $infile?biome)

output=california_biomes.ps

gmt psbasemap $boundsp -Jx$scale -B0 -P -K > $output

gmt grdimage -R -J $infile?biome -C$cpt -I$relief -nn -O -P -K >> $output

gmt psbasemap $boundsu+ue -Ja-115/45/$scale -B0 -O -P -K >> $output

gmt psxy $ocean -R -J -Gslategray1 -O -P -K >> $output

gmt psxy $rivers -R -J -Wthin,slategray1 -O -P -K >> $output

gmt psxy $lakes -R -J -Gslategray1 -O -P -K >> $output

gmt psxy $minorlakes -R -J -Gslategray1 -O -P -K >> $output

gmt psxy $states -R -J -W0.5p,black -O -P -K >> $output

gmt psxy $countries -R -J -W1p,black -O -P -K >> $output

# plot sample points (in geographic coordinates)
gmt psxy -R -J -Sa0.15i -Gyellow -Wthin,black -O -P -K << EOF >> $output
-124.024247126 41.3654057278
-119.909967252 38.5227851773
-122.002997515 39.6797625379
-119.734905687 37.1122185252
-116.836398522 36.2455129827

EOF

gmt pstext -R -J -F+f11p,Helvetica-Bold,black+jLM -D0.1i/0 -O -P -K << EOF >> $output
-116.836398522 36.2455129827 Badwater Basin
-119.909967252 38.5227851773 Carson Pass
-124.024247126 41.3654057278 Redwood NP
-122.002997515 39.6797625379 US-RGo
-119.734905687 37.1122185252 US-xSJ
EOF


gmt psbasemap -R -J -Ba2 -O -P -K >> $output

# legend
gmt gmtset FONT_ANNOT_PRIMARY 8p,Helvetica,black
sed "s/TITLE/$title/g" biome17.legend | gmt pslegend -Dx0.2i/0.2i+jBL+l1.5+w2.5i -F+pthinnest,black+gwhite -O -P >> $output


# Convert to PDF
gmt psconvert -A -Tf -Z $output