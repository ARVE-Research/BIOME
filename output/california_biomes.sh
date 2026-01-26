#!/usr/bin/env bash

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

ne=/work/kaplan_lab/datasets/naturalearth

ocean=$ne/ne_10m_ocean_blocks.gmt
rivers=$ne/ne_10m_rivers_lake_centerlines.gmt
lakes=$ne/ne_10m_lakes.gmt

cpt=/work/kaplan_lab/projects/may/BIOME/output/biome17.cpt

scale=1:4e6

infile=${1}

tmp=${infile##*/}

output=${tmp%%.*}_biomes.ps

title="BIOME1 Output"

boundsp=`gmt grdinfo -Ir $infile?biome`

boundsu=`gmt grdinfo -Io $infile?biome`

output=california_biomes.ps

gmt psbasemap $boundsp -Jx$scale -B0 -P -K > $output

gmt grdimage -R -J $infile?biome -C$cpt -nn -O -P -K >> $output

gmt psbasemap $boundsu+ue -Ja-100/50/$scale -B0 -O -P -K >> $output

gmt psxy $ocean -R -J -Gslategray1 -O -P -K >> $output

gmt psxy $rivers -R -J -Wthin,slategray1 -O -P -K >> $output

gmt psxy $lakes -R -J -Gslategray1 -O -P -K >> $output

gmt psbasemap -R -J -Ba2 -O -P -K >> $output

# legend
gmt gmtset FONT_ANNOT_PRIMARY 6p,Helvetica,black
sed "s/TITLE/$title/g" biome17.legend | gmt pslegend -Dx0.2i/0.2i+jBL+l1.5+w2.5i -F+pthinnest,black+gwhite -O -P >> $output


# Convert to PDF
gmt psconvert -A -Tf -Z $output
