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

boundsp=$(gmt grdinfo -Ir $infile?biome)
boundsu=$(gmt grdinfo -Io $infile?biome)

output=california_biomes.ps

gmt psbasemap $boundsp -Jx$scale -B0 -P -K > $output

gmt grdimage -R -J $infile?biome -C$cpt -nn -O -P -K >> $output

gmt psbasemap $boundsu+ue -Ja-115/45/$scale -B0 -O -P -K >> $output

gmt psxy $ocean -R -J -Gslategray1 -O -P -K >> $output

gmt psxy $rivers -R -J -Wthin,slategray1 -O -P -K >> $output

gmt psxy $lakes -R -J -Gslategray1 -O -P -K >> $output


# plot sample points (in geographic coordinates)
gmt psxy -R -J -Sa0.15i -Gyellow -Wthin,black -O -P -K << EOF >> $output
-124.024247126 41.3654057278
-122.274890264 40.5653406931
-117.185297647 36.4104758873
-119.909967252 38.5227851773
-119.758069727 36.6048553944
-118.2696034 33.9477746502
-121.4251 35.9358
-120.5514 35.2858

EOF

# label the points
gmt pstext -R -J -F+f9p,Helvetica,black+jLM -D0.1i/0 -O -P -K << EOF >> $output
-124.024247126 41.3654057278 Eureka
-122.274890264 40.5653406931 Redding
-117.185297647 36.4104758873 Death Valley
-119.909967252 38.5227851773 Carson Pass
-119.758069727 36.6048553944 Fresno
-118.2696034 33.9477746502 Los Angeles
-121.4251 35.9358 Big Sur
-120.5514 35.2858 SLO inland
EOF

gmt psbasemap -R -J -Ba2 -O -P -K >> $output

# legend
gmt gmtset FONT_ANNOT_PRIMARY 6p,Helvetica,black
sed "s/TITLE/$title/g" biome17.legend | gmt pslegend -Dx0.2i/0.2i+jBL+l1.5+w2.5i -F+pthinnest,black+gwhite -O -P >> $output


# Convert to PDF
gmt psconvert -A -Tf -Z $output