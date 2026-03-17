#!/usr/bin/env bash
# this script plots temperature/precip, ET, and snow using GMT Classic Mode

gmt gmtset PS_MEDIA a3
gmt gmtset FORMAT_DATE_MAP o
gmt gmtset FORMAT_TIME_PRIMARY_MAP Character
gmt gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black
gmt gmtset FONT_LABEL 11p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 13p,Helvetica-Bold,black

infile=dmet_-428047.90_-706770.64.txt
output=carsonpass_output.ps

title="Model Output: Carson Pass, Alpine County, CA"

# --- bounds for temperature/precip ---
bounds=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg %lg\n",$1,$2,$3,$4,$5}' $infile | gmt gmtinfo -I1/5/5/20 -C) )
t0=${bounds[0]}
t1=${bounds[1]}
tmin=${bounds[4]}
tmax=${bounds[3]}
pmax=${bounds[7]}

# Get minimum dewpoint and use it if lower than tmin
tdew_min=$(awk 'BEGIN{min=999} {if($13<min) min=$13} END{print int(min/5)*5 - 5}' $infile)
if (( $(echo "$tdew_min < $tmin" | bc -l) )); then
  tmin=$tdew_min
fi

# --- bounds for ET ---
pet_max=$(awk 'BEGIN{max=0} {if($17>max) max=$17} END{print int(max+1)}' $infile)

# --- bounds for snow ---
swe_max=$(awk 'BEGIN{max=0} {if($8>max) max=$8} END{print int(max/50+1)*50}' $infile)
snow_max=$(awk 'BEGIN{max=0} {if($6>max) max=$6; if($7>max) max=$7} END{print int(max+5)}' $infile)

# -------------------------------------------------------------------------
# GRAPH 1: Temperature and Precipitation
# -------------------------------------------------------------------------

gmt psbasemap -R$t0/$t1/$tmin/$tmax -JX19/8 -Bpxa1O -Bpya10f2+l"Temperature (C)" -BWSn+t"$title" -X4 -Y26 -P -K > $output

# tmin
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$29}' $infile | gmt psxy -R -J -Wthin,darkblue -O -P -K >> $output
# tmax
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$30}' $infile | gmt psxy -R -J -Wthin,darkred -O -P -K >> $output
# night
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$4}' $infile | gmt psxy -R -J -Wthin,blue -O -P -K >> $output
# day
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$3}' $infile | gmt psxy -R -J -Wthin,red -O -P -K >> $output
# dewpoint
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$13}' $infile | gmt psxy -R -J -Wthin,green,- -O -P -K >> $output

# precipitation on right axis
gmt psbasemap -R$t0/$t1/0/$pmax -JX19/8 -Bpya+l"Precipitation (mm)" -BE -O -P -K >> $output
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$5}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,dodgerblue -O -P -K >> $output

# legend (plotted last so it's on top of all data)
gmt pslegend -R -J -DjTR+w3c+o0.2c -F+p1p+gwhite -O -P -K >> $output << EOF
S 0.2c - 0.5c - thin,darkred 0.8c Tmax
S 0.2c - 0.5c - thin,red 0.8c Tday
S 0.2c - 0.5c - thin,blue 0.8c Tnight
S 0.2c - 0.5c - thin,darkblue 0.8c Tmin
S 0.2c - 0.5c - thin,green,- 0.8c Tdew
EOF

# -------------------------------------------------------------------------
# GRAPH 2: Actual and Potential Evapotranspiration
# -------------------------------------------------------------------------

gmt psbasemap -R$t0/$t1/0/$pet_max -JX19/6 -Bpxa1O -Bpya+l"mm" -BWSen+t"Actual and Potential Evapotranspiration" -Y-8 -P -O -K >> $output

# AET bars (column 14)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$14}' $infile | gmt psxy -R -J -Sb0.02 -Glightblue@30 -Wthin,lightblue -O -P -K >> $output

# PET line (column 17)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$17}' $infile | gmt psxy -R -J -Wthick,violetred1 -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w3.5c+o0.2c -F+p1p+gwhite -O -P -K >> $output << EOF
S 0.2c r 0.3c lightblue@30 thin,lightblue 0.8c AET
S 0.2c - 0.5c - thick,violetred1 0.8c PET
EOF

# -------------------------------------------------------------------------
# GRAPH 3: Snow Accumulation, Melt, and Albedo
# -------------------------------------------------------------------------

gmt psbasemap -R$t0/$t1/0/$snow_max -JX19/6 -Bpxa1O -Bpya+l"Snow (mm)" -BWSn+t"Daily Snow Accumulation, Melt, and Albedo" -Y-8 -P -O -K >> $output

# melt bars (transparent coral)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$7}' $infile | gmt psxy -R -J -Sb0.02 -Glightcoral@50 -Wthin,lightcoral -O -P -K >> $output

# accumulation bars (transparent blue)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$6}' $infile | gmt psxy -R -J -Sb0.02 -Gdodgerblue@50 -Wthin,dodgerblue -O -P -K >> $output

# albedo on right axis (0-1)
gmt psbasemap -R$t0/$t1/0/1 -JX19/6 -Bpya0.2f0.1+l"Albedo" -BE -O -P -K >> $output
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$11}' $infile | gmt psxy -R -J -Wthin,orange -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w3.5c+o0.2c -F+p1p+gwhite -O -P >> $output << EOF
S 0.2c r 0.3c dodgerblue@50 thin,dodgerblue 0.8c Snow accum.
S 0.2c r 0.3c lightcoral@50 thin,lightcoral 0.8c Snowmelt
S 0.2c - 0.5c - thin,orange 0.8c Albedo
EOF

# ---
# Convert to PDF
# ---

gmt psconvert -A+m0.5c -Tf -Z $output

echo "Plotting complete. PDF saved as carsonpass_output.pdf"
