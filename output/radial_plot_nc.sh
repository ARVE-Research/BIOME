# Radial plot of annual snow cycle from NetCDF

ncfile=california.nc
output=radial_snow.ps

gmt gmtset PS_MEDIA a4
gmt gmtset FONT_ANNOT_PRIMARY 10p,Helvetica,black
gmt gmtset FONT_LABEL 12p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 14p,Helvetica-Bold,black

# Extract data - specify to read all months for the single x,y point
gmt convert $ncfile?snow[0:11][0][0] > snow.txt
gmt convert $ncfile?mpet[0:11][0][0] > pet.txt
gmt convert $ncfile?swe[0:11][0][0] > swe.txt

# Convert month (1-12) to theta (0-330 degrees, 30° per month)
# Calculate monthly totals for snow and PET

# Calculate monthly total snowfall and PET (mpet already monthly total)
# First calculate snow monthly totals
awk '{
    month = NR
    days = (month==2) ? 28 : ((month==4||month==6||month==9||month==11) ? 30 : 31)
    snow_monthly = $1 * days
    print snow_monthly
}' snow.txt > snow_monthly.txt

# snow, PET, and SWE side by side
paste snow_monthly.txt pet.txt swe.txt > monthly_totals.txt

# Create polar coordinates for monthly total snowfall
# January at top (90°), going clockwise: 90, 60, 30, 0, 330, 300...
awk 'BEGIN {month=1} 
     {
       theta = 90 - (month-1)*30
       if (theta < 0) theta += 360
       print theta, $1
       month++
     }' monthly_totals.txt > snow_polar.txt
# Close polygon back to January value (not 0)
head -1 snow_polar.txt | awk '{print $1, $2}' >> snow_polar.txt

# Create polar coordinates for PET
awk 'BEGIN {month=1}
     {
       theta = 90 - (month-1)*30
       if (theta < 0) theta += 360
       print theta, $2
       month++
     }' monthly_totals.txt > pet_polar.txt
# Close the circle
head -1 pet_polar.txt >> pet_polar.txt

# Create polar coordinates for SWE (average daily SWE for the month)
awk 'BEGIN {month=1}
     {
       theta = 90 - (month-1)*30
       if (theta < 0) theta += 360
       days = (month==2) ? 28 : ((month==4||month==6||month==9||month==11) ? 30 : 31)
       avg_daily_swe = $3 / days
       print theta, avg_daily_swe
       month++
     }' monthly_totals.txt > swe_polar.txt
# Close the circle
head -1 swe_polar.txt >> swe_polar.txt

# Get location
lon=$(gmt convert $ncfile?lon[0][0] | tail -1)
lat=$(gmt convert $ncfile?lat[0][0] | tail -1)

# Create polar plot with proper annotations
gmt psbasemap -R0/360/0/400 -JP7i -Bxa30g30 -Bya50g50 -BWsne+t"Annual Snow Cycle - Lon: ${lon}, Lat: ${lat}" -P -K > $output

# Plot monthly total snowfall as line only (no fill)
gmt psxy snow_polar.txt -R -J -W2.5p,blue -O -K >> $output

# Plot monthly PET as line
gmt psxy pet_polar.txt -R -J -W2.5p,orange -O -K >> $output

# Plot SWE as line
gmt psxy swe_polar.txt -R -J -W2.5p,purple -O -K >> $output

# Add center point
echo "0 0" | gmt psxy -R -J -Sc0.15i -Gblack -O -K >> $output

# Add month labels on the plot at radius 135
cat > months.txt << EOF
90 135 JAN
60 135 FEB
30 135 MAR
0 135 APR
330 135 MAY
300 135 JUN
270 135 JUL
240 135 AUG
210 135 SEP
180 135 OCT
150 135 NOV
120 135 DEC
EOF
gmt pstext months.txt -R0/360/0/150 -J -F+f12p,Helvetica-Bold,blue+jCM -O -K >> $output

# Add radial scale labels
cat > radial_labels.txt << EOF
90 50 50
90 100 100
90 150 150
90 200 200
90 250 250
90 300 300
90 350 350
EOF
gmt pstext radial_labels.txt -R0/360/0/400 -J -F+f9p,Helvetica,gray50+jLM -Dj0.2i/0 -O -K >> $output

# Add axis label explanations
echo "0 -40 Monthly Snowfall (mm)" | gmt pstext -R0/360/-50/300 -J -F+f11p,Helvetica-Bold,blue+jCM -N -O -K >> $output
echo "120 -40 Potential ET (mm/month)" | gmt pstext -R -J -F+f11p,Helvetica-Bold,orange+jCM -N -O -K >> $output
echo "240 -40 Avg Daily SWE (mm)" | gmt pstext -R -J -F+f11p,Helvetica-Bold,purple+jCM -N -O >> $output

gmt psconvert -A -Tf -Z $output

echo "Radial plot complete: radial_snow.pdf"
echo "Location: ${lon}°E, ${lat}°N"

# Cleanup
rm -f snow.txt pet.txt swe.txt snow_monthly.txt monthly_totals.txt snow_polar.txt pet_polar.txt swe_polar.txt legend.txt months.txt radial_labels.txt


