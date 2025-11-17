#!/bin/csh

set PS=map.ps
#set Rmap=-Rg
#set Rmap=-R-135/40/-90/90
set Rmap=-R131/142/-32/-25 # Lake_eyre final
set Rmap=-R-170/-138/54/72 # AK
#set J=-JM8/22
set J=-Js-154/90/4i/60 #sa
#set J=-Jx0.023i/0.029i

set topogrd = /Users/keyser/Documents/mirage/etopo/ETOPO1_Ice_g_gmt5.grd
#gmt grdcut ETOPO1_Ice_g_gmt.grd -Getopo_b.grd -R126/147/-40/-22

set cpt = ~/Documents/cpts/spain.cpt

# gmt gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 6.5p MAP_TICK_LENGTH_PRIMARY .5p MAP_FRAME_PEN 0.8p
gmt gmtset  MAP_FRAME_TYPE fancy+ MAP_FRAME_WIDTH 1p FONT_ANNOT_PRIMARY 8.5p MAP_TICK_LENGTH_PRIMARY 2.5p MAP_FRAME_PEN 0.8p
gmt set MAP_ANNOT_OBLIQUE 0
# gmt gmtset COLOR_BACKGROUND lightgrey COLOR_FOREGROUND lightgrey
gmt psbasemap -BneWS -Bxa8f4 -Bya4f2 -V $J $Rmap -K >! $PS
gmt pscoast  $Rmap $J -B -N1/.05p -P -Sdodgerblue3 -Glightgrey -Dh -O -W.01p -K >> $PS

#pscoast -Rg -J -B15g15 -Dc -A10000 -Glightgrey -P -O -W.01p -K >> map.ps


gmt makecpt -A50 -C$cpt -Dlightgrey/white -T0/2500 > try.cpt # for spain, it was 2000
echo "B 154/192/205" >> try.cpt # 160/185/185

# psxy stationSa.txt -Sa.52 -h0 -W0.3+cf -Cabc.cpt $J $Rmap -O -V -K >> $PS
# for -Coleron_abyss.cpt, manually removed +ve values from 'gmt makecpt -Coleron.cpt -T-7000/7000/100 > oleron_abyss.cpt'

# gmt grdimage $topogrd $J $Rmap -Bx -By -Cfes_.cpt -I+nt.85 -K -O >> $PS # -I+nt.6 original..increased for extra contarst
gmt grdimage $topogrd $J $Rmap -Bx -By -Ctry.cpt -I+nt.35 -K -O >> $PS # -I+nt.6 original..increased for extra contarst

gmt pscoast $Rmap $J -Bx -By -N -A10000 -P -K -Di -O -W.01p >> $PS #-A10+l -Ia


# awk '{print $2,$3}' ~/Research/Lake_eyre_data/station/marla.txt | gmt psxy -: -Si.15 -GDimGray $J $Rmap -O -K >> $PS ### Marla
# gmt psxy ../AK_stations.txt -W.01 -Groyalblue3 -St.15  $J $Rmap -O -V -K >> $PS
# gmt psxy ../AT_stations.txt -W.01 -Groyalblue3 -St.15  $J $Rmap -O -V -K >> $PS
# gmt psxy ../AV_stations.txt -W.01 -Groyalblue3 -St.15  $J $Rmap -O -V -K >> $PS #gray39
# gmt psxy ../CN_stations.txt -W.01 -Groyalblue3 -St.15  $J $Rmap -O -V -K >> $PS

gmt psxy ../AK_stations.txt -W.01 -Ggray39@20 -St.14  $J $Rmap -O -V -K >> $PS
gmt psxy ../AT_stations.txt -W.01 -Ggray39@20 -St.14  $J $Rmap -O -V -K >> $PS
gmt psxy ../AV_stations.txt -W.01 -Ggray39@20 -St.14  $J $Rmap -O -V -K >> $PS #gray39
gmt psxy ../CN_stations.txt -W.01 -Ggray39@20 -St.14  $J $Rmap -O -V -K >> $PS

#awk '{print $8,$7}' STA_DISTANCE_LOC_gridnumber93.txt | gmt psxy -W.24,darkmagenta -Gwhite@20 -Sc.13  $J $Rmap -O -V -K >> $PS
awk '{print $8,$7}' STA_DISTANCE_LOC_gridnumber93.txt | gmt psxy -W.01 -Gdarkmagenta -St.17  $J $Rmap -O -V -K >> $PS

# awk '{print $1-.07,$2-.15,$3}' 5g_stations_LE.txt | gmt pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS

# awk '{print $1,$2}' 6k_stations.txt | gmt psxy -W.1 -Gorangered4 -St.3  $J $Rmap -O -V -K >> $PS
# awk '{print $1-.07,$2-.15,$3}' 6k_stations.txt | gmt pstext $Rmap $J -F+f2.5p,Helvetica-Bold -Gwhite -O -P -K >> $PS
gmt gmtset FONT_ANNOT_PRIMARY 10.5p MAP_FRAME_PEN .7p FONT_LABEL 12.5p

gmt psscale -Ctry.cpt -Dx6c/6.48c+w1.5c/.12c+jTR+e -Bxa500f250+l"Elevation (m)" -By -O -K -P >> $PS

gmt psimage ../eq_AK_all_conf.png -Dx6.9c/-.68c+w3.4c/2.6c+jBR -V -O -P >> $PS


# gmt pslegend -Dx.9c/1.1c+w.80c/1.1c+o-.15c/-.15c -F+gwhite+p.15 -O $J $Rmap << EOF >> $PS
# S 0.05c t 0.15c darkmagenta - 0.08i AK
# S 0.05c t 0.15c royalblue3 - 0.08i AT
# S 0.05c t 0.15c gray34 - 0.08i AV
# S 0.05c t 0.15c red3 - 0.08i CN
# #S 0.2c s 0.3c Black - 0.2i ANSN (AU)
# # S 0.2c c 0.22c firebrick - 0.2i Eq
# EOF


gmt ps2raster -A -Tj -E920 -P -Z -Vq $PS
open map.jpg
