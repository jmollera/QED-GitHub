#!/gnuplot
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 7    last modified 2019-05-29 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2018
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal windows 0 color solid butt enhanced
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border front linecolor rgb "blue"  linewidth 3.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0.00000, 0.00000 
set style ellipse size graph 0.05, 0.03, first 0.00000 angle 0 units xy
set dummy x, y
set format x "% .0f" 
set format y "% .0f" 
set format x2 "% h" 
set format y2 "% .0f" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set grid back nopolar xtics nomxtics ytics nomytics noztics nomztics
set grid nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid linecolor "blue" linetype 1 linewidth 1.00 
set raxis
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title ""
set key inside left bottom vertical Left reverse enhanced box
#set key at 1305, 30  Left reverse  box
set key noinvert samplen 1 spacing 1 width -2 height 1 
set key maxcolumns 3 maxrows 0
set key opaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
set decimalsign ','
set view 60, 30, 1, 1
set samples 2000
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set ticslevel 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set xtics border out scale 1,0.5 nomirror norotate  autojustify
set xtics 1200,50,1500.0 norangelimit textcolor rgb "blackz" 
set ytics border out scale 1,0.5 nomirror norotate  autojustify
set ytics 0.0,10.0,90.0 norangelimit textcolor rgb "black" 
set ztics border out scale 1,0.5 nomirror norotate  autojustify
set ztics autofreq  norangelimit textcolor rgb "black" 
unset x2tics
set y2tics border out scale 1,0.5 nomirror norotate  autojustify
set y2tics 0.0,5.0,45.0 norangelimit textcolor rgb "black"
set cbtics border out scale 1,0.5 nomirror norotate  autojustify
set cbtics autofreq  norangelimit textcolor rgb "black" 
set rtics axis out scale 1,0.5 nomirror norotate  autojustify
set rtics autofreq  norangelimit textcolor rgb "black" 
set paxis 1 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 1 tics autofreq  rangelimit textcolor rgb "black" 
set paxis 2 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 2 tics autofreq  rangelimit textcolor rgb "black" 
set paxis 3 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 3 tics autofreq  rangelimit textcolor rgb "black" 
set paxis 4 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 4 tics autofreq  rangelimit textcolor rgb "black" 
set paxis 5 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 5 tics autofreq  rangelimit textcolor rgb "black" 
set paxis 6 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 6 tics autofreq  rangelimit textcolor rgb "black" 
set paxis 7 tics border out scale 1,0.5 nomirror norotate  autojustify
set paxis 7 tics autofreq  rangelimit textcolor rgb "black" 
set title "" 
set title  font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel '$n\ped{m}\, / \,\unit{r/min}$'
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 1200.0000 : 1500.0000 ] noreverse nowriteback
set x2range [ 0.0000 : 5.0000 ] noreverse nowriteback
set ylabel  '$T\ped{m}\, / \,\unit{N.m}, \quad T\ped{load}\, / \,\unit{N.m}$'
set ylabel  font "" textcolor lt -1 rotate by -270
set y2label  '$I_1\, / \,\unit{A}$'
set y2label  font "" textcolor lt -1 rotate by -270
set yrange [ 0.000 : 90.000 ] noreverse nowriteback
set y2range [ 0.000 : 45.000 ] noreverse nowriteback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set paxis 1 range [ * : * ] noreverse nowriteback
set paxis 2 range [ * : * ] noreverse nowriteback
set paxis 3 range [ * : * ] noreverse nowriteback
set paxis 4 range [ * : * ] noreverse nowriteback
set paxis 5 range [ * : * ] noreverse nowriteback
set paxis 6 range [ * : * ] noreverse nowriteback
set paxis 7 range [ * : * ] noreverse nowriteback
set zero 1e-008
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "Catalan_Spain.1252"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
GNUTERM = "windows"
ARGC = 0
ARG0 = ""
set terminal cairolatex pdf size 16.5cm,10cm
set output "Cap-Motors-Induccio-MotTensRedSolEx-1-2.tex"
set object rectangle from 0,0 to 1500,150 fillstyle solid fillcolor rgb "#FFFFF0" behind
set arrow from 1431, 56.6 to 1431, 0 dt 2 linecolor rgb "dark-violet" linewidth 3.0 head noborder back
set label '{$\scriptstyle\qty{1431}{r/min}$}' at 1427, 15 center rotate by -270 front
set arrow from 1387.5, 53.4 to 1387.5, 0 dt 2 linecolor rgb "dark-violet" linewidth 3.0 head noborder back
set label '{$\scriptstyle\qty{1387,5}{r/min}$}' at 1383.7, 15 center rotate by -270 front
set arrow from 1431, 56.6 to 1200, 56.6 dt 2 linecolor rgb "dark-violet" linewidth 3.0 head noborder back
set label '{$\scriptstyle\qty{56,6}{N.m}$}' at 1230, 58 center front
set arrow from 1387.5, 53.4 to 1200, 53.4 dt 2 linecolor rgb "dark-violet" linewidth 3.0 head noborder back
set label '{$\scriptstyle\qty{53,4}{N.m}$}' at 1230, 51.5 center front
set arrow from 1431, second 16.6 to 1500, second 16.6 dt 2 linecolor rgb "dark-violet" linewidth 3.0 head noborder back
set label '{$\scriptstyle\qty{16,6}{A}$}' at 1480, second 17.3 center front
set arrow from 1387.5, second 19.6 to 1500, second 19.6 dt 2 linecolor rgb "dark-violet" linewidth 3.0 head noborder back
set label '{$\scriptstyle\qty{19,6}{A}$}' at 1480, second 20.6 center front
plot "../Python/Auxiliars/Cap-Motors-Induccio-MotTensRedSolEx-Aux-1.txt"  using 1:2 with lines axes x1y1 linecolor rgb "black" linewidth 6.000 title '$T\ped{load}$', "../Python/Auxiliars/Cap-Motors-Induccio-MotTensRedSolEx-Aux-1.txt"   using 1:3 with lines axes x1y1 linecolor rgb "red"  linewidth 6.000 title '$T\ped{m,\qty{100}{\%}}$', "../Python/Auxiliars/Cap-Motors-Induccio-MotTensRedSolEx-Aux-1.txt"  using 1:4 with lines axes x1y1 linecolor rgb "sea-green"  linewidth 6.000 title '$T\ped{m,\qty{80}{\%}}$', "../Python/Auxiliars/Cap-Motors-Induccio-MotTensRedSolEx-Aux-1.txt"  using 1:5 with lines axes x1y2 linecolor rgb "orange"  linewidth 6.000 title '$I\ped{1,\qty{100}{\%}}$', "../Python/Auxiliars/Cap-Motors-Induccio-MotTensRedSolEx-Aux-1.txt"  using 1:6 with lines axes x1y2 linecolor rgb "green"  linewidth 6.000 title '$I\ped{1,\qty{80}{\%}}$'
#    EOF
