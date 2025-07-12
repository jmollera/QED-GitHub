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
set format x "% h" 
set format y "% .0f" 
set format x2 "% h" 
set format y2 "% h" 
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
set key inside right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key opaque
unset key
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
set samples 800, 800
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
set xtics ('-$2\pi$' -2*pi, '-$\frac{3\pi}{2}$' -3*pi/2, '-$\pi$' -pi, '-$\frac{\pi}{2}$' -pi/2, 0, '$\frac{\pi}{2}$' pi/2, '$\pi$' pi, '$\frac{3\pi}{2}$' 3*pi/2, '$2\pi$' 2*pi) norangelimit textcolor rgb "black" 
set ytics border out scale 1,0.5 nomirror norotate  autojustify
set ytics -10.0, 2.0, 10.0 norangelimit textcolor rgb "black" 
set ztics border out scale 1,0.5 nomirror norotate  autojustify
set ztics autofreq  norangelimit textcolor rgb "black" 
unset x2tics
unset y2tics
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
set xlabel '$\alpha$'
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ -2*pi : 2*pi ] noreverse nowriteback
set x2range [ 0.00000 : 50.000 ] noreverse nowriteback
set ylabel '$\tan \alpha$' 
set ylabel  font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  font "" textcolor lt -1 rotate by -270
set yrange [ -10.0 : 10.0 ] noreverse nowriteback
set y2range [ -400.0000 : 400.0000 ] noreverse nowriteback
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
set terminal cairolatex pdf size 15cm,9.5cm
set output "Ape-Trigonometria-tan.tex"
set object rectangle from -2*pi,-10.0 to 2*pi,10.0 fillstyle solid fillcolor rgb "#FFFFF0" behind
plot tan(x) linecolor rgb "red"  linewidth 6.000
#    EOF