


set term postscript eps enhanced color size 20,10
# set 
set output 'SDSS_map.eps'
set multiplot layout 1,3

# set tmargin at screen 0.99
# set bmargin at screen 0.09
# set lmargin at screen 0.09
# set rmargin at screen 0.39

set xtics font ", 36" 
set ytics font ", 36" 

set size ratio 1
set xrange [ 0 : 520.55 ]
set yrange [ 0 : 520.55 ]
set table 'test.dat'
splot 'SDSS_1map.txt'
unset table
set xrange [ 0 : 520.55 ]
set yrange [ 0 : 520.55 ]
set xtics 100
set ytics 100
set xlabel " x [Mpc]" font  ",36"
set ylabel  " y [Mpc] " font  ",36"
# set key font ",36"
unset colorbox
set palette defined ( 0 '#FFFFFF',\
                      1 '#ff0000',\
                      2 '#ffff00',\
                      3 '#00ff00')
p 'test.dat' with image, 'SDSS_neighbours_1overlap.txt' pt 1 ps 0.5 lt 3  notitle

# set tmargin at screen 0.99
# set bmargin at screen 0.09
# set lmargin at screen 0.09
# set rmargin at screen 0.39


set size ratio 1
set xrange [ 0 : 520.55 ]
set yrange [ 0 : 520.55 ]
set table 'test.dat'
splot 'SDSS_2map.txt'
unset table
set xrange [ 0 : 520.55 ]
set yrange [ 0 : 520.55 ]
set xtics 100
set ytics 100
set xlabel " x [Mpc]" font  ",36"
set ylabel  " y [Mpc] " font  ",36"
# set key font ",36"
unset colorbox
set palette defined ( 0 '#FFFFFF',\
                      1 '#ff0000',\
                      2 '#ffff00',\
                      3 '#00ff00')
p 'test.dat' with image, 'SDSS_neighbours_2overlap.txt' pt 1 ps 0.5 lt 3  notitle

# set tmargin at screen 0.99
# set bmargin at screen 0.09
# set lmargin at screen 0.09
# set rmargin at screen 0.39


set size ratio 1
set xrange [ 0 : 520.55 ]
set yrange [ 0 : 520.55 ]
set table 'test.dat'
splot 'SDSS_3map.txt'
unset table
set xrange [ 0 : 520.55 ]
set yrange [ 0 : 520.55 ]
set xtics 100
set ytics 100
set xlabel " x [Mpc]" font  ",36"
set ylabel  " y [Mpc] " font  ",36"
# set key font ",36"
unset colorbox
set palette defined ( 0 '#FFFFFF',\
                      1 '#ff0000',\
                      2 '#ffff00',\
                      3 '#00ff00')
p 'test.dat' with image, 'SDSS_neighbours_3overlap.txt' pt 1 ps 0.5 lt 3  notitle


 unset multiplot
reset




