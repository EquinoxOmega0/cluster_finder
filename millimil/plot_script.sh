set term postscript eps enhanced color
set output 'velocities.eps'
set xrange [ -1500 : 1500 ]
set yrange [ 0 : 2000 ]
set xtics 500
set ytics 500
set xlabel " radial velocity [km/s] " font  ",20"
set ylabel " number of galaxies " font  ",20"
set mxtics 5
set mytics 5
f(x)=ncount/sqrt(2*pi)/smagg*exp(-((x-magg)**2)/(2*(smagg**2)))/20
set key font ",16"
unset key
set style line 1 lt 2 lw 5
set style line 2 lt 0 lw 5
set style line 3 lt 3 lw 5
set arrow from 232.5,0 to 232.5,2000 nohead ls 1 front
set arrow from -232.5,0 to -232.5,2000 nohead ls 1  front
set arrow from 0,0 to 0,2000 nohead ls 2  front
set arrow from 465,0 to 465,2000 nohead ls 3 front
set arrow from -465,0 to -465,2000 nohead ls 3 front
plot "vel_bin_numbers.txt" with histeps lw 5
reset












set term postscript eps enhanced color size 8in, 4in 
set output 'loglum_logmass_sdss.eps'
set multiplot layout 1,2

set tmargin at screen 0.99
set bmargin at screen 0.09
set lmargin at screen 0.09
set rmargin at screen 0.54
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
set table 'test.dat'
splot 'gal_lum_mass.txt'
unset table
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
set xtics 0.5
set ytics 0.5
set xlabel " log_1_0(L [L_{{/=12 O}&{/*-.52 O}{/=12 \267}} ])" font  ",20"
set ylabel " log_1_0(M [M_{{/=12 O}&{/*-.52 O}{/=12 \267}} ]) " font  ",20"
load 'plot_fit_lumass.txt' 
unset colorbox
unset key 
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image,f(x) lw 2



set tmargin at screen 0.99
set bmargin at screen 0.09
set lmargin at screen 0.54
set rmargin at screen 0.99
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
set table 'test.dat'
splot 'central_lum_mass.txt'
unset table
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
unset ylabel 
set xtics 0.5
unset ytics 
set xlabel " log_1_0(L [L_{{/=12 O}&{/*-.52 O}{/=12 \267}} ]) " font  ",20"
load 'plot_fit_lumass.txt' 
unset colorbox
unset key 
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, f(x) lw 2


unset multiplot
reset




set term postscript eps enhanced color size 8in, 4in 
set output 'loglum_logmass_2mass.eps'
set multiplot layout 1,2

set tmargin at screen 0.99
set bmargin at screen 0.09
set lmargin at screen 0.09
set rmargin at screen 0.54
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
set table 'test.dat'
splot 'gal_lum_mass_2mass.txt'
unset table
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
set xtics 0.5
set ytics 0.5
set xlabel " log_1_0(L [L_{{/=12 O}&{/*-.52 O}{/=12 \267}} ])" font  ",20"
set ylabel " log_1_0(M [M_{{/=12 O}&{/*-.52 O}{/=12 \267}} ]) " font  ",20"
load 'plot_fit_lumass_2mass.txt' 
unset colorbox
unset key 
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image,f(x) lw 2



set tmargin at screen 0.99
set bmargin at screen 0.09
set lmargin at screen 0.54
set rmargin at screen 0.99
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
set table 'test.dat'
splot 'central_lum_mass_2mass.txt'
unset table
set xrange [ 7.85 : 11.4 ]
 set yrange [ 10.3 : 14.3 ]
unset ylabel 
set xtics 0.5
unset ytics 
set xlabel " log_1_0(L [L_{{/=12 O}&{/*-.52 O}{/=12 \267}} ]) " font  ",20"
load 'plot_fit_lumass_2mass.txt' 
unset colorbox
unset key 
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, f(x) lw 2


unset multiplot
reset


