
 #-------------------------------------------------------------------






set term postscript eps enhanced color
set output 'fp_distance_multiplicity_dependence.eps'
set multiplot layout 2,2 rowsfirst
#set size ratio 0.8

set fit logfile 'fp_distance_multiplicity.txt'

f1(x)=k1*x+d1
f2(x)=k2*x+d2
f5(x)=k5*x+d5
f10(x)=k10*x+d10

fit f1(x) 'fp_distance_variations_1.txt' using 2:3 via k1,d1
fit f2(x) 'fp_distance_variations_2.txt' using 2:3 via k2,d2
fit f5(x) 'fp_distance_variations_5.txt' using 2:3 via k5,d5
fit f10(x) 'fp_distance_variations_10.txt' using 2:3 via k10,d10


set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.11
set rmargin at screen 0.55

set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set table 'test.dat'
splot 'fp_distance_map_1.txt'
unset table
set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set xtics 50
set ytics 0.1
unset xlabel #" m_g-m_J [mag]" font  ",14"
set ylabel " log_1_0(D_F_P/D_z)" font  ",20"
set key font ",20"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      7 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2 , f1(x) lt 1 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.55
set rmargin at screen 0.99

set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set table 'test.dat'
splot 'fp_distance_map_2.txt'
unset table
set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set xtics 50
set ytics 0.1
unset xlabel # " m_g-m_H [mag]" font  ",14"
unset ylabel 
unset ytics
# set ylabel "d (m_g-m_r) + e (m_r-m_i) + f " font  ",10"
set key font ",20"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      7 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2, f2(x) lt 1 lw 2



set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.11
set rmargin at screen 0.55

set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set table 'test.dat'
splot 'fp_distance_map_5.txt'
unset table
set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set xtics 50
set ytics 0.1
set xlabel " D_z [Mpc] " font  ",20"
set ylabel " log_1_0(D_F_P/D_z)" font  ",20"
set key font ",20"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      7 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2, f5(x) lt 1 lw 2

set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.55
set rmargin at screen 0.99

set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set table 'test.dat'
splot 'fp_distance_map_10.txt'
unset table
set xrange [ 0 : 490 ]
set yrange [ -0.35 : 0.35 ]
set xtics 50
set ytics 0.1
set xlabel " D_z [Mpc] " font  ",20"
unset ylabel #" {/Symbol D} (m_H_(_o_b_s_)-m_H_(_f_i_t_)) [mag]" font  ",14"
set key font ",20"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      7 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2, f10(x) lt 1 lw 2





unset multiplot
reset

 #-------------------------------------------------------------------







# set term postscript eps enhanced color
# set output 'optimal_2MRS+.eps'
# 
# set pm3d map
# set pm3d interpolate 2,2
# splot '2/2MRS/grid_rough.txt' matrix
# reset


set term postscript eps enhanced color
set output 'completness_function.eps'
load 'bin_plot_help.txt'
set xrange [ 0 : 500 ]
set yrange [ 0 : 1.1 ]
set xtics 50
set ytics 0.25
set xlabel " comoving distance [Mpc]" font  ",18" 
set ylabel " relative completness " font  ",18"
set mxtics 5 
set mytics 5
set key font ",16"
set key at 450,0.9
set arrow from survey_limit,0 to survey_limit,1.1 nohead lw 5  lt 4 front
set arrow from merge_limit,0 to merge_limit,1.1 nohead lw 5 lt 5 front
plot "vol_cor2_plot_sorted.txt" using 1:2 title '2MRS completness' with lines lw 5 lt 0, "vol_cor_plot_sorted.txt" using 1:4 title 'SDSS completness' with lines  lw 5 lt 1, "vol_cor_plot_sorted.txt" using 1:2 title 'SDSS saturation completness' with lines  lw 5 lt 2, "vol_cor_plot_sorted.txt" using 1:3 title 'SDSS Malmquist bias completness' with lines  lw 5 lt 3
reset
 #-------------------------------------------------------------------






set term postscript size 6.5,3.5 eps enhanced color
# set term postscript eps enhanced color
set output 'map_fi_first.eps'


set multiplot layout 2,3 rowsfirst
set size ratio 0.8


set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.09
set rmargin at screen 0.39



load 'fi_fit_para_1.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_first_1.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset colorbox
unset xlabel
unset key
unset xtics
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_first_1*x*x+a_first_2*x+a_first_3 lt 7  lw 2




set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.39
set rmargin at screen 0.69

load 'fi_fit_para_2.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_first_2.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset ytics
unset xtics
unset colorbox
unset ylabel 
unset xlabel
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_first_1*x*x+a_first_2*x+a_first_3 lt 7  lw 2



set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.69
set rmargin at screen 0.99

load 'fi_fit_para_3.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_first_3.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset ytics
unset xtics
unset colorbox
unset ylabel 
unset xlabel
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_first_1*x*x+a_first_2*x+a_first_3 lt 7  lw 2



set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.09
set rmargin at screen 0.39


load 'fi_fit_para_4.txt'
set xtics 1
set ytics 1
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_first_4.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset colorbox
unset key
set xtics 1
set ytics 1
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_first_1*x*x+a_first_2*x+a_first_3 lt 7  lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.39
set rmargin at screen 0.69



load 'fi_fit_para_5.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_first_5.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_first_1*x*x+a_first_2*x+a_first_3 lt 7  lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.69
set rmargin at screen 0.99



load 'fi_fit_para_6.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_first_6.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset colorbox
unset ytics
unset key
unset ylabel 
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_first_1*x*x+a_first_2*x+a_first_3 lt 7  lw 2



unset multiplot
reset

 #-------------------------------------------------------------------



set term postscript size 6.5,3.5 eps enhanced color
# set term postscript eps enhanced color
set output 'map_fi_final.eps'





set multiplot layout 2,3 rowsfirst
set size ratio 0.8


set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.09
set rmargin at screen 0.39



load 'fi_fit_para_1.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_final_1.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset colorbox
unset xlabel
unset key
unset xtics
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2




set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.39
set rmargin at screen 0.69

load 'fi_fit_para_2.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_final_2.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set xtics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset ytics
unset xtics
unset colorbox
unset ylabel 
unset xlabel
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2



set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.69
set rmargin at screen 0.99

load 'fi_fit_para_3.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_final_3.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set xtics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset ytics
unset xtics
unset colorbox
unset ylabel 
unset xlabel
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2



set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.09
set rmargin at screen 0.39


load 'fi_fit_para_4.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_final_4.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set ytics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.39
set rmargin at screen 0.69



load 'fi_fit_para_5.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_final_5.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set xtics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.69
set rmargin at screen 0.99



load 'fi_fit_para_6.txt'

set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set table 'test.dat'
splot 'map_fi_final_6.txt'
unset table
set xrange [ 10.5 : 15 ]
set yrange [ 10.5 : 15 ]
set xtics 1
set xtics 1
set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set key font ",20"
unset colorbox
unset ytics
unset key
unset ylabel 
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2



unset multiplot
reset

 #-------------------------------------------------------------------
# 
# set term postscript eps enhanced color
# set output 'map_fi_final.eps'
# load 'fi_fit_para.txt'
# 
# set xrange [ 10 : 15 ]
# set yrange [ 10 : 15 ]
# set table 'test.dat'
# splot 'map_fi_final.txt'
# unset table
# set xrange [ 10 : 15 ]
# set yrange [ 10 : 15 ]
# set xtics 1
# set xtics 1
# set xlabel "log_1_0({/Symbol S}_i M_h_a_l_o_,_i [M_{{/=12 O}&{/*-.76 O}{/=12 \267}} ])" font  ",14"
# set ylabel "log_1_0(M_f_i [M_{{/=12 O}&{/*-.76 O}{/=12 \267}} ])" font  ",14"
# set key font ",14"
# unset colorbox
# unset key
# # set palette defined ( 0 '#FFFFFF',\
# #                       1 '#000fff',\
# #                       2 '#0090ff',\
# #                       3 '#0fffee',\
# #                       5 '#90ff70',\
# #                       10 '#ffee00',\
# #                       12 '#ff7000',\
# #                       14 '#ee0000',\
# #                       16 '#7f0000')
# set palette defined ( 0 '#FFFFFF',\
#                       1 '#BBBBCC',\
#                       3 '#5555AA',\
#                       10 '#222288')
# p 'test.dat' with image, a_final_1*x*x+a_final_2*x+a_final_3 lt 7  lw 2
# reset

 #-------------------------------------------------------------------

# # 
# set term postscript eps enhanced color
# set output 'optimal_2MRS2.eps'
# 
# 
# set xrange [ 0.1 : 1 ]
# set yrange [ 0.1 : 1 ]
# 
# 
# # set isosample 500, 500
# 
# set table 'test0.dat'
# set dgrid3d 100,100,50
# splot 'grid_rough_2MRS.txt'  #,'2/2MRS/grid_fine.txt'
# unset table
# 
# set xrange [ 0.1 : 1 ]
# set yrange [ 0.1 : 1 ]
# set xtics 0.1
# set ytics 0.1
# set xlabel "{/Symbol a}_o_p_t" font  ",14"
# set ylabel "R_o_p_t" font  ",14"
# set key font ",14"
# unset colorbox
# unset key
# set pointsize 3
# set palette defined ( 0 '#0000ff',\
#                       1 '#0090ff',\
#                       2 '#0fffee',\
#                       3 '#90ff70',\
#                       5 '#ffee00',\
#                       9 '#ee0000',\
#                       10 '#4f0000')
#                
# 
# p 'test0.dat' with image #, '2/2MRS/bestfit.txt' with points lt 2 pt 2 
# 
# reset



set term postscript eps enhanced color
set output 'optimal_2MRS.eps'


set xrange [ 0.1 : 1 ]
set yrange [ 0.1 : 1 ]


# set isosample 500, 500

set table 'test0.dat'
set dgrid3d 100,100,50
splot 'grid_combo_2MRS.txt'  #,'2/2MRS/grid_fine.txt'
unset table

set xrange [ 0.1 : 1 ]
set yrange [ 0.1 : 1 ]
set xtics 0.1
set ytics 0.1
set xlabel "{/Symbol a}_o_p_t" font  ",20"
set ylabel "R_o_p_t" font  ",20"
set key font ",20"
unset colorbox
unset key
set pointsize 3
set palette defined ( 0 '#0000ff',\
                      1 '#0090ff',\
                      2 '#0fffee',\
                      3 '#90ff70',\
                      5 '#ffee00',\
                      9 '#ee0000',\
                      10 '#4f0000')
               

p 'test0.dat' with image, 'quality_2MRS/bestfit_plot.txt' with points lt 2 pt 2 

reset



set term postscript eps enhanced color
set output 'optimal_SDSS.eps'


set xrange [ 0.1 : 1 ]
set yrange [ 0.1 : 1 ]


# set isosample 500, 500

set table 'test0.dat'
set dgrid3d 100,100,50
splot 'grid_rough_SDSS.txt'  #,'2/2MRS/grid_fine.txt'
unset table

set xrange [ 0.1 : 1 ]
set yrange [ 0.1 : 1 ]
set xtics 0.1
set ytics 0.1
set xlabel "{/Symbol a}_o_p_t" font  ",20"
set ylabel "R_o_p_t" font  ",20"
set key font ",20"
unset colorbox
unset key
set pointsize 3
set palette defined ( 0 '#0000ff',\
                      1 '#0090ff',\
                      2 '#0fffee',\
                      3 '#90ff70',\
                      5 '#ffee00',\
                      9 '#ee0000',\
                      10 '#4f0000')
               

p 'test0.dat' with image, 'quality_SDSS/bestfit_plot.txt' with points lt 2 pt 2 

reset

 #-------------------------------------------------------------------


set term postscript eps enhanced color
set output 'density_function.eps'
load 'bin_plot_help.txt'
set xrange [ 0 : 500 ]
set yrange [ 0 : 1.3 ]
set xtics 50
set ytics 0.25
set xlabel " comoving distance [Mpc]" font  ",18" 
set ylabel " fraction of the total matter density " font  ",18"
set mxtics 5 
set mytics 5
set key font ",16"
set key center top
set arrow from survey_limit,0 to survey_limit,1.3 nohead lw 5  lt 4 front
set arrow from merge_limit,0 to merge_limit,1.3 nohead lw 5 lt 5 front
plot "bin_dens_combi.txt" title 'combined dataset' with histeps lw 5, "bin_dens_2mrs.txt" title '2MRS catalogue' with histeps lw 5, "bin_dens_sdss.txt" title 'SDSS catalogue' with histeps lw 5, 0.485 title 'average value of simulated FoF groups' lw 5 lt 7

reset
 #-------------------------------------------------------------------

# vorrübergehend deaktiviert 
# set term postscript eps enhanced color
# set output 'optimal_SDSS.eps'
# 
# 
# set xrange [ 0.1 : 1 ]
# set yrange [ 0.1 : 1 ]
# 
# 
# 
# set table 'test0.dat'
# set dgrid3d 100,100,50
# splot '2/SDSS/grid.txt'  #,'2/2MRS/grid_fine.txt'
# unset table
# 
# set xrange [ 0.1 : 1 ]
# set yrange [ 0.1 : 1 ]
# set xtics 0.1
# set ytics 0.1
# set xlabel "{/Symbol a}_o_p_t" font  ",14"
# set ylabel "R_o_p_t" font  ",14"
# set key font ",14"
# unset colorbox
# unset key
# set pointsize 3
# set palette defined ( 0 '#0000ff',\
#                       1 '#0090ff',\
#                       2 '#0fffee',\
#                       3 '#90ff70',\
#                       5 '#ffee00',\
#                       9 '#ee0000',\
#                       10 '#4f0000')
#                
# 
# p 'test0.dat' with image, '2/SDSS/bestfit.txt' with points lt 2 pt 2 
# 
# reset
# 
#  #-------------------------------------------------------------------
# 



set term postscript eps enhanced color
set output 'colour_deviation.eps'
set multiplot layout 2,1

set xrange [ -0.5 : 2 ]
set yrange [ -0.4 : 0.4 ]
set table 'test.dat'
splot 'colourshift_gr_map.txt'
unset table
set xrange [ -0.5 : 2 ]
set yrange [ -0.4 : 0.4 ]
set xtics 0.5
set ytics 0.1
set xlabel "corrected colour (g-r) [mag]" font  ",16"
set ylabel "colour deviation in (g-r) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



set xrange [ -0.5 : 2 ]
set yrange [ -0.2 : 0.2 ]
set table 'test.dat'
splot 'colourshift_JK_map.txt'
unset table
set xrange [ -0.5 : 2 ]
set yrange [ -0.2 : 0.2 ]
set xtics 0.5
set ytics 0.05
set xlabel "corrected colour (J-K_s) [mag]" font  ",16"
set ylabel "colour deviation in (J-K_s) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


unset multiplot
reset

 #-------------------------------------------------------------------






set term postscript eps enhanced color
set output '2mass_sdss_transformation_fit.eps'
set multiplot layout 2,3 rowsfirst
set size ratio 0.8


set tmargin at screen 0.99
set bmargin at screen 0.54
set lmargin at screen 0.09
set rmargin at screen 0.39

set xrange [ 1.55 : 4.95 ]
set yrange [ 1.5 : 4.5 ]
set table 'test.dat'
splot 'edgeon_J.txt'
unset table
set xrange [ 1.55 : 4.95 ]
set yrange [ 1.5 : 4.5 ]
set xtics 0.5
set ytics 0.5
set xlabel " m_g-m_J [mag]" font  ",16"
set ylabel "d (m_g-m_r) + e (m_r-m_i) + f [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, x lt 7 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.54
set lmargin at screen 0.39
set rmargin at screen 0.69

set xrange [ 1.55 : 4.95 ]
set yrange [ 1.5 : 4.5 ]
set table 'test.dat'
splot 'edgeon_H.txt'
unset table
set xrange [ 1.55 : 4.95 ]
set yrange [ 1.5 : 4.5 ]
set xtics 0.5
set ytics 0.5
set xlabel " m_g-m_H [mag]" font  ",16"
unset ylabel 
unset ytics
# set ylabel "d (m_g-m_r) + e (m_r-m_i) + f " font  ",10"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, x lt 7 lw 2


set tmargin at screen 0.99
set bmargin at screen 0.54
set lmargin at screen 0.69
set rmargin at screen 0.99


set xrange [ 1.55 : 4.95 ]
set yrange [ 1.5 : 4.5 ]
set table 'test.dat'
splot 'edgeon_Ks.txt'
unset table
set xrange [ 1.55 : 4.95 ]
set yrange [ 1.5 : 4.5 ]
set xtics 0.5
set ytics 0.5
set xlabel " m_g-m_K_s [mag]" font  ",16"
unset ytics
# set ylabel "d (m_g-m_r) + e (m_r-m_i) + f " font  ",14"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, x lt 7 lw 2


set tmargin at screen 0.50
set bmargin at screen 0.09
set lmargin at screen 0.09
set rmargin at screen 0.39

set xrange [ 8.1 : 12.9 ]
set yrange [ -1 : 1 ]
set table 'test.dat'
splot 'delta_J_real_map.txt'
unset table
set xrange [ 8.1 : 12.9 ]
set yrange [ -1 : 1 ]
set xtics 1
set ytics 0.5
set xlabel " m_J [mag]" font  ",16"
set ylabel " {/Symbol D} (m_o_b_s-m_f_i_t) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.50
set bmargin at screen 0.09
set lmargin at screen 0.39
set rmargin at screen 0.69

set xrange [ 8.1 : 12.9 ]
set yrange [ -1 : 1 ]
set table 'test.dat'
splot 'delta_H_real_map.txt'
unset table
set xrange [ 8.1 : 12.9 ]
set yrange [ -1 : 1 ]
set xtics 1
set ytics 0.5
set xlabel " m_H [mag]" font  ",16"
set ylabel " {/Symbol D} (m_H_(_o_b_s_)-m_H_(_f_i_t_)) [mag]" font  ",16"
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



set tmargin at screen 0.50
set bmargin at screen 0.09
set lmargin at screen 0.69
set rmargin at screen 0.99

set xrange [ 8.1 : 12.9 ]
set yrange [ -1 : 1 ]
set table 'test.dat'
splot 'delta_Ks_real_map.txt'
unset table
set xrange [ 8.1 : 12.9 ]
set yrange [ -1 : 1 ]
set xtics 1
set ytics 0.5
set xlabel " m_K_s [mag]" font  ",16"
unset ytics

set ylabel " {/Symbol D} (m_K_s_(_o_b_s_)-m_K_s_(_f_i_t_)) [mag]" font  ",16"
set key font ",14"
unset ylabel 
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2




unset multiplot
reset




#-------------------------------------------------------------



set term postscript eps enhanced color
set output 'comparision_of_fits.eps'
set multiplot layout 2,2

set xrange [ 8.1 : 12.9 ]
set yrange [ -3 : 3 ]
set table 'test.dat'
splot 'delta_Ks_real_map.txt'
unset table
set xrange [ 8.1 : 12.9 ]
set yrange [ -3 : 3 ]
set xtics 1
set ytics 1
set xlabel " m_K_s [mag]" font  ",16"

set ylabel " {/Symbol D} (m_K_s_(_o_b_s_)-m_K_s_(_f_i_t_)) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



set xrange [ 8.1 : 12.9 ]
set yrange [ -3 : 3 ]
set table 'test.dat'
splot 'delta_Ks_bilir_map.txt'
unset table
set xrange [ 8.1 : 12.9 ]
set yrange [ -3 : 3 ]
set xtics 1
set ytics 1
set xlabel " m_K_s [mag]" font  ",16"

set ylabel " {/Symbol D} (m_K_s_(_o_b_s_)-m_K_s_(_f_i_t_)) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set xrange [ 8.1 : 12.9 ]
set yrange [ -3 : 3 ]
set table 'test.dat'
splot 'delta_Ks_mill_map.txt'
unset table
set xrange [ 8.1 : 12.9 ]
set yrange [ -3 : 3 ]
set xtics 1
set ytics 1
set xlabel " m_K_s [mag]" font  ",16"

set ylabel " {/Symbol D} (m_K_s_(_o_b_s_)-m_K_s_(_f_i_t_)) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



set xrange [ -12.1 : -25.9 ]
set yrange [ -3 : 3 ]
set table 'test.dat'
splot 'delta_Ks_bilirmill_map.txt'
unset table
set xrange [ -12.1 : -25.9 ]
set yrange [ -3 : 3 ]
set xtics 2
set ytics 1
set xlabel " M_K_s [mag]" font  ",16"

set ylabel " {/Symbol D} (m_K_s_(_o_b_s_)-m_K_s_(_f_i_t_)) [mag]" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



unset multiplot
reset










 #-------------------------------------------------------------------

# 

set term postscript eps enhanced color
set output 'distance_distribution.eps'
set multiplot layout 2,1
load 'dist_limit.txt'
set xrange [ 6 : 9 ]
set yrange [ -15 : -30 ]
set table 'test.dat'
splot 'sdss_absmag_vs_d.txt'
unset table
set xrange [ 6 : 9 ]
set yrange [ -15 : -30 ]
set xtics 0.5
set ytics 5
set xlabel " luminosity distance log_1_0(d_L [pc]) " font  ",18"
set ylabel " abs. magnitude r band [mag] " font  ",16"
set key font ",12"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
set arrow from dist_limit,-15 to dist_limit,-30 nohead ls 1 lt 4 front
f1(x)=-5.0*x+17.77+5
p 'test.dat' with image,  f1(x)  lt 1 lw 2



set xrange [ 6 : 9 ]
set yrange [ -15 : -30 ]
set table 'test.dat'
splot '2mass_absmag_vs_d.txt'
unset table
set xrange [ 6 : 9 ]
set yrange [ -15 : -30 ]
set xtics 0.5
set ytics 5
set xlabel " luminosity distance log_1_0(d_L [pc]) " font  ",18"
set ylabel " abs. magnitude K_s band [mag] " font  ",16"
set key font ",12"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
set arrow from dist_limit,-15 to dist_limit,-30 nohead ls 1 lt 4 front
f2(x)=-5.0*x+11.75+5
p 'test.dat' with image,  f2(x)  lt 1 lw 2


unset multiplot
reset




# 
#-------------------------------------------------------------------

set term postscript eps enhanced color
set output 'luminosity_function.eps'
set multiplot layout 2,1
load 'plot_luminosity_function_2MRS.txt'
load 'plot_luminosity_function_SDSS.txt'

set xrange [ -15 : -25 ]
set yrange [ -7 : -2 ]
set xtics 2
set ytics 1
set xlabel " absolute magnitude r band [mag]" font  ",18" 
set ylabel " logarithm of weighted numbers " font  ",16"
set mxtics 2 
set mytics 2
set key font ",16"
set key bottom
set key left
plot "log_luminosity_function_SDSS.txt" using 1:2 title 'corrected observational data' with histeps lw 5, "log_luminosity_function_SDSS.txt" using 1:3 title 'corrected biased mock catalogues' with histeps lw 5, "log_luminosity_function_SDSS.txt" using 1:4 title 'uncorrected unbiased mock catalogues' with histeps lw 5 #, f1(x) title 'fit on real data',f2(x) title 'fit on mock', f3(x) title 'fit on all galaxies mock' lt 7


set xrange [ -18 : -27 ]
set yrange [ -6 : -2 ]
set xtics 2
set ytics 1
set xlabel " absolute magnitude K_s band [mag]" font  ",18" 
set ylabel " logarithm of weighted numbers " font  ",16"
set mxtics 2 
set mytics 2
set key bottom
set key left
set key font ",16"
plot "log_luminosity_function_2MRS.txt" using 1:2 title 'corrected observational data' with histeps lw 5, "log_luminosity_function_2MRS.txt" using 1:3 title 'corrected biased mock catalogues' with histeps lw 5, "log_luminosity_function_2MRS.txt" using 1:4 title 'uncorrected unbiased mock catalogues' with histeps lw 5 #, f1(x) title 'fit on real data',f2(x) title 'fit on mock', f3(x) title 'fit on all galaxies mock' lt 7
# reset


unset multiplot
reset







 #-------------------------------------------------------------------




set term postscript size 7,4 eps enhanced color
set output 'all_residuals.eps'

set multiplot layout 2,3 rowsfirst
# set size ratio 0.5


set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.09
set rmargin at screen 0.39


set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/deviation_2mrs_s.txt'
unset table
set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
unset xtics
# set xtics 0.5
set ytics 0.25
unset xlabel
# set xlabel " log_1_0(M) " font  ",14" 
set ylabel " residuals {/Symbol D}_M" font  ",20"
# set mxtics 5
set mytics 5
set key font ",16"
set key bottom
set key left
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image,  0 lt 7 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.39
set rmargin at screen 0.69


set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/deviation_2mrs_ml.txt'
unset table
set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
unset ytics
unset ylabel
# set xtics 0.5
# set ytics 0.25
# set xlabel " log_1_0(M) " font  ",14" 
# set ylabel " residuals " font  ",14"
# set mxtics 5
# set mytics 5
set key font ",16"
set key bottom
set key left
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image,  0 lt 7 lw 2
# 
# set tmargin at screen 0.54
# set bmargin at screen 0.09
# set lmargin at screen 0.13
# set rmargin at screen 0.56
# 
# set xrange [ 10.7 : 15.3 ]
# set yrange [ -1.1 : 1.1 ]
# set table 'test.dat'
# splot 'mass/deviation_sdss_s.txt'
# unset table
# set xrange [ 10.7 : 15.3 ]
# set yrange [ -1.1 : 1.1 ]
# set xtics 0.5
# set ytics 0.25
# set xlabel " log_1_0(M_g_r_o_u_p [M_{{/=12 O}&{/*-.76 O}{/=12 \267}} ]) " font  ",14" 
# set ylabel " residuals {/Symbol D}_M" font  ",14"
# set mxtics 5
# set mytics 5
# set key font ",16"
# set key bottom
# set key left
# unset colorbox
# unset key
# set palette defined ( 0 '#FFFFFF',\
#                       1 '#BBBBCC',\
#                       30 '#5555AA',\
#                       100 '#222288')
# p 'test.dat' with image,  0 lt 7 lw 2
# # 
# set tmargin at screen 0.54
# set bmargin at screen 0.09
# set lmargin at screen 0.56
# set rmargin at screen 0.99
# 
# set xrange [ 10.7 : 15.3 ]
# set yrange [ -1.1 : 1.1 ]
# set table 'test.dat'
# splot 'mass/deviation_sdss_m.txt'
# unset table
# set xrange [ 10.7 : 15.3 ]
# set yrange [ -1.1 : 1.1 ]
# unset ytics
# set xtics 0.5
# # set ytics 0.25
# set xlabel " log_1_0(M_g_r_o_u_p [M_{{/=12 O}&{/*-.76 O}{/=12 \267}} ]) " font  ",14" 
# # set ylabel " residuals " font  ",14"
# unset ylabel
# set mxtics 5
# # set mytics 5
# set key font ",16"
# set key bottom
# set key left
# unset colorbox
# unset key
# set palette defined ( 0 '#FFFFFF',\
#                       1 '#BBBBCC',\
#                       30 '#5555AA',\
#                       100 '#222288')
# p 'test.dat' with image,  0 lt 7 lw 2
# 

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.69
set rmargin at screen 0.99


set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/deviation_2mrs_mh.txt'
unset table
set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]

unset xtics
#set xlabel " m_g-m_K_s [mag]" font  ",14"
unset ytics
# set ylabel "d (m_g-m_r) + e (m_r-m_i) + f " font  ",14"
set key font ",16"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.09
set rmargin at screen 0.39




set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/deviation_sdss_s.txt'
unset table
set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set xtics 1
set ytics 0.25
set xlabel " log_1_0(M_g_r_o_u_p_,_m_o_c_k [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ]) " font  ",20" 
set ylabel " residuals {/Symbol D}_M" font  ",20"
set mxtics 5
set mytics 5

set key font ",16"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.39
set rmargin at screen 0.69

set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/deviation_sdss_ml.txt'
unset table
set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set xtics 1
set ytics 0.25
set xlabel " log_1_0(M_g_r_o_u_p_,_m_o_c_k [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ]) " font  ",20" 
unset ylabel # " {/Symbol D} (m_H_(_o_b_s_)-m_H_(_f_i_t_)) [mag]" font  ",14"
set key font ",16"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.69
set rmargin at screen 0.99

set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/deviation_sdss_mh.txt'
unset table
set xrange [ 10.7 : 15.3 ]
set yrange [ -1.1 : 1.1 ]
set xtics 1
set ytics 0.25
set xlabel " log_1_0(M_g_r_o_u_p_,_m_o_c_k [M_{{/=12 O}&{/*-.55 O}{/=12 \267}} ]) " font  ",20" 
unset ytics

#set ylabel " {/Symbol D} (m_K_s_(_o_b_s_)-m_K_s_(_f_i_t_)) [mag]" font  ",14"
set key font ",16"
unset ylabel 
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2







unset multiplot
reset






# 
# #-------------------------------------------------------------
# 
# 
# 



set term postscript size 5,3.5 eps enhanced color
set output 'res_single.eps'
set multiplot layout 2,2 rowsfirst
#set size ratio 0.8



set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.11
set rmargin at screen 0.55


set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_lum_2mrs_s.txt'
unset table
set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset xlabel 
set ylabel "residuals {/Symbol D}_M" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.55
set rmargin at screen 0.99

set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_dist_2mrs_s.txt'
unset table
set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2



set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.11
set rmargin at screen 0.55


set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_lum_sdss_s.txt'
unset table
set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel " log_1_0(L_t_o_t [L_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel " residuals {/Symbol D}_M" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.55
set rmargin at screen 0.99

set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_dist_sdss_s.txt'
unset table
set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.25
set ytics 0.5
set xlabel "log_1_0(D_L [Mpc])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2





unset multiplot
reset

 #-------------------------------------------------------------------


 
set term postscript size 10,4 eps enhanced color
set output 'res_multi_low.eps'
set multiplot layout 2,4 rowsfirst
#set size ratio 0.8




set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.09
set rmargin at screen 0.31

set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_lum_2mrs_ml.txt'
unset table
set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset xlabel 
set ylabel "residuals {/Symbol D}_M" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.31
set rmargin at screen 0.53



set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_dist_2mrs_ml.txt'
unset table
set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2
set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.53
set rmargin at screen 0.75

set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_sigma_2mrs_ml.txt'
unset table
set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.75
set rmargin at screen 0.97


set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_rad_2mrs_ml.txt'
unset table
set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2










set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.09
set rmargin at screen 0.31




set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_lum_sdss_ml.txt'
unset table
set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel " log_1_0(L_t_o_t [L_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel " residuals {/Symbol D}_M" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.31
set rmargin at screen 0.53



set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_dist_sdss_ml.txt'
unset table
set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.25
set ytics 0.5
set xlabel "log_1_0(D_L [Mpc])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.53
set rmargin at screen 0.75



set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_sigma_sdss_ml.txt'
unset table
set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel "log_1_0({/Symbol s}_g_r_o_u_p [km/s])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.75
set rmargin at screen 0.97

set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_rad_sdss_ml.txt'
unset table
set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel "log_1_0(R_g_r_o_u_p [kpc])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

unset multiplot
reset






# 
# #-------------------------------------------------------------
# 
# 
# 



set term postscript size 12,5 eps enhanced color
set output 'res_multi_high.eps'
set multiplot layout 2,5 rowsfirst
#set size ratio 0.8



set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.09
set rmargin at screen 0.27


set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_lum_2mrs_mh.txt'
unset table
set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset xlabel 
set ylabel "residuals {/Symbol D}_M" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.27
set rmargin at screen 0.45

set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_dist_2mrs_mh.txt'
unset table
set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.45
set rmargin at screen 0.63

set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_sigma_2mrs_mh.txt'
unset table
set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.63
set rmargin at screen 0.81

set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_rad_2mrs_mh.txt'
unset table
set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2




set tmargin at screen 0.99
set bmargin at screen 0.55
set lmargin at screen 0.81
set rmargin at screen 0.99

set xrange [ 0.6 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_nvis_2mrs_mh.txt'
unset table
set xrange [ 0.6 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
unset xlabel 
unset xtics
unset ylabel 
unset ytics
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2












set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.09
set rmargin at screen 0.27


set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_lum_sdss_mh.txt'
unset table
set xrange [ 9 : 12.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel " log_1_0(L_t_o_t [L_{{/=12 O}&{/*-.55 O}{/=12 \267}} ])" font  ",20"
set ylabel " residuals {/Symbol D}_M" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.27
set rmargin at screen 0.45

set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_dist_sdss_mh.txt'
unset table
set xrange [ 0.9 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.25
set ytics 0.5
set xlabel "log_1_0(D_L [Mpc])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.45
set rmargin at screen 0.63

set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_sigma_sdss_mh.txt'
unset table
set xrange [ 1.4 : 3.4 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel "log_1_0({/Symbol s}_g_r_o_u_p [km/s])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2


set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.63
set rmargin at screen 0.81

set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_rad_sdss_mh.txt'
unset table
set xrange [ 0.1 : 3.9 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel "log_1_0(R_g_r_o_u_p [kpc])" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2

set tmargin at screen 0.55
set bmargin at screen 0.11
set lmargin at screen 0.81
set rmargin at screen 0.99

set xrange [ 0.6 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set table 'test.dat'
splot 'mass/res_nvis_sdss_mh.txt'
unset table
set xrange [ 0.6 : 2.8 ]
set yrange [ -1.1 : 1.1 ]
set xtics 0.5
set ytics 0.5
set xlabel "log_1_0(N_F_O_F)" font  ",20"
unset ylabel
set key font ",14"
unset ytics
unset colorbox
unset ylabel 
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      30 '#5555AA',\
                      100 '#222288')
p 'test.dat' with image, 0 lt 7 lw 2




unset multiplot
reset






#-------------------------------------------------------------

set term postscript eps enhanced color
set output 'bin_clusters_2MRS.eps'

set xrange [ 0.5 : 100.5 ]
set yrange [ 0.5 : 35000 ]
set xtics 10
set logscale y
# set ytics 1
set ylabel " number of groups" font  ",20" 
set xlabel " N_F_o_F " font  ",20"
set mxtics 2 
set mytics 5
set key bottom
set key left
set key font ",16"
unset key
plot "catalogues/binclusterfof_2MRS.txt" using 1:2 with histeps lw 5
reset






#-------------------------------------------------------------

set term postscript eps enhanced color
set output 'bin_clusters_SDSS.eps'

set xrange [ 0.5 : 100.5 ]
set yrange [ 0.5 : 250000 ]
set xtics 10
set logscale y
# set ytics 1
set ylabel " number of groups" font  ",20" 
set xlabel " N_F_o_F " font  ",20"
set mxtics 2 
set mytics 5
set key bottom
set key left
set key font ",16"
unset key
plot "catalogues/binclusterfof_SDSS.txt" using 1:2 with histeps lw 5
reset




#-------------------------------------------------------------

set term postscript eps enhanced color
set output 'bin_clusters_combo.eps'

set xrange [ 0.5 : 100.5 ]
set yrange [ 0.5 : 250000 ]
set xtics 10
set logscale y
# set ytics 1
set ylabel " number of groups" font  ",20" 
set xlabel " N_F_o_F " font  ",20"
set mxtics 2 
set mytics 5
set key bottom
set key left
set key font ",16"
unset key
plot "catalogues/binclusterfof_SDSS.txt" using 1:2 with histeps lw 5, "catalogues/binclusterfof_2MRS.txt" using 1:2 with histeps lw 5
reset



#-------------------------------------------------------------
set term postscript size 6,2 eps enhanced color
set output 'bin_ndes.eps'
set multiplot layout 1,2 rowsfirst

# set term postscript eps enhanced color
# set output 'bin_ndes.eps'

set xrange [ 0 : 0.11 ]
set yrange [ 0 : 0.05 ]
set xtics 0.02
# set logscale y
set ytics 0.01
set ylabel " galaxies / Mpc^3" font  ",16" 
set xlabel " redshift " font  ",16"
set mxtics 2 
set mytics 2
# set key bottom
# set key left
set key font ",16"
# unset key
plot   "bin_ndens_2MRS.txt" using 1:3 with histeps lw 3 lt 2 title 'mock catalogues', "bin_ndens_2MRS.txt" using 1:4 with histeps lw 3 lt 2 notitle, "bin_ndens_2MRS.txt" using 1:5 with histeps lw 3 lt 2 notitle, "bin_ndens_2MRS.txt" using 1:6 with histeps lw 3 lt 2 notitle, "bin_ndens_2MRS.txt" using 1:7 with histeps lw 3 lt 2 notitle, "bin_ndens_2MRS.txt" using 1:8 with histeps lw 3 lt 2 notitle, "bin_ndens_2MRS.txt" using 1:9 with histeps lw 3 lt 2 notitle, "bin_ndens_2MRS.txt" using 1:10 with histeps lw 3 lt 2  notitle, "bin_ndens_2MRS.txt" using 1:2 with histeps lw 5 lt 1 title 'observational data'







set xrange [ 0 : 0.11 ]
set yrange [ 0 : 0.1 ]
set xtics 0.02
# set logscale y
set ytics 0.02
set ylabel " galaxies / Mpc^3" font  ",16" 
set xlabel " redshift " font  ",16"
set mxtics 2 
set mytics 2
# set key bottom
# set key left
set key font ",16"
# unset key
plot   "bin_ndens_SDSS.txt" using 1:3 with histeps lw 3 lt 2 title 'mock catalogues', "bin_ndens_SDSS.txt" using 1:4 with histeps lw 3 lt 2 notitle, "bin_ndens_SDSS.txt" using 1:5 with histeps lw 3 lt 2 notitle, "bin_ndens_SDSS.txt" using 1:6 with histeps lw 3 lt 2 notitle, "bin_ndens_SDSS.txt" using 1:7 with histeps lw 3 lt 2 notitle, "bin_ndens_SDSS.txt" using 1:8 with histeps lw 3 lt 2 notitle, "bin_ndens_SDSS.txt" using 1:9 with histeps lw 3 lt 2 notitle, "bin_ndens_SDSS.txt" using 1:10 with histeps lw 3 lt 2  notitle, "bin_ndens_SDSS.txt" using 1:2 with histeps lw 5 lt 1 title 'observational data'



unset multiplot

reset




#-------------------------------------------------------------

set term postscript eps enhanced color
set output 'bin_earlytypes_percluster.eps'

set xrange [ 0.5 : 100.5 ]
set yrange [ 0.5 : 50000 ]
set xtics 10
set logscale y
# set ytics 1
set ylabel " number of groups" font  ",20" 
set xlabel " N_E_T_G / group" font  ",20"
set mxtics 2 
set mytics 5
set key bottom
set key left
set key font ",16"
unset key
plot "catalogues/bin_earlytypes_percluster.txt" using 1:2 with histeps lw 5
reset



#-------------------------------------------------------------



set term postscript eps enhanced color
set output 'fp_distance_err.eps'
# set multiplot layout 1,3

set xrange [ 0.5 : 30.5 ]
set yrange [ -0.35 : 0.35 ]
set table 'test.dat'
splot 'fp_distance_err_m.txt'
unset table
set xrange [ 0.5 : 30.5 ]
set yrange [ -0.35 : 0.35 ]
set xtics 5
set ytics 0.1
set xlabel  " N_E_T_G / group" font  ",20"
set ylabel " log_1_0(D_F_P/D_z)" font  ",20"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      7 '#5555AA',\
                      10 '#222288')
# set palette defined ( 0 '#ffffff',\
#                       1 '#0000ff',\
#                       2 '#0090ff',\
#                       3 '#0fffee',\
#                       4 '#90ff70',\
#                       5 '#ffee00',\
#                       7 '#ff7000',\
#                       9 '#ee0000',\
#                       10 '#7f0000')
p 'test.dat' with image, 0 with lines lt 4 lw 1, 'fp_distance_statistics_m.txt' using 1:2  with lines lt 7 lw 2, 'fp_distance_statistics_m.txt' using 1:4  with lines lt 1 lw 2,'fp_distance_statistics_m.txt' using 1:5  with lines lt 1 lw 2, 0.0920/sqrt(x) lt 2 lw 2, -0.0920/sqrt(x) lt 2 lw 2  

# 
# set xrange [ 0 : 30 ]
# set yrange [ -100 : 100 ]
# set table 'test.dat'
# splot 'fp_distance_err_s.txt'
# unset table
# set xrange [ 0 : 30 ]
# set yrange [ -100 : 100 ]
# set xtics 10
# set ytics 0.2
# set xlabel  " N_E_T_G / cluster" font  ",14"
# set ylabel " relative distance error [%]" font  ",14"
# set key font ",14"
# unset colorbox
# unset key
# set palette defined ( 0 '#ffffff',\
#                       1 '#0000ff',\
#                       2 '#0090ff',\
#                       3 '#0fffee',\
#                       5 '#90ff70',\
#                       10 '#ffee00',\
#                       30 '#ff7000',\
#                       75 '#ee0000',\
#                       100 '#7f0000')
# p 'test.dat' with image, 0 lt 7
# 
# 
# set xrange [ 0 : 30 ]
# set yrange [ -100 : 100 ]
# set table 'test.dat'
# splot 'fp_distance_err_p.txt'
# unset table
# set xrange [ 0 : 30 ]
# set yrange [ -100 : 100 ]
# set xtics 10
# set ytics 0.2
# set xlabel  " N_E_T_G / cluster" font  ",14"
# set ylabel " relative distance error [%]" font  ",14"
# set key font ",14"
# unset colorbox
# unset key
# set palette defined ( 0 '#ffffff',\
#                       1 '#0000ff',\
#                       2 '#0090ff',\
#                       3 '#0fffee',\
#                       5 '#90ff70',\
#                       10 '#ffee00',\
#                       30 '#ff7000',\
#                       40 '#ee0000',\
#                       50 '#7f0000')
# p 'test.dat' with image, 0 lt 7
# 
# 
# 
# unset multiplot
reset










set term postscript size 8.5,3 eps enhanced color
# set term postscript eps enhanced color
set output 'resdist_SDSS.eps'



set multiplot layout 2,4 rowsfirst
set size ratio 0.65


set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13 ]
set table 'test.dat'
splot 'distdep_lum_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(L_t_o_t [L_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13 ]
set table 'test.dat'
splot 'distdep_lumobs_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(L_o_b_s [L_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 1.4 : 3 ]
set table 'test.dat'
splot 'distdep_sigma_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 1.4 : 3 ]
set xtics 100
set ytics 0.5
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0({/Symbol s}_g_r_o_u_p [km/s])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 0 : 3.5 ]
set table 'test.dat'
splot 'distdep_rad_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 0 : 3.5 ]
set xtics 100
set ytics 0.5
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(R_g_r_o_u_p [kpc])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image




set xrange [ 0 : 500 ]
set yrange [ 10.9 : 15 ]
set table 'test.dat'
splot 'distdep_mass_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 10.9 : 15 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(M_t_o_t [M_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 8 : 14 ]
set table 'test.dat'
splot 'distdep_mstar_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 8 : 14 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(M_* [M_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image







set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13]
set table 'test.dat'
splot 'distdep_mdyn_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(M_d_y_n [M_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image





set xrange [ 0 : 500 ]
set yrange [ -0.1 : 2 ]
set table 'test.dat'
splot 'distdep_nvis_SDSS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ -0.1 : 2 ]
set xtics 100
set ytics 0.5
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(N_F_O_F)" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image


unset multiplot
reset

 #-------------------------------------------------------------------









set term postscript size 8.5,3 eps enhanced color
# set term postscript eps enhanced color
set output 'resdist_2MRS.eps'



set multiplot layout 2,4 rowsfirst
set size ratio 0.65


set xrange [ 0 : 500 ]
set yrange [ 8.5 : 15.5 ]
set table 'test.dat'
splot 'distdep_lum_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 8.5 : 15.5 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(L_t_o_t [L_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image

set xrange [ 0 : 500 ]
set yrange [ 8 : 13.5 ]
set table 'test.dat'
splot 'distdep_lumobs_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 8 : 13.5 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(L_o_b_s [L_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image

set xrange [ 0 : 500 ]
set yrange [ 1.4 : 3 ]
set table 'test.dat'
splot 'distdep_sigma_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 1.4 : 3 ]
set xtics 100
set ytics 0.5
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0({/Symbol s}_g_r_o_u_p [km/s])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 0 : 3.5 ]
set table 'test.dat'
splot 'distdep_rad_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 0 : 3.5 ]
set xtics 100
set ytics 0.5
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(R_g_r_o_u_p [kpc])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 10.5 : 15.5 ]
set table 'test.dat'
splot 'distdep_mass_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 10.5 : 15.5 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(M_t_o_t [M_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image



set xrange [ 0 : 500 ]
set yrange [ 8.5 : 14 ]
set table 'test.dat'
splot 'distdep_mstar_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 8.5 : 14 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(M_* [M_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image





set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13]
set table 'test.dat'
splot 'distdep_mdyn_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ 7.5 : 13 ]
set xtics 100
set ytics 1
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(M_d_y_n [M_{{/=12 O}&{/*-.65 O}{/=12 \267}} ])" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image





set xrange [ 0 : 500 ]
set yrange [ -0.1 : 2 ]
set table 'test.dat'
splot 'distdep_nvis_2MRS.txt'
unset table
set xrange [ 0 : 500 ]
set yrange [ -0.1 : 2 ]
set xtics 100
set ytics 0.5
set xlabel "D_L [Mpc]" font  ",16"
set ylabel "log_1_0(N_F_O_F)" font  ",16"
set key font ",14"
unset colorbox
unset key
set palette defined ( 0 '#FFFFFF',\
                      1 '#BBBBCC',\
                      3 '#5555AA',\
                      10 '#222288')
p 'test.dat' with image


unset multiplot
reset

 #-------------------------------------------------------------------






































