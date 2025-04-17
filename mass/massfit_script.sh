
set fit logfile 'logmass_logfile.txt'
 FIT_LIMIT = 1e-6
 
fss(x,y) = sa1*x+sa2*x*x+sa3*x*x*x+sa4*y+sa5*y*y+sa6*y*y*y+sa7
fit fss(x,y) 'log_dep_single_SDSS_combi.txt' using 2:3:1:(1) via sa1,sa2,sa3,sa4,sa5,sa6,sa7

fms(x,y) = ma1*x+ma2*x*x+ma3*x*x*x+ma4*y+ma5*y*y+ma6*y*y*y+ma7
fit fms(x,y) 'log_dep_single_2MRS_combi.txt' using 2:3:1:(1) via ma1,ma2,ma3,ma4,ma5,ma6,ma7


fsml(x,y,t,u) = sc1*x+sc2*x*x+sc3*x*x*x+sc4*y+sc5*y*y+sc6*y*y*y+sc7*t+sc8*u+sc9
fit fsml(x,y,t,u) 'log_dep_multi_low_SDSS_combi.txt' using 2:3:5:6:1:(1) via sc1,sc2,sc3,sc4,sc5,sc6,sc7,sc8,sc9

fmml(x,y,t,u) = mc1*x+mc2*x*x+mc3*x*x*x+mc4*y+mc5*y*y+mc6*y*y*y+mc7*t+mc8*u+mc9
fit fmml(x,y,t,u) 'log_dep_multi_low_2MRS_combi.txt' using 2:3:5:6:1:(1) via mc1,mc2,mc3,mc4,mc5,mc6,mc7,mc8,mc9


fsmh(x,y,t,u,v) = sd1*x+sd2*x*x+sd3*x*x*x+sd4*y+sd5*y*y+sd6*y*y*y+sd7*t+sd8*u+sd9*v+sd10
fit fsmh(x,y,t,u,v) 'log_dep_multi_high_SDSS_combi.txt' using 2:3:5:6:7:1:(1) via sd1,sd2,sd3,sd4,sd5,sd6,sd7,sd8,sd9,sd10

fmmh(x,y,t,u,v) = md1*x+md2*x*x+md3*x*x*x+md4*y+md5*y*y+md6*y*y*y+md7*t+md8*u+md9*v+md10
fit fmmh(x,y,t,u,v) 'log_dep_multi_high_2MRS_combi.txt' using 2:3:5:6:7:1:(1) via md1,md2,md3,md4,md5,md6,md7,md8,md9,md10






