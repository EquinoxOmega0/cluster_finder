PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer(kind=8) :: io_err,n,i,ii,iii,n2,counter,old_fof_n,fof_n,n_dreieck,help_count,help_index,central_index
character(200) :: filename,intial_values,listname,fulllistname
integer(kind=8) :: n_l,n_fl,calibrationmode,n_truegroups_mock,n_truegroups_all,n2plus_truemock,n2plus_trueall
logical :: flipvar
integer(kind=8) :: n_g2_fof,n_g2_mock,n_g2_bij,n_multi_cut,n_help,n_used,n_fp

logical, allocatable :: group_bij(:),samecluster(:),oldact(:)
double precision, allocatable ::  P_mock(:),P_fof(:)

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor,R_factor,basic_link
double precision :: d_red,d_ang_dist,angular_sep,delta_z,hreal,l_next,hl_z,hl_a,l_old,mabs_vollim
double precision :: avlogmass,avmass,mag_sol_r,help_lum,s_power,mag_sol_g,weight,dummy_var,dummy_var2

double precision, allocatable ::  ra(:),dec(:),z(:),mag_g(:),mag_r(:)
logical, allocatable :: active(:),fof(:),fof_recursive(:),iter_help(:)
double precision, allocatable :: ra_clustercore(:),dec_clustercore(:),z_clustercore(:),lumtot_clustercore(:)
double precision, allocatable :: sigma_clustercore(:),dummy_core(:),radius_clustercore(:),angradius_clustercore(:)
double precision, allocatable :: lumdistance_clustercore(:),mass_clustercore(:),colour_clustercore(:)
double precision, allocatable :: lumotherband(:),lum_g(:),sigma_gap(:),lumobs_clustercore(:),mdyn_clustercore(:)
integer(kind=8), allocatable :: members_clustercore(:),clusterindex(:),clustergalnmember(:),members_clustercore2(:)
integer(kind=8), allocatable :: id1(:),id2(:),numberofgal(:)
double precision, dimension(1:7) :: mass_coeff_s
double precision, dimension(1:9) :: mass_coeff_ml
double precision, dimension(1:10) :: mass_coeff_mh
double precision :: loglum,logmdyn,logcolour,loglumdist,logsigma,logmember,lograd,loglumobs

double precision, allocatable ::  angular_dist(:),ra_rad(:),dec_rad(:),lum(:),volume_weight(:)
double precision, allocatable ::  mass_err(:),lumtot_err(:)
double precision, allocatable ::  lum_dist(:),vol_int(:),b_eff(:),R_eff(:)
double precision, allocatable ::  basic_link_a(:),basic_link_R(:),vred(:),lum_gal(:),lum_int(:)
double precision :: mag_limit,D_limit,t1,t2,t3,z_limit,V_limit,V_sum,binhelp,mag_vislim,magmin,magmax
double precision :: deltamag,strechdistant,av_rad_pec_vel,angdist_help,area_cover
double precision :: mock_cover,lum_power,dist_help,red_measure,mhelpg,mhelpr,lum_sum
double precision, allocatable :: ra_clustercore2(:),dec_clustercore2(:),z_clustercore2(:)
double precision, allocatable :: lumtot_clustercore2(:),lumtot_err2(:),lumobs_clustercore2(:)
double precision, allocatable :: mass_clustercore2(:),mass_err2(:),mdyn_clustercore2(:)
double precision, allocatable :: sigma_clustercore2(:),radius_clustercore2(:),angradius_clustercore2(:)
double precision, allocatable :: lumdistance_clustercore2(:)

double precision, allocatable :: fp_ra_clustercore(:),fp_dec_clustercore(:),fp_z_clustercore(:)
double precision, allocatable :: cluster_fp_adist_m(:),err_fpd_m(:),cluster_fp_cdist_m(:)
double precision, allocatable :: cluster_fp_ldist_m(:),cluster_z_adist_m(:),err_red_m(:)
double precision, allocatable :: cluster_z_cdist_m(:),cluster_z_ldist_m(:),fp_err_z_m(:)
integer(kind=8), allocatable ::  n_fp_in_cluster(:),fp_members_clustercore(:),fp_id(:)

double precision, allocatable :: sdss_gal_b(:),sdss_gal_l(:),twomass_gal_b(:),twomass_gal_l(:)
double precision, allocatable :: fp_gal_b(:),fp_gal_l(:)



double precision :: ra_help,dec_help,xhelp1,xhelp2,xhelp3
double precision :: sigma_lumtot_s,sigma_lumtot_ml,sigma_lumtot_mh,z_help
real :: v_cmb,l_cmb,b_cmb,z_cmb,gal_b,gal_l,z_cmb_x,z_cmb_y,z_cmb_z,z_gal_x,z_gal_y,z_gal_z,arcsec
real :: v_mws,l_mws,b_mws,z_mws,gal_b_mws,gal_l_mws,z_mws_x,z_mws_y,z_mws_z

double precision :: E_fof,E_mock,E_tot,Q_fof,Q_mock,Q_tot,S_tot,P_help,rms_s,rms_ml,rms_mh

double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y
integer, dimension(1:102) :: binclustermock,binclusterfof

character(20) :: appendix


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme CLUSTER FINDER started'
WRITE(*,*) '============================================================'



OPEN(50,file='input_cf.txt')

READ(50,*) Omega_m,Omega_l,H0
READ(50,*) b_factor,R_factor,lum_power
READ(50,*) appendix
 CLOSE(50)
PI=ACOS(-1.D0)

q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)


! get length of file
OPEN(50,file='catalogues/cluster_list_all_SDSS.txt')
io_err=0
n=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n=n+1
END DO
 CLOSE(50)
n=n-1

WRITE(*,*) n,'galaxies will be used'


! 
allocate(ra_clustercore(1:n))
allocate(dec_clustercore(1:n))
allocate(z_clustercore(1:n))
allocate(lumtot_clustercore(1:n))
allocate(lumtot_err(1:n))
allocate(lumobs_clustercore(1:n))
allocate(mass_clustercore(1:n))
allocate(mass_err(1:n))
allocate(mdyn_clustercore(1:n))
allocate(sigma_clustercore(1:n))
allocate(radius_clustercore(1:n))
allocate(angradius_clustercore(1:n))
allocate(lumdistance_clustercore(1:n))
allocate(members_clustercore(1:n))
allocate(sdss_gal_b(1:n))
allocate(sdss_gal_l(1:n))


OPEN(50,file='catalogues/cluster_list_all_SDSS.txt')
DO i=1,n
READ(50,*) iii,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
lumtot_clustercore(i),lumtot_err(i),lumobs_clustercore(i),&
mass_clustercore(i),mass_err(i),mdyn_clustercore(i),&
sigma_clustercore(i),radius_clustercore(i),angradius_clustercore(i),lumdistance_clustercore(i),&
members_clustercore(i)
! WRITE(*,*) ra_clustercore(i),dec_clustercore(i)
END DO
 CLOSE(50)
 

! get length of file
OPEN(50,file='catalogues/cluster_list_all_2MRS.txt')
io_err=0
n2=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n2=n2+1
END DO
 CLOSE(50)
n2=n2-1

WRITE(*,*) n2,'galaxies will be used'


! 
allocate(ra_clustercore2(1:n2))
allocate(dec_clustercore2(1:n2))
allocate(z_clustercore2(1:n2))
allocate(lumtot_clustercore2(1:n2))
allocate(lumtot_err2(1:n2))
allocate(lumobs_clustercore2(1:n2))
allocate(mass_clustercore2(1:n2))
allocate(mass_err2(1:n2))
allocate(mdyn_clustercore2(1:n2))
allocate(sigma_clustercore2(1:n2))
allocate(radius_clustercore2(1:n2))
allocate(angradius_clustercore2(1:n2))
allocate(lumdistance_clustercore2(1:n2))
allocate(members_clustercore2(1:n2))
allocate(twomass_gal_b(1:n))
allocate(twomass_gal_l(1:n))


OPEN(50,file='catalogues/cluster_list_all_2MRS.txt')
DO i=1,n2
READ(50,*) iii,ra_clustercore2(i),dec_clustercore2(i),z_clustercore2(i),&
lumtot_clustercore2(i),lumtot_err2(i),lumobs_clustercore2(i),&
mass_clustercore2(i),mass_err2(i),mdyn_clustercore2(i),&
sigma_clustercore2(i),radius_clustercore2(i),angradius_clustercore2(i),lumdistance_clustercore2(i),&
members_clustercore2(i)
END DO
 CLOSE(50)
 
 
 
 
 
 
 
 
OPEN(50,file='catalogues/cluster_fp_distances_SDSS_m.txt')
io_err=0
n_fp=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_fp=n_fp+1
END DO
 CLOSE(50)
n_fp=n_fp-1


WRITE(*,*) n_fp,'galaxies will be used'


allocate(fp_ra_clustercore(1:n_fp))
allocate(fp_dec_clustercore(1:n_fp))
allocate(fp_z_clustercore(1:n_fp))
allocate(fp_err_z_m(1:n_fp))
allocate(cluster_fp_adist_m(1:n_fp))
allocate(err_fpd_m(1:n_fp))
allocate(cluster_fp_cdist_m(1:n_fp))
allocate(cluster_fp_ldist_m(1:n_fp))
allocate(cluster_z_adist_m(1:n_fp))
allocate(err_red_m(1:n_fp))
allocate(cluster_z_cdist_m(1:n_fp))
allocate(cluster_z_ldist_m(1:n_fp))
allocate(n_fp_in_cluster(1:n_fp))
allocate(fp_members_clustercore(1:n_fp))
allocate(fp_id(1:n_fp))
allocate(fp_gal_b(1:n))
allocate(fp_gal_l(1:n))


OPEN(50,file='catalogues/cluster_fp_distances_SDSS_m.txt')
DO i=1,n_fp

READ(50,*) fp_id(i),fp_ra_clustercore(i),fp_dec_clustercore(i),fp_z_clustercore(i),fp_err_z_m(i),&
 cluster_fp_adist_m(i),err_fpd_m(i),&
 cluster_fp_cdist_m(i),dummy_var,&
 cluster_fp_ldist_m(i),dummy_var,&
 cluster_z_adist_m(i),err_red_m(i),&
 cluster_z_cdist_m(i),dummy_var,&
 cluster_z_ldist_m(i),dummy_var,&
 n_fp_in_cluster(i),fp_members_clustercore(i)
!  WRITE(*,*) i
err_fpd_m(i)=err_fpd_m(i)/cluster_fp_adist_m(i)
err_red_m(i)=err_red_m(i)/cluster_z_adist_m(i)
END DO
 CLOSE(50)

 
 WRITE(*,*) 'all data read in'
 
 
 
 
 
 
  v_mws=-220.D0
  l_mws=270.D0
  b_mws=0.D0
  
  z_mws=v_mws/light
gal_b_mws=b_mws*PI/180.D0
gal_l_mws=l_mws*PI/180.D0
z_mws_x=z_mws*COS(gal_b_mws)*COS(gal_l_mws)
z_mws_y=z_mws*COS(gal_b_mws)*SIN(gal_l_mws)
z_mws_z=z_mws*SIN(gal_b_mws)


 
 
!prepare for calculation of correction for CMB motion
! values for the motion of the sun relativ to the CMB
v_cmb=-369.D0 ! ANTI-correction
l_cmb=263.99D0
b_cmb=48.26D0

! calculate sun's motion in z space relativ to the CMB
z_cmb=v_cmb/light
gal_b=b_cmb*PI/180.D0
gal_l=l_cmb*PI/180.D0
z_cmb_x=z_cmb*COS(gal_b)*COS(gal_l)
z_cmb_y=z_cmb*COS(gal_b)*SIN(gal_l)
z_cmb_z=z_cmb*SIN(gal_b)


xhelp1=62.6*PI/180.D0
xhelp2=282.25*PI/180.D0
xhelp3=33.0*PI/180.D0

DO i=1,n
dec_help=dec_clustercore(i)*PI/180.D0
ra_help=ra_clustercore(i)*PI/180.D0
! WRITE(*,*) dec_help,ra_help

dummy_var=SIN(dec_help)*COS(xhelp1)-COS(dec_help)*SIN(ra_help-xhelp2)*SIN(xhelp1)
sdss_gal_b(i)=ASIN(dummy_var)

dummy_var=COS(dec_help)*COS(ra_help-xhelp2)/COS(sdss_gal_b(i))
sdss_gal_l(i)=(ACOS(dummy_var)+xhelp3)

sdss_gal_b(i)=180.D0/PI*sdss_gal_b(i)
sdss_gal_l(i)=180.D0/PI*sdss_gal_l(i)
! WRITE(*,*) sdss_gal_b(i),sdss_gal_l(i)
END DO


DO i=1,n2
dec_help=dec_clustercore2(i)*PI/180.D0
ra_help=ra_clustercore2(i)*PI/180.D0

dummy_var=SIN(dec_help)*COS(xhelp1)-COS(dec_help)*SIN(ra_help-xhelp2)*SIN(xhelp1)
twomass_gal_b(i)=ASIN(dummy_var)

dummy_var=COS(dec_help)*COS(ra_help-xhelp2)/COS(twomass_gal_b(i))
twomass_gal_l(i)=(ACOS(dummy_var)+xhelp3)

twomass_gal_b(i)=180.D0/PI*twomass_gal_b(i)
twomass_gal_l(i)=180.D0/PI*twomass_gal_l(i)
END DO


DO i=1,n_fp
dec_help=fp_dec_clustercore(i)*PI/180.D0
ra_help=fp_ra_clustercore(i)*PI/180.D0

dummy_var=SIN(dec_help)*COS(xhelp1)-COS(dec_help)*SIN(ra_help-xhelp2)*SIN(xhelp1)
fp_gal_b(i)=ASIN(dummy_var)

dummy_var=COS(dec_help)*COS(ra_help-xhelp2)/COS(fp_gal_b(i))
fp_gal_l(i)=(ACOS(dummy_var)+xhelp3)

fp_gal_b(i)=180.D0/PI*fp_gal_b(i)
fp_gal_l(i)=180.D0/PI*fp_gal_l(i)
END DO



DO i=1,n
z_help=z_clustercore(i)

gal_b=sdss_gal_b(i)*PI/180.D0
gal_l=sdss_gal_l(i)*PI/180.D0
z_gal_x=z_help*COS(gal_b)*COS(gal_l)
z_gal_y=z_help*COS(gal_b)*SIN(gal_l)
z_gal_z=z_help*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
z_clustercore(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


DO i=1,n2
z_help=z_clustercore2(i)

gal_b=twomass_gal_b(i)*PI/180.D0
gal_l=twomass_gal_l(i)*PI/180.D0
z_gal_x=z_help*COS(gal_b)*COS(gal_l)
z_gal_y=z_help*COS(gal_b)*SIN(gal_l)
z_gal_z=z_help*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
z_clustercore2(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


 DO i=1,n_fp
z_help=fp_z_clustercore(i)

gal_b=fp_gal_b(i)*PI/180.D0
gal_l=fp_gal_l(i)*PI/180.D0
z_gal_x=z_help*COS(gal_b)*COS(gal_l)
z_gal_y=z_help*COS(gal_b)*SIN(gal_l)
z_gal_z=z_help*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
fp_z_clustercore(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
    END DO
 
 
 
 
DO i=1,n
z_help=z_clustercore(i)

gal_b=sdss_gal_b(i)*PI/180.D0
gal_l=sdss_gal_l(i)*PI/180.D0
z_gal_x=z_help*COS(gal_b)*COS(gal_l)
z_gal_y=z_help*COS(gal_b)*SIN(gal_l)
z_gal_z=z_help*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_mws_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_mws_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_mws_z))-1.D0
z_clustercore(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


DO i=1,n2
z_help=z_clustercore2(i)

gal_b=twomass_gal_b(i)*PI/180.D0
gal_l=twomass_gal_l(i)*PI/180.D0
z_gal_x=z_help*COS(gal_b)*COS(gal_l)
z_gal_y=z_help*COS(gal_b)*SIN(gal_l)
z_gal_z=z_help*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_mws_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_mws_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_mws_z))-1.D0
z_clustercore2(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


 DO i=1,n_fp
z_help=fp_z_clustercore(i)

gal_b=fp_gal_b(i)*PI/180.D0
gal_l=fp_gal_l(i)*PI/180.D0
z_gal_x=z_help*COS(gal_b)*COS(gal_l)
z_gal_y=z_help*COS(gal_b)*SIN(gal_l)
z_gal_z=z_help*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_mws_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_mws_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_mws_z))-1.D0
fp_z_clustercore(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
    END DO
 
 
 
 
 DO i=1,n
!luminosity distance in Mpc!
IF (z_clustercore(i)>0.D0) THEN
divisor=(SQRT(1.D0+2.D0*q0*z_clustercore(i))+1.D0+q0*z_clustercore(i))
lumdistance_clustercore(i)=light/H0*z_clustercore(i)*(1.D0+((z_clustercore(i)*(1.D0-q0))/divisor))
ELSE
lumdistance_clustercore(i)=0.D0
END IF
END DO

!  
  DO i=1,n2
!luminosity distance in Mpc!
IF (z_clustercore2(i)>0.D0) THEN
divisor=(SQRT(1.D0+2.D0*q0*z_clustercore2(i))+1.D0+q0*z_clustercore2(i))
lumdistance_clustercore2(i)=light/H0*z_clustercore2(i)*(1.D0+((z_clustercore2(i)*(1.D0-q0))/divisor))
ELSE
lumdistance_clustercore2(i)=0.D0
END IF
END DO



    DO i=1,n_fp
   IF (fp_z_clustercore(i)>0.D0) THEN 
 divisor=(SQRT(1.D0+2.D0*q0*fp_z_clustercore(i))+1.D0+q0*fp_z_clustercore(i))
 cluster_z_ldist_m(i)=light/H0*fp_z_clustercore(i)*(1.D0+((fp_z_clustercore(i)*(1.D0-q0))/divisor))
 cluster_z_cdist_m(i)=cluster_z_ldist_m(i)*((1.D0+fp_z_clustercore(i))**(-1.D0))
 cluster_z_adist_m(i)=cluster_z_ldist_m(i)*((1.D0+fp_z_clustercore(i))**(-2.D0)) 
 ELSE
 cluster_z_ldist_m(i)=0.D0
 cluster_z_cdist_m(i)=0.D0
 cluster_z_adist_m(i)=0.D0
 END IF
 END DO 
 



OPEN(50,file='ts_sor/catalogues/cluster_list_all_SDSS.txt')
DO i=1,n
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
lumtot_clustercore(i),lumtot_err(i),lumobs_clustercore(i),&
mass_clustercore(i),mass_err(i),mdyn_clustercore(i),&
sigma_clustercore(i),radius_clustercore(i),angradius_clustercore(i),lumdistance_clustercore(i),&
members_clustercore(i)
END DO
 CLOSE(50)
 
 
 
 


OPEN(50,file='ts_sor/catalogues/cluster_list_all_2MRS.txt')
DO i=1,n2
WRITE(50,*) i,ra_clustercore2(i),dec_clustercore2(i),z_clustercore2(i),&
lumtot_clustercore2(i),lumtot_err2(i),lumobs_clustercore2(i),&
mass_clustercore2(i),mass_err2(i),mdyn_clustercore2(i),&
sigma_clustercore2(i),radius_clustercore2(i),angradius_clustercore2(i),lumdistance_clustercore2(i),&
members_clustercore2(i)
END DO
 CLOSE(50)
 
 
 
 
 
 OPEN(50,file='ts_sor/catalogues/cluster_fp_distances_SDSS_m.txt')
DO i=1,n_fp
WRITE(50,*) fp_id(i),fp_ra_clustercore(i),fp_dec_clustercore(i),fp_z_clustercore(i),fp_err_z_m(i),&
 cluster_fp_adist_m(i),(err_fpd_m(i)*cluster_fp_adist_m(i)),&
 cluster_fp_cdist_m(i),(err_fpd_m(i)*cluster_fp_cdist_m(i)),&
 cluster_fp_ldist_m(i),(err_fpd_m(i)*cluster_fp_ldist_m(i)),&
 cluster_z_adist_m(i),(err_red_m(i)*cluster_z_adist_m(i)),&
 cluster_z_cdist_m(i),(err_red_m(i)*cluster_z_cdist_m(i)),&
 cluster_z_ldist_m(i),(err_red_m(i)*cluster_z_ldist_m(i)),&
n_fp_in_cluster(i),fp_members_clustercore(i)
END DO
 CLOSE(50)

 
 
 
 

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


 
 
 
 
 
