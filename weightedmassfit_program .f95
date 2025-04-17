PROGRAM recalibration
! declaration of variables
IMPLICIT NONE
! basic constants
integer, parameter :: dim_s=7
integer, parameter :: dim_ml=11
integer, parameter :: dim_mh=13

integer :: io_err,n_2mrs_s,n_2mrs_m,n_sdss_s,n_sdss_m,i,ii,iii,hc1,hc2,uuu,n_sdss_tot,n_2mrs_tot
integer :: n_2mrs_mh,n_sdss_mh,n_2mrs_ml,n_sdss_ml

double precision, allocatable :: mass_2mrs_s(:),lum_2mrs_s(:),lum2_2mrs_s(:),dist_2mrs_s(:)
double precision, allocatable :: mass_sdss_s(:),lum_sdss_s(:),lum2_sdss_s(:),dist_sdss_s(:)
double precision, allocatable :: mass_2mrs_m(:),lum_2mrs_m(:),lum2_2mrs_m(:),dist_2mrs_m(:)
double precision, allocatable :: sigma_2mrs_m(:),rad_2mrs_m(:),nvis_2mrs_m(:)
double precision, allocatable :: mass_sdss_m(:),lum_sdss_m(:),lum2_sdss_m(:),dist_sdss_m(:)
double precision, allocatable :: sigma_sdss_m(:),rad_sdss_m(:),nvis_sdss_m(:)
double precision, allocatable :: colour_sdss_m(:),colour_sdss_s(:),colour_2mrs_m(:),colour_2mrs_s(:)
double precision, allocatable :: lumcor_2mrs_s(:),lumcor_sdss_s(:),lumcor_2mrs_m(:),lumcor_sdss_m(:)
double precision, allocatable :: mdyn_2mrs(:),mdyn_sdss(:)
double precision, dimension(1:dim_s) :: solution_sdss_s,solution_2mrs_s,errorbar_sdss_s,errorbar_2mrs_s
double precision, dimension(1:dim_mh) :: solution_sdss_mh,solution_2mrs_mh,errorbar_sdss_mh,errorbar_2mrs_mh
double precision, dimension(1:dim_ml) :: solution_sdss_ml,solution_2mrs_ml,errorbar_sdss_ml,errorbar_2mrs_ml
double precision, allocatable :: datalist_sdss_s(:,:),datalist_2mrs_s(:,:)
double precision, allocatable :: datalist_sdss_mh(:,:),datalist_2mrs_mh(:,:)
double precision, allocatable :: datalist_sdss_ml(:,:),datalist_2mrs_ml(:,:)
double precision :: rms_sdss_s,rms_2mrs_s,rms_sdss_mh,rms_2mrs_mh,rms_sdss_ml,rms_2mrs_ml

double precision, allocatable :: deviation_sdss_s(:),deviation_2mrs_s(:)
double precision, allocatable :: otherside_sdss_s(:),otherside_2mrs_s(:)
double precision, allocatable :: deviation_sdss_mh(:),deviation_2mrs_mh(:)
double precision, allocatable :: otherside_sdss_mh(:),otherside_2mrs_mh(:)
double precision, allocatable :: deviation_sdss_ml(:),deviation_2mrs_ml(:)
double precision, allocatable :: otherside_sdss_ml(:),otherside_2mrs_ml(:)

double precision :: mass_fit_sum_sdss_s,mass_fit_sum_2mrs_s,mass_orig_sum_sdss_s,mass_orig_sum_2mrs_s
double precision :: mass_fit_sum_sdss_ml,mass_fit_sum_2mrs_ml,mass_orig_sum_sdss_ml,mass_orig_sum_2mrs_ml
double precision :: mass_fit_sum_sdss_mh,mass_fit_sum_2mrs_mh,mass_orig_sum_sdss_mh,mass_orig_sum_2mrs_mh

double precision :: mass_fit_sum_sdss_tot,mass_fit_sum_2mrs_tot
double precision :: mass_orig_sum_sdss_tot,mass_orig_sum_2mrs_tot
double precision :: mag_sol_r,mag_sol_g,mag_sol_Ks,mag_sol_J,area_cover

double precision, allocatable :: ra_sdss(:),dec_sdss(:),z_sdss(:),mag_g_sdss(:),mag_r_sdss(:)
double precision, allocatable :: ra_2mrs(:),dec_2mrs(:),z_2mrs(:),mag_g_2mrs(:),mag_r_2mrs(:)
logical, allocatable :: active_sdss(:),active_2mrs(:)
double precision :: V_sum_sdss,V_sum_2mrs,D_limit,magmax_sdss,magmin_sdss,magmax_2mrs,magmin_2mrs
double precision :: mag_limit_sdss,mag_limit_2mrs,z_limit,V_limit,mabs_vollim,lum_sum_sdss,lum_sum_2mrs
double precision, allocatable :: volume_weight_sdss(:),volume_weight_2mrs(:),lum_dist_sdss(:),lum_dist_2mrs(:)
double precision, allocatable :: vol_int_sdss(:),vol_int_2mrs(:),lum_gal_sdss(:),lum_gal_2mrs(:)
double precision, allocatable :: lum_int_sdss(:),lum_int_2mrs(:),lum_true_sdss(:),lum_true_2mrs(:)
double precision ::t1,t2,t3

double precision, dimension(0:500,0:500) :: mapmap
real :: help_x2,help_y2
integer :: help_x,help_y,calculate_lumweights
character(20) :: appendix,mockstr
character(200) :: filename_sdss,filename_2mrs
integer, dimension(1:8) :: n_mock_sdss,n_mock_2mrs
integer :: n_mock_sdss_all,n_mock_2mrs_all


double precision :: Omega_m,Omega_l,H0,q0,light,h,G,rho_crit,vol,m_expected,z_max,d_com,divisor,PI


! ! define constants
PI=ACOS(-1.D0)

Omega_m=0.25D0
Omega_l=0.75D0
H0=73.D0


 calculate_lumweights=0

! sdss_cover=9274.D0/41253.D0
! twomrs_cover=0.91D0

q0=Omega_m/2.D0-Omega_l
light=3.D5
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
 
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)

z_max=0.11D0


divisor=(SQRT(1.D0+2.D0*q0*z_max)+1.D0+q0*z_max)
d_com=light/H0*z_max*(1.D0+((z_max*(1.D0-q0))/divisor))
d_com=d_com/(1.D0+z_max)*1.D6

mag_sol_r=4.71D0
mag_sol_g=5.31D0     

mag_sol_Ks=3.28D0 !Ks
mag_sol_J=3.64D0 !J

vol=4.D0*PI/3.D0*(d_com**3)

m_expected=vol*rho_crit*Omega_m

area_cover=1.D0


!  cV_mill=(380.D0/h)**3
! m_expected=cV_mill*rho_crit*Omega_m

! WRITE(*,*) m_expected
! 
! mock_cover=1.D0/8.D0
! 

! get length of file
OPEN(50,file='mass/log_mass_dependences_multi_2MRS_combi.txt')
io_err=0
n_2mrs_m=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mrs_m=n_2mrs_m+1
END DO
 CLOSE(50)
n_2mrs_m=n_2mrs_m-1



! get length of file
OPEN(50,file='mass/log_mass_dependences_single_2MRS_combi.txt')
io_err=0
n_2mrs_s=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mrs_s=n_2mrs_s+1
END DO
 CLOSE(50)
n_2mrs_s=n_2mrs_s-1



! get length of file
OPEN(50,file='mass/log_mass_dependences_multi_SDSS_combi.txt')
io_err=0
n_sdss_m=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss_m=n_sdss_m+1
END DO
 CLOSE(50)
n_sdss_m=n_sdss_m-1



! get length of file
OPEN(50,file='mass/log_mass_dependences_single_SDSS_combi.txt')
io_err=0
n_sdss_s=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss_s=n_sdss_s+1
END DO
 CLOSE(50)
n_sdss_s=n_sdss_s-1


! allocate arrays
n_sdss_tot=n_sdss_s+n_sdss_m
n_2mrs_tot=n_2mrs_s+n_2mrs_m

allocate(mass_2mrs_s(1:n_2mrs_s))
allocate(lum_2mrs_s(1:n_2mrs_s))
allocate(dist_2mrs_s(1:n_2mrs_s))
allocate(colour_2mrs_s(1:n_2mrs_s))


allocate(mass_sdss_s(1:n_sdss_s))
allocate(lum_sdss_s(1:n_sdss_s))
allocate(dist_sdss_s(1:n_sdss_s))
allocate(colour_sdss_s(1:n_sdss_s))


allocate(mass_2mrs_m(1:n_2mrs_m))
allocate(lum_2mrs_m(1:n_2mrs_m))
allocate(dist_2mrs_m(1:n_2mrs_m))
allocate(sigma_2mrs_m(1:n_2mrs_m))
allocate(rad_2mrs_m(1:n_2mrs_m))
allocate(nvis_2mrs_m(1:n_2mrs_m))
allocate(colour_2mrs_m(1:n_2mrs_m))


allocate(mass_sdss_m(1:n_sdss_m))
allocate(lum_sdss_m(1:n_sdss_m))
allocate(dist_sdss_m(1:n_sdss_m))
allocate(sigma_sdss_m(1:n_sdss_m))
allocate(rad_sdss_m(1:n_sdss_m))
allocate(nvis_sdss_m(1:n_sdss_m))
allocate(colour_sdss_m(1:n_sdss_m))


allocate(mdyn_2mrs(1:n_2mrs_m))
allocate(mdyn_sdss(1:n_sdss_m))

allocate(lumcor_2mrs_s(1:n_2mrs_s))
allocate(lumcor_sdss_s(1:n_sdss_s))
allocate(lumcor_2mrs_m(1:n_2mrs_m))
allocate(lumcor_sdss_m(1:n_sdss_m))

!  dim_s=6
!  dim_m=9
!  

! read file
OPEN(50,file='mass/log_mass_dependences_single_2MRS_combi.txt')
DO i=1,n_2mrs_s
READ(50,*) mass_2mrs_s(i),lum_2mrs_s(i),dist_2mrs_s(i),colour_2mrs_s(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_2mrs_s,' galaxies'


! read file
OPEN(50,file='mass/log_mass_dependences_single_SDSS_combi.txt')
DO i=1,n_sdss_s
READ(50,*) mass_sdss_s(i),lum_sdss_s(i),dist_sdss_s(i),colour_sdss_s(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_sdss_s,' galaxies'


! read file
OPEN(50,file='mass/log_mass_dependences_multi_2MRS_combi.txt')
DO i=1,n_2mrs_m
READ(50,*) mass_2mrs_m(i),lum_2mrs_m(i),dist_2mrs_m(i),colour_2mrs_m(i),&
sigma_2mrs_m(i),rad_2mrs_m(i),nvis_2mrs_m(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_2mrs_m,' galaxies'


! read file
OPEN(50,file='mass/log_mass_dependences_multi_SDSS_combi.txt')
DO i=1,n_sdss_m
READ(50,*) mass_sdss_m(i),lum_sdss_m(i),dist_sdss_m(i),colour_sdss_m(i),&
sigma_sdss_m(i),rad_sdss_m(i),nvis_sdss_m(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_sdss_m,' galaxies'


 



DO uuu=1,8 

IF (uuu==1) THEN
mockstr='mock1'
END IF
IF (uuu==2) THEN
mockstr='mock2'
END IF
IF (uuu==3) THEN
mockstr='mock3'
END IF
IF (uuu==4) THEN
mockstr='mock4'
END IF
IF (uuu==5) THEN
mockstr='mock5'
END IF
IF (uuu==6) THEN
mockstr='mock6'
END IF
IF (uuu==7) THEN
mockstr='mock7'
END IF
IF (uuu==8) THEN
mockstr='mock8'
END IF
 
  mockstr=TRIM(adjustl(mockstr))

appendix='SDSS'
appendix=TRIM(adjustl(appendix))
filename_sdss=TRIM(appendix)//'_'//TRIM(mockstr)//'/'//TRIM(appendix)//'_galaxies_'//TRIM(mockstr)//'_final.txt'

appendix='2MRS'
appendix=TRIM(adjustl(appendix))
filename_2mrs=TRIM(appendix)//'_'//TRIM(mockstr)//'/'//TRIM(appendix)//'_galaxies_'//TRIM(mockstr)//'_final.txt'




! get length of file
OPEN(50,file=filename_sdss)
io_err=0
n_mock_sdss(uuu)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock_sdss(uuu)=n_mock_sdss(uuu)+1
END DO
 CLOSE(50)
n_mock_sdss(uuu)=n_mock_sdss(uuu)-1

! get length of file
OPEN(50,file=filename_2mrs)
io_err=0
n_mock_2mrs(uuu)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock_2mrs(uuu)=n_mock_2mrs(uuu)+1
END DO
 CLOSE(50)
n_mock_2mrs(uuu)=n_mock_2mrs(uuu)-1

END DO

n_mock_sdss_all=0
n_mock_2mrs_all=0

DO uuu=1,8
n_mock_sdss_all=n_mock_sdss_all+n_mock_sdss(uuu)
n_mock_2mrs_all=n_mock_2mrs_all+n_mock_2mrs(uuu)
END DO



allocate(ra_sdss(1:n_mock_sdss_all))
allocate(dec_sdss(1:n_mock_sdss_all))
allocate(z_sdss(1:n_mock_sdss_all))
allocate(mag_g_sdss(1:n_mock_sdss_all))
allocate(mag_r_sdss(1:n_mock_sdss_all))
allocate(ra_2mrs(1:n_mock_2mrs_all))
allocate(dec_2mrs(1:n_mock_2mrs_all))
allocate(z_2mrs(1:n_mock_2mrs_all))
allocate(mag_g_2mrs(1:n_mock_2mrs_all))
allocate(mag_r_2mrs(1:n_mock_2mrs_all))


hc1=0
hc2=0

DO uuu=1,8 

IF (uuu==1) THEN
mockstr='mock1'
END IF
IF (uuu==2) THEN
mockstr='mock2'
END IF
IF (uuu==3) THEN
mockstr='mock3'
END IF
IF (uuu==4) THEN
mockstr='mock4'
END IF
IF (uuu==5) THEN
mockstr='mock5'
END IF
IF (uuu==6) THEN
mockstr='mock6'
END IF
IF (uuu==7) THEN
mockstr='mock7'
END IF
IF (uuu==8) THEN
mockstr='mock8'
END IF
 
  mockstr=TRIM(adjustl(mockstr))

appendix='SDSS'
appendix=TRIM(adjustl(appendix))
filename_sdss=TRIM(appendix)//'_'//TRIM(mockstr)//'/'//TRIM(appendix)//'_galaxies_'//TRIM(mockstr)//'_final.txt'

appendix='2MRS'
appendix=TRIM(adjustl(appendix))
filename_2mrs=TRIM(appendix)//'_'//TRIM(mockstr)//'/'//TRIM(appendix)//'_galaxies_'//TRIM(mockstr)//'_final.txt'



OPEN(50,file=filename_sdss)
DO i=1,n_mock_sdss(uuu)
READ(50,*) ra_sdss(i+hc1),dec_sdss(i+hc1),z_sdss(i+hc1),mag_g_sdss(i+hc1),mag_r_sdss(i+hc1)
END DO
 CLOSE(50)



OPEN(50,file=filename_2mrs)
DO i=1,n_mock_2mrs(uuu)
READ(50,*) ra_2mrs(i+hc2),dec_2mrs(i+hc2),z_2mrs(i+hc2),mag_g_2mrs(i+hc2),mag_r_2mrs(i+hc2)
END DO
 CLOSE(50)



hc1=hc1+n_mock_sdss(uuu)
hc2=hc2+n_mock_2mrs(uuu)




END DO

WRITE(*,*) 'all data loaded'
allocate(active_sdss(1:n_mock_sdss_all))
allocate(active_2mrs(1:n_mock_2mrs_all))
allocate(volume_weight_sdss(1:n_mock_sdss_all))
allocate(volume_weight_2mrs(1:n_mock_2mrs_all))

allocate(vol_int_sdss(1:n_mock_sdss_all))
allocate(vol_int_2mrs(1:n_mock_2mrs_all))
allocate(lum_gal_sdss(1:n_mock_sdss_all))
allocate(lum_gal_2mrs(1:n_mock_2mrs_all))

allocate(lum_dist_sdss(1:n_sdss_tot))
allocate(lum_dist_2mrs(1:n_2mrs_tot))

allocate(lum_int_sdss(1:n_sdss_tot))
allocate(lum_int_2mrs(1:n_2mrs_tot))


! 
magmax_sdss=-30.D0
magmin_sdss=-15.D0
magmax_2mrs=-30.D0
magmin_2mrs=-18.D0
mag_limit_sdss=17.77D0
mag_limit_2mrs=11.75D0



! calculate volume weights to corret for Malmquist bias
V_sum_sdss=0.D0
lum_sum_sdss=0.D0
DO i=1,n_mock_sdss_all

lum_gal_sdss(i)=10.D0**(-0.4D0*(mag_r_sdss(i)-mag_sol_r)) 

active_sdss(i)=.TRUE.
volume_weight_sdss(i)=0.D0

IF (mag_r_sdss(i)<magmax_sdss) THEN
active_sdss(i)=.FALSE.
END IF

IF (mag_r_sdss(i)>magmin_sdss) THEN
active_sdss(i)=.FALSE.
END IF

IF (mag_g_sdss(i)<(magmax_sdss-1.5D0)) THEN
active_sdss(i)=.FALSE.
END IF

IF (mag_g_sdss(i)>(magmin_sdss+1.5D0)) THEN
active_sdss(i)=.FALSE.
END IF


IF ((mag_g_sdss(i)-mag_r_sdss(i))>2.D0) THEN
active_sdss(i)=.FALSE.
END IF

IF ((mag_g_sdss(i)-mag_r_sdss(i))<-1.D0) THEN
active_sdss(i)=.FALSE.
END IF


IF (active_sdss(i)) THEN
D_limit=-0.2D0*mag_r_sdss(i)+((mag_limit_sdss+5.D0)/5.D0)
D_limit=10.D0**D_limit
D_limit=D_limit/1.D6



t1=(light**2)*q0-(light**2)
t2=(light**4)*(q0**2)-2.D0*(light**4)*q0+(light**4)+2.D0*(light**3)*D_limit*H0*(q0**2)
t2=t2-4.D0*(light**3)*D_limit*H0*q0+2.D0*(light**3)*D_limit*H0
t2=SQRT(t2)
t3=light*D_limit*H0*q0

z_limit=1/(light**2)*(t1+t2+t3)
 
D_limit=D_limit*((1.D0+z_limit)**(-1.D0))!comoving diameter distance


IF (z_limit>0.11D0) THEN
z_limit=0.11D0
divisor=(SQRT(1.D0+2.D0*q0*z_limit)+1.D0+q0*z_limit)
D_limit=light/H0*z_limit*(1.D0+((z_limit*(1.D0-q0))/divisor))
D_limit=D_limit*((1.D0+z_limit)**(-1.D0))

END IF


V_limit=4.D0/3.D0*PI*(D_limit**3)

V_limit=V_limit*area_cover


volume_weight_sdss(i)=1.D0/V_limit

V_sum_sdss=V_sum_sdss+volume_weight_sdss(i)

lum_gal_sdss(i)=lum_gal_sdss(i)*volume_weight_sdss(i)
lum_sum_sdss=lum_sum_sdss+lum_gal_sdss(i)

END IF

END DO 





! calculate volume weights to corret for Malmquist bias
V_sum_2mrs=0.D0
lum_sum_2mrs=0.D0
DO i=1,n_mock_2mrs_all

lum_gal_2mrs(i)=10.D0**(-0.4D0*(mag_r_2mrs(i)-mag_sol_Ks)) 
active_2mrs(i)=.TRUE.
volume_weight_2mrs(i)=0.D0

IF (mag_r_2mrs(i)<magmax_2mrs) THEN
active_2mrs(i)=.FALSE.
END IF

IF (mag_r_2mrs(i)>magmin_2mrs) THEN
active_2mrs(i)=.FALSE.
END IF

IF (mag_g_2mrs(i)<(magmax_2mrs-1.5D0)) THEN
active_2mrs(i)=.FALSE.
END IF

IF (mag_g_2mrs(i)>(magmin_2mrs+1.5D0)) THEN
active_2mrs(i)=.FALSE.
END IF


IF ((mag_g_2mrs(i)-mag_r_2mrs(i))>2.D0) THEN
active_2mrs(i)=.FALSE.
END IF

IF ((mag_g_2mrs(i)-mag_r_2mrs(i))<-1.D0) THEN
active_2mrs(i)=.FALSE.
END IF


IF (active_2mrs(i)) THEN
D_limit=-0.2D0*mag_r_2mrs(i)+((mag_limit_2mrs+5.D0)/5.D0)
D_limit=10.D0**D_limit
D_limit=D_limit/1.D6



t1=(light**2)*q0-(light**2)
t2=(light**4)*(q0**2)-2.D0*(light**4)*q0+(light**4)+2.D0*(light**3)*D_limit*H0*(q0**2)
t2=t2-4.D0*(light**3)*D_limit*H0*q0+2.D0*(light**3)*D_limit*H0
t2=SQRT(t2)
t3=light*D_limit*H0*q0

z_limit=1/(light**2)*(t1+t2+t3)
 
D_limit=D_limit*((1.D0+z_limit)**(-1.D0))!comoving diameter distance


IF (z_limit>0.11D0) THEN
z_limit=0.11D0
divisor=(SQRT(1.D0+2.D0*q0*z_limit)+1.D0+q0*z_limit)
D_limit=light/H0*z_limit*(1.D0+((z_limit*(1.D0-q0))/divisor))
D_limit=D_limit*((1.D0+z_limit)**(-1.D0))

END IF


V_limit=4.D0/3.D0*PI*(D_limit**3)

V_limit=V_limit*area_cover


volume_weight_2mrs(i)=1.D0/V_limit

V_sum_2mrs=V_sum_2mrs+volume_weight_2mrs(i)

lum_gal_2mrs(i)=lum_gal_2mrs(i)*volume_weight_2mrs(i)
lum_sum_2mrs=lum_sum_2mrs+lum_gal_2mrs(i)

END IF

END DO 



DO i=1,n_mock_sdss_all
IF (active_sdss(i)) THEN
lum_gal_sdss(i)=lum_gal_sdss(i)/lum_sum_sdss
END IF
END DO

DO i=1,n_mock_2mrs_all
IF (active_2mrs(i)) THEN
lum_gal_2mrs(i)=lum_gal_2mrs(i)/lum_sum_2mrs
END IF
END DO

WRITE(*,*) 'weights calculated'
 
 
IF (calculate_lumweights==1) THEN
 
DO i=1,n_sdss_s
lum_dist_sdss(i)=dist_sdss_s(i)+6.D0
END DO

DO i=1,n_sdss_m
lum_dist_sdss(i+n_sdss_s)=dist_sdss_m(i)+6.D0
END DO


DO i=1,n_2mrs_s
lum_dist_2mrs(i)=dist_2mrs_s(i)+6.D0
END DO

DO i=1,n_2mrs_m
lum_dist_2mrs(i+n_2mrs_s)=dist_2mrs_m(i)+6.D0
END DO


DO i=1,n_sdss_tot
WRITE(*,*) i,'SDSS'
lum_int_sdss(i)=0.D0
mabs_vollim=mag_limit_sdss+5.D0-5.D0*lum_dist_sdss(i)
DO ii=1,n_mock_sdss_all
IF (active_sdss(ii)) THEN
IF (mag_r_sdss(ii)<mabs_vollim) THEN
lum_int_sdss(i)=lum_int_sdss(i)+lum_gal_sdss(ii)
END IF
END IF
END DO


END DO


DO i=1,n_2mrs_tot
WRITE(*,*) i,'2MRS'
lum_int_2mrs(i)=0.D0 
mabs_vollim=mag_limit_2mrs+5.D0-5.D0*lum_dist_2mrs(i)
DO ii=1,n_mock_2mrs_all
IF (active_2mrs(ii)) THEN
IF (mag_r_2mrs(ii)<mabs_vollim) THEN
lum_int_2mrs(i)=lum_int_2mrs(i)+lum_gal_2mrs(ii)
END IF
END IF
END DO


END DO


  OPEN(61,file='mass/lum_int_2mrs.txt')
DO i=1,n_2mrs_tot
WRITE(61,*) lum_int_2mrs(i),lum_dist_2mrs(i)
END DO
 CLOSE(61)
  
  OPEN(61,file='mass/lum_int_sdss.txt')
DO i=1,n_sdss_tot
WRITE(61,*) lum_int_sdss(i),lum_dist_sdss(i)
END DO
 CLOSE(61)

 
 ELSE
 
  OPEN(61,file='mass/lum_int_2mrs.txt')
DO i=1,n_2mrs_tot
READ(61,*) lum_int_2mrs(i),lum_dist_2mrs(i)
END DO
 CLOSE(61)
  
  OPEN(61,file='mass/lum_int_sdss.txt')
DO i=1,n_sdss_tot
READ(61,*) lum_int_sdss(i),lum_dist_sdss(i)
END DO
 CLOSE(61)
 
END IF





DO i=1,n_sdss_s
lumcor_sdss_s(i)=lum_int_sdss(i)
END DO

DO i=1,n_sdss_m
lumcor_sdss_m(i)=lum_int_sdss(i+n_sdss_s)
END DO

DO i=1,n_2mrs_s
lumcor_2mrs_s(i)=lum_int_2mrs(i)
END DO

DO i=1,n_2mrs_m
lumcor_2mrs_m(i)=lum_int_2mrs(i+n_2mrs_s)
END DO


DO i=1,n_2mrs_s
lum_2mrs_s(i)=(10.D0**lum_2mrs_s(i))/lumcor_2mrs_s(i)
lum_2mrs_s(i)=LOG10(lum_2mrs_s(i))
END DO

DO i=1,n_sdss_s
lum_sdss_s(i)=(10.D0**lum_sdss_s(i))/lumcor_sdss_s(i)
lum_sdss_s(i)=LOG10(lum_sdss_s(i))
END DO

DO i=1,n_2mrs_m
lum_2mrs_m(i)=(10.D0**lum_2mrs_m(i))/lumcor_2mrs_m(i)
lum_2mrs_m(i)=LOG10(lum_2mrs_m(i))
END DO

DO i=1,n_sdss_m
lum_sdss_m(i)=(10.D0**lum_sdss_m(i))/lumcor_sdss_m(i)
lum_sdss_m(i)=LOG10(lum_sdss_m(i))
END DO


WRITE(*,*) 'luminosity rescaled'


DO i=1,n_2mrs_m
mdyn_2mrs(i)=3.D0*((10.D0**sigma_2mrs_m(i))**2)*(10.D0**rad_2mrs_m(i))/G
END DO


DO i=1,n_sdss_m
mdyn_sdss(i)=3.D0*((10.D0**sigma_sdss_m(i))**2)*(10.D0**rad_sdss_m(i))/G
END DO



DO i=1,n_2mrs_m
mdyn_2mrs(i)=LOG10(mdyn_2mrs(i))
END DO


DO i=1,n_sdss_m
mdyn_sdss(i)=LOG10(mdyn_sdss(i))
END DO

WRITE(*,*) 'dynamical masses caluclated'



  OPEN(61,file='mass/newdata_2mrs_s.txt')
DO i=1,n_2mrs_s
WRITE(61,*) mass_2mrs_s(i),lum_2mrs_s(i),dist_2mrs_s(i)!,colour_2mrs_s(i) 
END DO
 CLOSE(61)
  
  OPEN(61,file='mass/newdata_sdss_s.txt')
DO i=1,n_sdss_s
WRITE(61,*) mass_sdss_s(i),lum_sdss_s(i),dist_sdss_s(i)!,colour_sdss_s(i)
END DO
 CLOSE(61)
  
  OPEN(61,file='mass/newdata_2mrs_m.txt')
DO i=1,n_2mrs_m
WRITE(61,*) mass_2mrs_m(i),lum_2mrs_m(i),dist_2mrs_m(i),mdyn_2mrs(i),nvis_2mrs_m(i)!,colour_2mrs_m(i)
END DO
 CLOSE(61)

  OPEN(61,file='mass/newdata_sdss_m.txt')
DO i=1,n_sdss_m
WRITE(61,*) mass_sdss_m(i),lum_sdss_m(i),dist_sdss_m(i),mdyn_sdss(i),nvis_sdss_m(i)!,colour_sdss_m(i)
END DO
 CLOSE(61)


 n_2mrs_mh=0
 n_sdss_mh=0
 n_2mrs_ml=0
 n_sdss_ml=0
 
 
   OPEN(61,file='mass/newdata_2mrs_m_low.txt')
DO i=1,n_2mrs_m
IF ((10**(nvis_2mrs_m(i)))<4.5) THEN
WRITE(61,*) mass_2mrs_m(i),lum_2mrs_m(i),dist_2mrs_m(i),mdyn_2mrs(i),nvis_2mrs_m(i)
n_2mrs_ml=n_2mrs_ml+1
END IF
END DO
 CLOSE(61)

  OPEN(61,file='mass/newdata_sdss_m_low.txt')
DO i=1,n_sdss_m
IF ((10**(nvis_sdss_m(i)))<4.5) THEN
WRITE(61,*) mass_sdss_m(i),lum_sdss_m(i),dist_sdss_m(i),mdyn_sdss(i),nvis_sdss_m(i)
n_sdss_ml=n_sdss_ml+1
END IF
END DO
 CLOSE(61)

  OPEN(61,file='mass/newdata_2mrs_m_high.txt')
DO i=1,n_2mrs_m
IF ((10**(nvis_2mrs_m(i)))>4.5) THEN
WRITE(61,*) mass_2mrs_m(i),lum_2mrs_m(i),dist_2mrs_m(i),mdyn_2mrs(i),nvis_2mrs_m(i)
n_2mrs_mh=n_2mrs_mh+1
END IF
END DO
 CLOSE(61)

  OPEN(61,file='mass/newdata_sdss_m_high.txt')
DO i=1,n_sdss_m
IF ((10**(nvis_sdss_m(i)))>4.5) THEN
WRITE(61,*) mass_sdss_m(i),lum_sdss_m(i),dist_sdss_m(i),mdyn_sdss(i),nvis_sdss_m(i)
n_sdss_mh=n_sdss_mh+1
END IF
END DO
 CLOSE(61)


 
 WRITE(*,*) 'data saved'

 

allocate(datalist_sdss_s(1:(dim_s+1),1:n_sdss_s))
allocate(datalist_2mrs_s(1:(dim_s+1),1:n_2mrs_s))



allocate(deviation_sdss_s(1:n_sdss_s))
allocate(deviation_2mrs_s(1:n_2mrs_s))
allocate(otherside_sdss_s(1:n_sdss_s))
allocate(otherside_2mrs_s(1:n_2mrs_s))


allocate(datalist_sdss_mh(1:(dim_mh+1),1:n_sdss_mh))
allocate(datalist_2mrs_mh(1:(dim_mh+1),1:n_2mrs_mh))
allocate(datalist_sdss_ml(1:(dim_ml+1),1:n_sdss_ml))
allocate(datalist_2mrs_ml(1:(dim_ml+1),1:n_2mrs_ml))

allocate(deviation_sdss_mh(1:n_sdss_mh))
allocate(deviation_2mrs_mh(1:n_2mrs_mh))
allocate(otherside_sdss_mh(1:n_sdss_mh))
allocate(otherside_2mrs_mh(1:n_2mrs_mh))

allocate(deviation_sdss_ml(1:n_sdss_ml))
allocate(deviation_2mrs_ml(1:n_2mrs_ml))
allocate(otherside_sdss_ml(1:n_sdss_ml))
allocate(otherside_2mrs_ml(1:n_2mrs_ml))
 
 
 
  WRITE(*,*) 'arrays allocated'

 

 
 
 

 DO i=1,n_sdss_s
 datalist_sdss_s(1,i)=mass_sdss_s(i)
 datalist_sdss_s(2,i)=lum_sdss_s(i)
 datalist_sdss_s(3,i)=lum_sdss_s(i)**2
 datalist_sdss_s(4,i)=lum_sdss_s(i)**3
 datalist_sdss_s(5,i)=dist_sdss_s(i)
 datalist_sdss_s(6,i)=dist_sdss_s(i)**2
 datalist_sdss_s(7,i)=dist_sdss_s(i)**3
 datalist_sdss_s(8,i)=1.D0
 END DO
 
 DO i=1,n_2mrs_s
 datalist_2mrs_s(1,i)=mass_2mrs_s(i)
 datalist_2mrs_s(2,i)=lum_2mrs_s(i) 
 datalist_2mrs_s(3,i)=lum_2mrs_s(i)**2
 datalist_2mrs_s(4,i)=lum_2mrs_s(i)**3
 datalist_2mrs_s(5,i)=dist_2mrs_s(i)
 datalist_2mrs_s(6,i)=dist_2mrs_s(i)**2
 datalist_2mrs_s(7,i)=dist_2mrs_s(i)**3
 datalist_2mrs_s(8,i)=1.D0
 END DO
 

  
 ii=0
 DO i=1,n_sdss_m
 IF ((10**(nvis_sdss_m(i)))<4.5) THEN
ii=ii+1
 datalist_sdss_ml(1,ii)=mass_sdss_m(i)
 datalist_sdss_ml(2,ii)=lum_sdss_m(i)
 datalist_sdss_ml(3,ii)=lum_sdss_m(i)**2
 datalist_sdss_ml(4,ii)=lum_sdss_m(i)**3
 datalist_sdss_ml(5,ii)=lum_sdss_m(i)**4
 datalist_sdss_ml(6,ii)=dist_sdss_m(i)
 datalist_sdss_ml(7,ii)=dist_sdss_m(i)**2
 datalist_sdss_ml(8,ii)=dist_sdss_m(i)**3
 datalist_sdss_ml(9,ii)=mdyn_sdss(i)
 datalist_sdss_ml(10,ii)=mdyn_sdss(i)**2
 datalist_sdss_ml(11,ii)=mdyn_sdss(i)**3
 datalist_sdss_ml(12,ii)=1.D0
 END IF
 END DO

 
 ii=0
 DO i=1,n_2mrs_m
 IF ((10**(nvis_2mrs_m(i)))<4.5) THEN
ii=ii+1
 datalist_2mrs_ml(1,ii)=mass_2mrs_m(i)
 datalist_2mrs_ml(2,ii)=lum_2mrs_m(i)
 datalist_2mrs_ml(3,ii)=lum_2mrs_m(i)**2
 datalist_2mrs_ml(4,ii)=lum_2mrs_m(i)**3
 datalist_2mrs_ml(5,ii)=lum_2mrs_m(i)**4
 datalist_2mrs_ml(6,ii)=dist_2mrs_m(i)
 datalist_2mrs_ml(7,ii)=dist_2mrs_m(i)**2
 datalist_2mrs_ml(8,ii)=dist_2mrs_m(i)**3
 datalist_2mrs_ml(9,ii)=mdyn_2mrs(i)
 datalist_2mrs_ml(10,ii)=mdyn_2mrs(i)**2
 datalist_2mrs_ml(11,ii)=mdyn_2mrs(i)**3
 datalist_2mrs_ml(12,ii)=1.D0
 END IF
 END DO
 
 
 
 ii=0
 DO i=1,n_sdss_m
 IF ((10**(nvis_sdss_m(i)))>4.5) THEN
ii=ii+1
 datalist_sdss_mh(1,ii)=mass_sdss_m(i)
 datalist_sdss_mh(2,ii)=lum_sdss_m(i)
 datalist_sdss_mh(3,ii)=lum_sdss_m(i)**2
 datalist_sdss_mh(4,ii)=lum_sdss_m(i)**3
 datalist_sdss_mh(5,ii)=lum_sdss_m(i)**4
 datalist_sdss_mh(6,ii)=lum_sdss_m(i)**5
 datalist_sdss_mh(7,ii)=dist_sdss_m(i)
 datalist_sdss_mh(8,ii)=dist_sdss_m(i)**2
 datalist_sdss_mh(9,ii)=dist_sdss_m(i)**3
 datalist_sdss_mh(10,ii)=mdyn_sdss(i)
 datalist_sdss_mh(11,ii)=mdyn_sdss(i)**2
 datalist_sdss_mh(12,ii)=mdyn_sdss(i)**3
 datalist_sdss_mh(13,ii)=nvis_sdss_m(i)
 datalist_sdss_mh(14,ii)=1.D0
 END IF
 END DO

 
 ii=0
 DO i=1,n_2mrs_m
 IF ((10**(nvis_2mrs_m(i)))>4.5) THEN
ii=ii+1
 datalist_2mrs_mh(1,ii)=mass_2mrs_m(i)
 datalist_2mrs_mh(2,ii)=lum_2mrs_m(i)
 datalist_2mrs_mh(3,ii)=lum_2mrs_m(i)**2
 datalist_2mrs_mh(4,ii)=lum_2mrs_m(i)**3
 datalist_2mrs_mh(5,ii)=lum_2mrs_m(i)**4
 datalist_2mrs_mh(6,ii)=lum_2mrs_m(i)**5
 datalist_2mrs_mh(7,ii)=dist_2mrs_m(i)
 datalist_2mrs_mh(8,ii)=dist_2mrs_m(i)**2
 datalist_2mrs_mh(9,ii)=dist_2mrs_m(i)**3
 datalist_2mrs_mh(10,ii)=mdyn_2mrs(i)
 datalist_2mrs_mh(11,ii)=mdyn_2mrs(i)**2
 datalist_2mrs_mh(12,ii)=mdyn_2mrs(i)**3
 datalist_2mrs_mh(13,ii)=nvis_2mrs_m(i) 
 datalist_2mrs_mh(14,ii)=1.D0
 END IF
 END DO
  

  
  
    WRITE(*,*) 'data assigned'

 

  
    
  OPEN(61,file='mass/weightedmfit_sdss_s.txt')
READ(61,*) rms_sdss_s
DO i=1,dim_s
READ(61,*) solution_sdss_s(i),errorbar_sdss_s(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/weightedmfit_2mrs_s.txt')
READ(61,*) rms_2mrs_s
DO i=1,dim_s
READ(61,*) solution_2mrs_s(i),errorbar_2mrs_s(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/weightedmfit_sdss_ml.txt')
READ(61,*) rms_sdss_ml
DO i=1,dim_ml
READ(61,*) solution_sdss_ml(i),errorbar_sdss_ml(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/weightedmfit_2mrs_ml.txt')
READ(61,*) rms_2mrs_ml
DO i=1,dim_ml
READ(61,*) solution_2mrs_ml(i),errorbar_2mrs_ml(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/weightedmfit_sdss_mh.txt')
READ(61,*) rms_sdss_mh
DO i=1,dim_mh
READ(61,*) solution_sdss_mh(i),errorbar_sdss_mh(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/weightedmfit_2mrs_mh.txt')
READ(61,*) rms_2mrs_mh
DO i=1,dim_mh
READ(61,*) solution_2mrs_mh(i),errorbar_2mrs_mh(i)
END DO
 CLOSE(61)
  
  
  
  
      WRITE(*,*) 'fit data loaded'

 
  
  
  

DO i=1,n_sdss_s
otherside_sdss_s(i)=0.D0
DO ii=2,(dim_s+1)
otherside_sdss_s(i)=otherside_sdss_s(i)+datalist_sdss_s(ii,i)*solution_sdss_s(ii-1)
END DO
END DO

DO i=1,n_2mrs_s
otherside_2mrs_s(i)=0.D0
DO ii=2,(dim_s+1)
otherside_2mrs_s(i)=otherside_2mrs_s(i)+datalist_2mrs_s(ii,i)*solution_2mrs_s(ii-1)
END DO
END DO
! WRITE(*,*) otherside_2mrs_s

DO i=1,n_sdss_ml
otherside_sdss_ml(i)=0.D0
DO ii=2,(dim_ml+1)
otherside_sdss_ml(i)=otherside_sdss_ml(i)+datalist_sdss_ml(ii,i)*solution_sdss_ml(ii-1)
END DO
END DO

DO i=1,n_2mrs_ml
otherside_2mrs_ml(i)=0.D0
DO ii=2,(dim_ml+1)
otherside_2mrs_ml(i)=otherside_2mrs_ml(i)+datalist_2mrs_ml(ii,i)*solution_2mrs_ml(ii-1)
END DO
END DO

DO i=1,n_sdss_mh
otherside_sdss_mh(i)=0.D0
DO ii=2,(dim_mh+1)
otherside_sdss_mh(i)=otherside_sdss_mh(i)+datalist_sdss_mh(ii,i)*solution_sdss_mh(ii-1)
END DO
END DO

DO i=1,n_2mrs_mh
otherside_2mrs_mh(i)=0.D0
DO ii=2,(dim_mh+1)
otherside_2mrs_mh(i)=otherside_2mrs_mh(i)+datalist_2mrs_mh(ii,i)*solution_2mrs_mh(ii-1)
END DO
END DO





DO i=1,n_sdss_s
deviation_sdss_s(i)=otherside_sdss_s(i)-datalist_sdss_s(1,i)
END DO
DO i=1,n_2mrs_s
deviation_2mrs_s(i)=otherside_2mrs_s(i)-datalist_2mrs_s(1,i)
END DO
DO i=1,n_sdss_ml
deviation_sdss_ml(i)=otherside_sdss_ml(i)-datalist_sdss_ml(1,i)
END DO
DO i=1,n_2mrs_ml
deviation_2mrs_ml(i)=otherside_2mrs_ml(i)-datalist_2mrs_ml(1,i)
END DO
DO i=1,n_sdss_mh
deviation_sdss_mh(i)=otherside_sdss_mh(i)-datalist_sdss_mh(1,i)
END DO
DO i=1,n_2mrs_mh
deviation_2mrs_mh(i)=otherside_2mrs_mh(i)-datalist_2mrs_mh(1,i)
END DO


WRITE(*,*) 'residuals caluclated'

mass_fit_sum_sdss_s=0.D0
mass_fit_sum_2mrs_s=0.D0
mass_orig_sum_sdss_s=0.D0
mass_orig_sum_2mrs_s=0.D0

mass_fit_sum_sdss_ml=0.D0
mass_fit_sum_2mrs_ml=0.D0
mass_orig_sum_sdss_ml=0.D0
mass_orig_sum_2mrs_ml=0.D0

mass_fit_sum_sdss_mh=0.D0
mass_fit_sum_2mrs_mh=0.D0
mass_orig_sum_sdss_mh=0.D0
mass_orig_sum_2mrs_mh=0.D0




DO i=1,n_sdss_s
mass_fit_sum_sdss_s=mass_fit_sum_sdss_s+(10.D0**otherside_sdss_s(i))
mass_orig_sum_sdss_s=mass_orig_sum_sdss_s+(10.D0**datalist_sdss_s(1,i))
END DO

DO i=1,n_sdss_ml
mass_fit_sum_sdss_ml=mass_fit_sum_sdss_ml+(10.D0**otherside_sdss_ml(i))
mass_orig_sum_sdss_ml=mass_orig_sum_sdss_ml+(10.D0**datalist_sdss_ml(1,i))
END DO

DO i=1,n_sdss_mh
mass_fit_sum_sdss_mh=mass_fit_sum_sdss_mh+(10.D0**otherside_sdss_mh(i))
mass_orig_sum_sdss_mh=mass_orig_sum_sdss_mh+(10.D0**datalist_sdss_mh(1,i))
END DO

DO i=1,n_2mrs_s
mass_fit_sum_2mrs_s=mass_fit_sum_2mrs_s+(10.D0**otherside_2mrs_s(i))
mass_orig_sum_2mrs_s=mass_orig_sum_2mrs_s+(10.D0**datalist_2mrs_s(1,i))
END DO

DO i=1,n_2mrs_ml
mass_fit_sum_2mrs_ml=mass_fit_sum_2mrs_ml+(10.D0**otherside_2mrs_ml(i))
mass_orig_sum_2mrs_ml=mass_orig_sum_2mrs_ml+(10.D0**datalist_2mrs_ml(1,i))
END DO

DO i=1,n_2mrs_mh
mass_fit_sum_2mrs_mh=mass_fit_sum_2mrs_mh+(10.D0**otherside_2mrs_mh(i))
mass_orig_sum_2mrs_mh=mass_orig_sum_2mrs_mh+(10.D0**datalist_2mrs_mh(1,i))
END DO


mass_fit_sum_sdss_tot=mass_fit_sum_sdss_s+mass_fit_sum_sdss_ml+mass_fit_sum_sdss_mh
mass_fit_sum_2mrs_tot=mass_fit_sum_2mrs_s+mass_fit_sum_2mrs_ml+mass_fit_sum_2mrs_mh
mass_orig_sum_sdss_tot=mass_orig_sum_sdss_s+mass_orig_sum_sdss_ml+mass_orig_sum_sdss_mh
mass_orig_sum_2mrs_tot=mass_orig_sum_2mrs_s+mass_orig_sum_2mrs_ml+mass_orig_sum_2mrs_mh



OPEN(61,file='mass/weightedsums.txt')
WRITE(61,*) mass_fit_sum_sdss_s,mass_orig_sum_sdss_s
WRITE(61,*) mass_fit_sum_sdss_ml,mass_orig_sum_sdss_ml
WRITE(61,*) mass_fit_sum_sdss_mh,mass_orig_sum_sdss_mh
WRITE(61,*) mass_fit_sum_2mrs_s,mass_orig_sum_2mrs_s
WRITE(61,*) mass_fit_sum_2mrs_ml,mass_orig_sum_2mrs_ml
WRITE(61,*) mass_fit_sum_2mrs_mh,mass_orig_sum_2mrs_mh
WRITE(61,*)
WRITE(61,*) mass_fit_sum_sdss_tot,mass_orig_sum_sdss_tot
WRITE(61,*) mass_fit_sum_2mrs_tot,mass_orig_sum_2mrs_tot
WRITE(61,*)
WRITE(61,*) (mass_fit_sum_sdss_tot/m_expected),(mass_orig_sum_sdss_tot/m_expected)
WRITE(61,*) (mass_fit_sum_2mrs_tot/m_expected),(mass_orig_sum_2mrs_tot/m_expected)
 CLOSE(61)





WRITE(*,*) 'masses summed up'




DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_ml
help_x=NINT((datalist_2mrs_ml(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_ml.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_mh
help_x=NINT((datalist_2mrs_mh(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_mh.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_ml
help_x=NINT((datalist_sdss_ml(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_ml.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
  
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_mh
help_x=NINT((datalist_sdss_mh(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_mh.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 WRITE(*,*) 'mass residuals written'
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_ml
help_x=NINT((datalist_2mrs_ml(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_ml.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_mh
help_x=NINT((datalist_2mrs_mh(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_mh.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)



DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_s.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_ml
help_x=NINT((datalist_sdss_ml(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_ml.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_mh
help_x=NINT((datalist_sdss_mh(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_mh.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_s.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

 
 
 
 
  WRITE(*,*) 'luminosity residuals written'
 
 
 
 
 
 
 
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_ml
help_x=NINT((datalist_2mrs_ml(6,i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_ml.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_mh
help_x=NINT((datalist_2mrs_mh(7,i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_mh.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(6,i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_ml
help_x=NINT((datalist_sdss_ml(6,i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_ml.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_mh
help_x=NINT((datalist_sdss_mh(7,i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_mh.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)



DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(6,i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

 
 
 
  WRITE(*,*) 'distance residuals written'
 
 
 

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_ml
help_x=NINT((datalist_2mrs_ml(9,i)-1.4D0)*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_2mrs_ml.txt')
DO i=0,100
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.4D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_mh
help_x=NINT((datalist_2mrs_mh(10,i)-1.4D0)*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_2mrs_mh.txt')
DO i=0,100
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.4D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_ml
help_x=NINT((datalist_sdss_ml(9,i)-1.4D0)*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_sdss_ml.txt')
DO i=0,100
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.4D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_mh
help_x=NINT((datalist_sdss_mh(10,i)-1.4D0)*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_sdss_mh.txt')
DO i=0,100
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.4D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 

 
 
 
  WRITE(*,*) 'dynamical mass residuals written'
 

 
  

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_mh
help_x=NINT((datalist_2mrs_mh(13,i))*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<100.D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+0.05D0
END IF
END IF
END IF
END DO

OPEN(61,file='mass/res_nvis_2mrs_mh.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_mh
help_x=NINT((datalist_sdss_mh(13,i))*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<100.D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+0.05D0
END IF
END IF
END IF
END DO

OPEN(61,file='mass/res_nvis_sdss_mh.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 
 
 
 
  WRITE(*,*) 'number residuals written'
 
 
 
 
 
 
 
 
 
 
 
 
 OPEN(61,file='mass/all_residuals_sdss_s.txt')
DO i=1,n_sdss_s
WRITE(61,*) deviation_sdss_s(i),datalist_sdss_s(2:dim_s,i)
END DO
 CLOSE(61)

 
  OPEN(61,file='mass/all_residuals_sdss_ml.txt')
DO i=1,n_sdss_ml
WRITE(61,*) deviation_sdss_ml(i),datalist_sdss_ml(2:dim_ml,i)
END DO
 CLOSE(61)
 
 
   OPEN(61,file='mass/all_residuals_sdss_mh.txt')
DO i=1,n_sdss_mh
WRITE(61,*) deviation_sdss_mh(i),datalist_sdss_mh(2:dim_mh,i)
END DO
 CLOSE(61)
  
 
 OPEN(61,file='mass/all_residuals_2mrs_s.txt')
DO i=1,n_2mrs_s
WRITE(61,*) deviation_2mrs_s(i),datalist_2mrs_s(2:dim_s,i)
END DO
 CLOSE(61)

 
  OPEN(61,file='mass/all_residuals_2mrs_ml.txt')
DO i=1,n_2mrs_ml
WRITE(61,*) deviation_2mrs_ml(i),datalist_2mrs_ml(2:dim_ml,i)
END DO
 CLOSE(61)
 
 
   OPEN(61,file='mass/all_residuals_2mrs_mh.txt')
DO i=1,n_2mrs_mh
WRITE(61,*) deviation_2mrs_mh(i),datalist_2mrs_mh(2:dim_mh,i)
END DO
 CLOSE(61)
 
 
 
 
 WRITE(*,*) 'all residuals written'



END PROGRAM 








