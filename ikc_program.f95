PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,n_2mass,n_sdss,i,ii,iii,iiii,active_count,hcount,nedge,n_gal_max,n_unified,n_halo_max,n_hunified
integer ::  n_qso,n_mock_final,n_sdss_used,n_2mass_used,n_unified2,n_mock_final2,uuu,n_mockcat,n_dummy
integer, dimension(1:6) :: n_gal,n_halo
double precision :: dummy_id,divisor,helpflip,av_sigma,sigma_sigma,binhelp,zmax,zmin,disthelp1,disthelp2,cmv
character(19), allocatable :: sdss_objID(:)
double precision, allocatable ::sdss_ra(:),sdss_dec(:),sdss_z(:),sdss_gl(:),sdss_gb(:),z_cor(:)
double precision, allocatable ::  cModelMag_g(:),cModelMag_r(:)
double precision, allocatable ::  extinction_g(:),extinction_r(:)
double precision, allocatable :: cModelMagErr_g(:),cModelMagErr_r(:)
character(200) :: read_string
character(20) :: helpstring
 character(20) :: number_mockset
double precision, allocatable :: twomass_ID1(:),twomass_ID2(:),twomass_ra(:),twomass_dec(:)
double precision, allocatable :: twomass_gl(:),twomass_gb(:)
double precision, allocatable :: Ktmag(:),Jtmag(:),Kmagerr(:),Jmagerr(:)
double precision, allocatable :: E_BV(:),twomass_cz(:),twomass_ecz(:) 

double precision :: sdss_cover,twomass_cover,mock_cover,ratio,help_ratio,qso_ratio,av_red_2mass,av_J_err,av_Ks_err



double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp

double precision, dimension(1:6) :: mtot,m_ratio,partsum

logical, allocatable :: twomass_active(:),sdss_active(:),all_active(:)
double precision :: v_cmb,l_cmb,b_cmb,z_cmb,gal_b,gal_l,z_cmb_x,z_cmb_y,z_cmb_z,z_gal_x,z_gal_y,z_gal_z,arcsec
double precision :: limitmag_g,limitmag_r,limitmag_i,limitmag_z,limitmag_J,limitmag_H,limitmag_Ks
double precision :: limitmag_SDSS_official,limitmag_2MRS_official,completeness_2mrs
double precision, allocatable ::  helplist_J(:),helplist_H(:),helplist_Ks(:)

!K correction coefficients
double precision, dimension(0:5,0:3) :: K_coeff_r,K_coeff_i,K_coeff_z,K_coeff_J,K_coeff_Ks
double precision, dimension(0:7,0:3) :: K_coeff_g,K_coeff_H
double precision :: dummy_coeffient,d_J,e_J,f_J,d_Ks,e_Ks,f_Ks
double precision, allocatable :: K_g(:),K_r(:)

 double precision, allocatable :: magapp_g(:),magapp_r(:)
 double precision, allocatable :: K_J(:),K_Ks(:),z_2mass(:)
 double precision, allocatable :: magapp_J(:),magapp_Ks(:)
 double precision, allocatable :: distance_sdss(:),distance_2mass(:)
  double precision, allocatable :: loglumdist_sdss(:),loglumdist_2mass(:)
 double precision, allocatable :: magabs_g(:),magabs_r(:)
 double precision, allocatable :: magabs_J(:),magabs_Ks(:)
 double precision, allocatable :: xpos_sdss(:),ypos_sdss(:),zpos_sdss(:)
  double precision, allocatable :: xpos_2mass(:),ypos_2mass(:),zpos_2mass(:)
   double precision, allocatable :: magorig_g(:),magorig_r(:),magorig_J(:),magorig_Ks(:)
  
   double precision ::  fiber_coll_r,decollided_sampling,loss_coll,final_sampling,angular_sep
  
integer, dimension(0:500,0:400) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y

double precision :: z_dummy,color_dummy
double precision, allocatable ::  Kcor_g_dummy(:),Kcor_r_dummy(:),Kcor_J_dummy(:),Kcor_Ks_dummy(:)
double precision, allocatable ::  z_dummy_array(:),colour_gr_dummy(:),colour_JK_dummy(:)

double precision, dimension(0:110) :: bin_sdss,bin_2mass



OPEN(33,file='logfile.txt')


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme INPUT K-CORRECTION started'
WRITE(*,*) '============================================================'
WRITE(33,*) '============================================================'
WRITE(33,*) '    programme INPUT K-CORRECTION started'
WRITE(33,*) '============================================================'

! random seed
 CALL SYSTEM_CLOCK(hcount)
 hcount=hcount-INT(hcount/100000)*100000
 CALL srand(hcount) 


! define constants
PI=ACOS(-1.D0)
arcsec=2.D0/3600.D0

n_mockcat=8

n_qso=1889
sdss_cover=9274.D0/41253.D0
twomass_cover=0.91D0
mock_cover=1.D0/8.D0

fiber_coll_r=55.D0/3600.D0
decollided_sampling=0.99D0
loss_coll=0.D0
final_sampling=0.92D0
 completeness_2mrs=0.976D0

limitmag_SDSS_official=17.77D0
limitmag_2MRS_official=11.75D0

Omega_m=0.272D0
Omega_l=0.728D0

q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
H0=70.4D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)
 cV_mill=(350.D6/h)**3
m_expected=cV_mill*rho_crit*Omega_m
part_exp=(2160.D0**3)*((350.D0/500.D0)**3)







! get length of file
OPEN(50,file='2MRS.txt')
io_err=0
n_2mass=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mass=n_2mass+1
END DO
 CLOSE(50)
n_2mass=n_2mass-3 ! remove header

! get length of file
OPEN(50,file='SDSS_DR10.txt')
io_err=0
n_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss=n_sdss+1
END DO
 CLOSE(50)
n_sdss=n_sdss-2 ! remove header




allocate(sdss_objID(1:n_sdss))
allocate(sdss_ra(1:n_sdss))
allocate(sdss_dec(1:n_sdss))
allocate(sdss_gl(1:n_sdss))
allocate(sdss_gb(1:n_sdss))
allocate(sdss_z(1:n_sdss))
allocate(cModelMag_g(1:n_sdss))
allocate(cModelMag_r(1:n_sdss))
allocate(cModelMagErr_g(1:n_sdss))
allocate(cModelMagErr_r(1:n_sdss))
allocate(extinction_g(1:n_sdss))
allocate(extinction_r(1:n_sdss))
allocate(sdss_active(1:n_sdss))

allocate(twomass_ID1(1:n_2mass))
allocate(twomass_ID2(1:n_2mass))
allocate(twomass_ra(1:n_2mass))
allocate(twomass_dec(1:n_2mass))
allocate(twomass_gl(1:n_2mass))
allocate(twomass_gb(1:n_2mass))
allocate(Ktmag(1:n_2mass))
allocate(Jtmag(1:n_2mass))
allocate(Kmagerr(1:n_2mass))
allocate(Jmagerr(1:n_2mass))
allocate(E_BV(1:n_2mass))
allocate(twomass_cz(1:n_2mass))
allocate(twomass_ecz(1:n_2mass))
allocate(twomass_active(1:n_2mass))


allocate(K_g(1:n_sdss))
allocate(K_r(1:n_sdss))
allocate(z_cor(1:n_sdss))
allocate(magapp_g(1:n_sdss))
allocate(magapp_r(1:n_sdss))
allocate(distance_sdss(1:n_sdss))
allocate(magabs_g(1:n_sdss))
allocate(magabs_r(1:n_sdss))
allocate(loglumdist_sdss(1:n_sdss))


allocate(K_J(1:n_2mass))
allocate(K_Ks(1:n_2mass))
allocate(z_2mass(1:n_2mass))
allocate(magapp_J(1:n_2mass))
allocate(magapp_Ks(1:n_2mass))
allocate(distance_2mass(1:n_2mass))
allocate(magabs_J(1:n_2mass))
allocate(magabs_Ks(1:n_2mass))
allocate(loglumdist_2mass(1:n_2mass))

allocate(helplist_J(1:n_2mass))
allocate(helplist_Ks(1:n_2mass))

allocate(xpos_sdss(1:n_sdss))
allocate(ypos_sdss(1:n_sdss))
allocate(zpos_sdss(1:n_sdss))
allocate(xpos_2mass(1:n_2mass))
allocate(ypos_2mass(1:n_2mass))
allocate(zpos_2mass(1:n_2mass))


allocate(magorig_g(1:n_sdss))
allocate(magorig_r(1:n_sdss))
allocate(magorig_J(1:n_2mass))
allocate(magorig_Ks(1:n_2mass))



! read file
OPEN(50,file='SDSS_DR10.txt')
READ(50,*)
DO i=1,n_sdss

READ(50,*) dummy_id,sdss_ra(i),sdss_dec(i),sdss_gb(i),sdss_gl(i),sdss_z(i),&
 cModelMag_g(i),cModelMag_r(i),dummy_id,cModelMagErr_g(i),cModelMagErr_r(i),dummy_id,&
 extinction_g(i),extinction_r(i),dummy_id

sdss_active(i)=.TRUE.

END DO

 CLOSE(50)


! read file
OPEN(50,file='SDSS_DR10.txt')
READ(50,*)
DO i=1,n_sdss 
READ(50,*) read_string
sdss_objID(i)=read_string(1:19)
END DO

 CLOSE(50)


WRITE(*,*) n_sdss,'galaxies from SDSS DR10'
WRITE(33,*) n_sdss,'galaxies from SDSS DR10'


! read file
OPEN(50,file='2MRS.txt')

READ(50,*)
READ(50,*)

DO i=1,n_2mass 
READ(50,*) twomass_ID1(i),twomass_ID2(i),twomass_ra(i),twomass_dec(i),&
twomass_gl(i),twomass_gb(i),Ktmag(i),Kmagerr(i),Jtmag(i),Jmagerr(i),&
E_BV(i),twomass_cz(i),twomass_ecz(i)

twomass_active(i)=.TRUE.

END DO

 CLOSE(50)
WRITE(*,*) n_2mass,'galaxies from 2MRS'
WRITE(33,*) n_2mass,'galaxies from 2MRS'



 











!prepare for calculation of correction for CMB motion
! values for the motion of the sun relativ to the CMB
v_cmb=369.D0
l_cmb=263.99D0
b_cmb=48.26D0

! calculate sun's motion in z space relativ to the CMB
z_cmb=v_cmb/light
gal_b=b_cmb*PI/180.D0
gal_l=l_cmb*PI/180.D0
z_cmb_x=z_cmb*COS(gal_b)*COS(gal_l)
z_cmb_y=z_cmb*COS(gal_b)*SIN(gal_l)
z_cmb_z=z_cmb*SIN(gal_b)


DO i=1,n_sdss
gal_b=sdss_gb(i)*PI/180.D0
gal_l=sdss_gl(i)*PI/180.D0
z_gal_x=sdss_z(i)*COS(gal_b)*COS(gal_l)
z_gal_y=sdss_z(i)*COS(gal_b)*SIN(gal_l)
z_gal_z=sdss_z(i)*SIN(gal_b)
z_gal_x=z_gal_x+z_cmb_x
z_gal_y=z_gal_y+z_cmb_y
z_gal_z=z_gal_z+z_cmb_z
z_cor(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


DO i=1,n_2mass
gal_b=twomass_gb(i)*PI/180.D0
gal_l=twomass_gl(i)*PI/180.D0
z_gal_x=twomass_cz(i)/light*COS(gal_b)*COS(gal_l)
z_gal_y=twomass_cz(i)/light*COS(gal_b)*SIN(gal_l)
z_gal_z=twomass_cz(i)/light*SIN(gal_b)
z_gal_x=z_gal_x+z_cmb_x
z_gal_y=z_gal_y+z_cmb_y
z_gal_z=z_gal_z+z_cmb_z
z_2mass(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO



!get parameters for K correction  
OPEN(52,file='K_correction_g.txt')
DO i=0,7
READ(52,*) K_coeff_g(i,0:3)
END DO
 CLOSE(52)
OPEN(52,file='K_correction_r.txt')
DO i=0,5
READ(52,*) K_coeff_r(i,0:3)
END DO
 CLOSE(52)


OPEN(52,file='K_correction_J.txt')
DO i=0,5
READ(52,*) K_coeff_J(i,0:3)
END DO
 CLOSE(52)

OPEN(52,file='K_correction_Ks.txt')
DO i=0,5
READ(52,*) K_coeff_Ks(i,0:3)
END DO
 CLOSE(52)

 OPEN(52,file='limit_magnitudes.txt')
READ(52,*) dummy_id,limitmag_g,limitmag_r,limitmag_i,dummy_id
 CLOSE(52)
 
 
 
 OPEN(52,file='2mass_sdss_transformation3.txt')
READ(52,*) 
READ(52,*) d_J,dummy_coeffient,d_Ks
READ(52,*) e_J,dummy_coeffient,e_Ks
READ(52,*) f_J,dummy_coeffient,f_Ks
 CLOSE(52)
 

 
 
 
 
 
 
 
 
!correct for galactic extinction
DO i=1,n_sdss
magapp_g(i)=cModelMag_g(i)-extinction_g(i)
magapp_r(i)=cModelMag_r(i)-extinction_r(i)

END DO


DO i=1,n_2mass
magapp_J(i)=Jtmag(i)
magapp_Ks(i)=Ktmag(i)
END DO


DO i=1,n_sdss
magorig_g(i)=magapp_g(i)
magorig_r(i)=magapp_r(i)
END DO

DO i=1,n_2mass
magorig_J(i)=magapp_J(i)
magorig_Ks(i)=magapp_Ks(i)
END DO


!K-correction
DO i=1,n_sdss

K_g(i)=0.D0
DO ii=0,7
DO iii=0,3
K_g(i)=K_g(i)+K_coeff_g(ii,iii)*(sdss_z(i)**ii)*((magapp_g(i)-magapp_r(i))**iii)
END DO
END DO

K_r(i)=0.D0
DO ii=0,5
DO iii=0,3
K_r(i)=K_r(i)+K_coeff_r(ii,iii)*(sdss_z(i)**ii)*((magapp_g(i)-magapp_r(i))**iii)
END DO
END DO



END DO 



!K-correction
DO i=1,n_2mass

K_J(i)=0.D0
DO ii=0,5
DO iii=0,3
K_J(i)=K_J(i)+K_coeff_J(ii,iii)*((twomass_cz(i)/light)**ii)*((magapp_J(i)-magapp_Ks(i))**iii)
END DO
END DO


K_Ks(i)=0.D0
DO ii=0,5
DO iii=0,3
K_Ks(i)=K_Ks(i)+K_coeff_Ks(ii,iii)*((twomass_cz(i)/light)**ii)*((magapp_J(i)-magapp_Ks(i))**iii)
END DO
END DO

! WRITE(*,*) K_Ks(i),K_J(i)

END DO


! corrected apparent magnitudes
DO i=1,n_sdss

magapp_g(i)=magapp_g(i)-K_g(i) 
magapp_r(i)=magapp_r(i)-K_r(i) 


END DO



! corrected apparent magnitudes
DO i=1,n_2mass

magapp_J(i)=magapp_J(i)-K_J(i) 

magapp_Ks(i)=magapp_Ks(i)-K_Ks(i) 

END DO


DO i=1,n_sdss
IF (z_cor(i)>0.D0) THEN

divisor=(SQRT(1.D0+2.D0*q0*z_cor(i))+1.D0+q0*z_cor(i))
!luminosity distance
distance_sdss(i)=light/H0*z_cor(i)*(1.D0+((z_cor(i)*(1.D0-q0))/divisor))
!distance in kpc
distance_sdss(i)=distance_sdss(i)*1000.D0

ELSE
sdss_active(i)=.FALSE.
z_cor(i)=10000.D0
distance_sdss(i)=1000000.D0
END IF
END DO


DO i=1,n_sdss
IF (sdss_active(i)) THEN
loglumdist_sdss(i)=LOG10(distance_sdss(i)*1000.D0)
END IF
END DO

DO i=1,n_sdss
IF (sdss_active(i)) THEN
magabs_g(i)=magapp_g(i)-5.D0*LOG10(distance_sdss(i)*1000.D0)+5.D0
magabs_r(i)=magapp_r(i)-5.D0*LOG10(distance_sdss(i)*1000.D0)+5.D0

END IF
END DO


DO i=1,n_sdss
IF (sdss_active(i)) THEN
IF (magabs_r(i)>-15.D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO

DO i=1,n_sdss
IF (sdss_active(i)) THEN
IF (z_cor(i)>0.11D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO

DO i=1,n_sdss
IF (sdss_active(i)) THEN
IF (z_cor(i)<0.D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO

 DO i=1,n_sdss
 IF (sdss_active(i)) THEN
IF (cModelMagErr_g(i)>1.D0) THEN
sdss_active(i)=.FALSE.
END IF
IF (cModelMagErr_g(i)<0.D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO

 DO i=1,n_sdss
 IF (sdss_active(i)) THEN
IF (cModelMagErr_r(i)>1.D0) THEN
sdss_active(i)=.FALSE.
END IF
IF (cModelMagErr_r(i)<0.D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO


DO i=1,n_2mass
IF (twomass_cz(i)<0.D0) THEN
twomass_active(i)=.FALSE.
END IF
END DO

DO i=1,n_2mass
IF (z_2mass(i)>0.D0) THEN

divisor=(SQRT(1.D0+2.D0*q0*z_2mass(i))+1.D0+q0*z_2mass(i))
!luminosity distance
distance_2mass(i)=light/H0*z_2mass(i)*(1.D0+((z_2mass(i)*(1.D0-q0))/divisor))
!distance in kpc
distance_2mass(i)=distance_2mass(i)*1000.D0


ELSE
twomass_active(i)=.FALSE.
z_2mass(i)=10000.D0
distance_2mass(i)=1000000.D0
END IF
END DO

DO i=1,n_2mass
IF (twomass_active(i)) THEN
loglumdist_2mass(i)=LOG10(distance_2mass(i)*1000.D0)
END IF
END DO

DO i=1,n_2mass
IF (twomass_active(i)) THEN
magabs_J(i)=magapp_J(i)-5.D0*LOG10(distance_2mass(i)*1000.D0)+5.D0
magabs_Ks(i)=magapp_Ks(i)-5.D0*LOG10(distance_2mass(i)*1000.D0)+5.D0
END IF
END DO

DO i=1,n_sdss
distance_sdss(i)=distance_sdss(i)*((1.D0+z_cor(i))**(-1.D0))
END DO

DO i=1,n_2mass
distance_2mass(i)=distance_2mass(i)*((1.D0+z_2mass(i))**(-1.D0))
END DO


DO i=1,n_sdss




IF (magapp_r(i)>(limitmag_SDSS_official+0.5D0)) THEN
sdss_active(i)=.FALSE.
END IF


END DO




DO i=1,n_2mass



IF (magapp_Ks(i)>(limitmag_2MRS_official+0.5D0)) THEN
twomass_active(i)=.FALSE.
END IF

                 

END DO





n_dummy=150*250+1


allocate(Kcor_g_dummy(1:n_dummy))
allocate(Kcor_r_dummy(1:n_dummy))
allocate(Kcor_J_dummy(1:n_dummy))
allocate(Kcor_Ks_dummy(1:n_dummy))
allocate(z_dummy_array(1:n_dummy))
allocate(colour_gr_dummy(1:n_dummy))
allocate(colour_JK_dummy(1:n_dummy))


DO i=1,150
DO ii=1,250
z_dummy=DBLE(i-1)/1000.D0
 color_dummy=DBLE(ii-50)/100.D0
 
z_dummy_array(i+(ii-1)*150)=z_dummy
 colour_gr_dummy(i+(ii-1)*150)=color_dummy
 colour_JK_dummy(i+(ii-1)*150)=color_dummy 
 
END DO
END DO


!K-correction
DO i=1,(n_dummy-1)

Kcor_g_dummy(i)=0.D0
DO ii=0,7
DO iii=0,3
Kcor_g_dummy(i)=Kcor_g_dummy(i)+K_coeff_g(ii,iii)*(z_dummy_array(i)**ii)*(colour_gr_dummy(i)**iii)
END DO
END DO

Kcor_r_dummy(i)=0.D0
DO ii=0,5
DO iii=0,3
Kcor_r_dummy(i)=Kcor_r_dummy(i)+K_coeff_r(ii,iii)*(z_dummy_array(i)**ii)*(colour_gr_dummy(i)**iii)
END DO
END DO


Kcor_J_dummy(i)=0.D0
DO ii=0,5
DO iii=0,3
Kcor_J_dummy(i)=Kcor_J_dummy(i)+K_coeff_J(ii,iii)*(z_dummy_array(i)**ii)*(colour_JK_dummy(i)**iii)
END DO
END DO


Kcor_Ks_dummy(i)=0.D0
DO ii=0,5
DO iii=0,3
Kcor_Ks_dummy(i)=Kcor_Ks_dummy(i)+K_coeff_Ks(ii,iii)*(z_dummy_array(i)**ii)*(colour_JK_dummy(i)**iii)
END DO
END DO

END DO

WRITE(*,*) 'writing output for',n_dummy,'artificial datapoints'
WRITE(33,*) 'writing output for',n_dummy,'artificial datapoints'





DO i=1,n_sdss
IF (sdss_active(i)) THEN
IF ((magapp_g(i)-magapp_r(i))>3.D0) THEN
sdss_active(i)=.FALSE.
END IF
IF ((magapp_g(i)-magapp_r(i))<-1.D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO


DO i=1,n_2mass
IF (twomass_active(i)) THEN
IF ((magapp_J(i)-magapp_Ks(i))>3.D0) THEN
twomass_active(i)=.FALSE.
END IF
IF ((magapp_J(i)-magapp_Ks(i))<-1.D0) THEN
twomass_active(i)=.FALSE.
END IF
END IF
END DO

DO i=1,n_sdss
IF (sdss_active(i)) THEN
IF ((magorig_g(i)-magorig_r(i))>3.D0) THEN
sdss_active(i)=.FALSE.
END IF
IF ((magorig_g(i)-magorig_r(i))<-1.D0) THEN
sdss_active(i)=.FALSE.
END IF
END IF
END DO


DO i=1,n_2mass
IF (twomass_active(i)) THEN
IF ((magorig_J(i)-magorig_Ks(i))>3.D0) THEN
twomass_active(i)=.FALSE.
END IF
IF ((magorig_J(i)-magorig_Ks(i))<-1.D0) THEN
twomass_active(i)=.FALSE.
END IF
END IF
END DO


active_count=0
DO i=1,n_sdss
IF (sdss_active(i)) THEN
active_count=active_count+1
END IF
END DO
n_sdss_used=active_count

WRITE(*,*) 'writing output for',n_sdss_used,'SDSS galaxies '
WRITE(33,*) 'writing output for',n_sdss_used,'SDSS galaxies '


active_count=0
DO i=1,n_2mass
IF (twomass_active(i)) THEN
active_count=active_count+1
END IF
END DO
n_2mass_used=active_count
WRITE(*,*) 'writing output for',n_2mass_used,'2MRS galaxies '
WRITE(33,*) 'writing output for',n_2mass_used,'2MRS galaxies '


OPEN(50,file='basis_K_cor_SDSS_g.txt')
DO i=1,n_sdss
IF (sdss_active(i)) THEN
WRITE(50,*) sdss_z(i),(magapp_g(i)-magapp_r(i)),K_g(i)
END IF
END DO
DO i=1,(n_dummy-1)
WRITE(50,*) z_dummy_array(i),colour_gr_dummy(i),Kcor_g_dummy(i)
END DO
CLOSE(50)

OPEN(50,file='basis_K_cor_SDSS_r.txt')
DO i=1,n_sdss
IF (sdss_active(i)) THEN
WRITE(50,*) sdss_z(i),(magapp_g(i)-magapp_r(i)),K_r(i)
END IF
END DO
DO i=1,(n_dummy-1)
WRITE(50,*) z_dummy_array(i),colour_gr_dummy(i),Kcor_r_dummy(i)
END DO
CLOSE(50)

OPEN(50,file='basis_K_cor_2MRS_J.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) (twomass_cz(i)/light),(magapp_J(i)-magapp_Ks(i)),K_J(i)
END IF
END DO
DO i=1,(n_dummy-1)
WRITE(50,*) z_dummy_array(i),colour_JK_dummy(i),Kcor_J_dummy(i)
END DO
CLOSE(50)

OPEN(50,file='basis_K_cor_2MRS_Ks.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) (twomass_cz(i)/light),(magapp_J(i)-magapp_Ks(i)),K_Ks(i)
END IF
END DO
DO i=1,(n_dummy-1)
WRITE(50,*) z_dummy_array(i),colour_JK_dummy(i),Kcor_Ks_dummy(i)
END DO
CLOSE(50)

 
 
 
 
 
 OPEN(50,file='colourshift_gr.txt')
DO i=1,n_sdss
IF (sdss_active(i)) THEN
WRITE(50,*) (magapp_g(i)-magapp_r(i)),((magorig_g(i)-magorig_r(i))-(magapp_g(i)-magapp_r(i)))
END IF
END DO
CLOSE(50) 
 
 
  OPEN(50,file='colourshift_JK.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) (magapp_J(i)-magapp_Ks(i)),((magorig_J(i)-magorig_Ks(i))-(magapp_J(i)-magapp_Ks(i)))
END IF
END DO
CLOSE(50) 
 


 
 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,n_sdss
IF (sdss_active(i)) THEN
help_x=NINT((magapp_g(i)-magapp_r(i)+0.5D0)*100.D0)
help_y=NINT(((magorig_g(i)-magorig_r(i))-(magapp_g(i)-magapp_r(i))+0.5D0)*100.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<100)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END IF


END DO

 DO i=0,500
DO ii=0,400
IF (mapmap(i,ii)>0) THEN
mapmap(i,ii)=mapmap(i,ii)+100
END IF
END DO
END DO

OPEN(61,file='colourshift_gr_map.txt')
DO i=0,250
DO ii=0,100
help_x2=(DBLE(i)/100.D0)-0.5D0
help_y2=(DBLE(ii)/100.D0)-0.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


  
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,n_2mass
IF (twomass_active(i)) THEN
help_x=NINT((magapp_J(i)-magapp_Ks(i)+0.5)*100.D0)
help_y=NINT(((magorig_J(i)-magorig_Ks(i))-(magapp_J(i)-magapp_Ks(i))+0.20D0)*250.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<100)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END IF

END DO


 DO i=0,500
DO ii=0,400
IF (mapmap(i,ii)>0) THEN
mapmap(i,ii)=mapmap(i,ii)+5
END IF
END DO
END DO

OPEN(61,file='colourshift_JK_map.txt')
DO i=0,250
DO ii=0,100
help_x2=(DBLE(i)/100.D0)-0.5D0
help_y2=(DBLE(ii)/250.D0)-0.2D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'
WRITE(33,*) '============================================================'
WRITE(33,*) '    programme complete'
WRITE(33,*) '============================================================'

CLOSE(33)


END PROGRAM