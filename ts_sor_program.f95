PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE

! include 'mpif.h'

real :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,n_2mass,n_sdss,i,ii,iii,iiii,active_count,hcount,nedge,n_gal_max,n_unified,n_halo_max,n_hunified
integer ::  n_qso,n_mock_final,n_sdss_used,n_2mass_used,n_unified2,n_mock_final2,uuu,n_mockcat
integer, dimension(1:6) :: n_gal,n_halo
real :: dummy_id,divisor,helpflip,av_sigma,sigma_sigma,binhelp,zmax,zmin,disthelp1,disthelp2,cmv,z_pec
character(19), allocatable :: sdss_objID(:)
real, allocatable ::sdss_ra(:),sdss_dec(:),sdss_z(:),sdss_gl(:),sdss_gb(:),z_cor(:)
real, allocatable ::  cModelMag_g(:),cModelMag_r(:)
real, allocatable ::  extinction_g(:),extinction_r(:)
real, allocatable :: cModelMagErr_g(:),cModelMagErr_r(:)
character(200) :: read_string
character(20) :: helpstring
 character(20) :: number_mockset
real, allocatable :: twomass_ID1(:),twomass_ID2(:),twomass_ra(:),twomass_dec(:)
real, allocatable :: twomass_gl(:),twomass_gb(:)
real, allocatable :: Ktmag(:),Jtmag(:),Kmagerr(:),Jmagerr(:)
real, allocatable :: E_BV(:),twomass_cz(:),twomass_ecz(:) 

real :: sdss_cover,twomass_cover,mock_cover,ratio,help_ratio,qso_ratio,av_red_2mass,av_J_err,av_Ks_err

real, allocatable :: gal_id(:,:),gx(:,:),gy(:,:),gz(:,:),gvx(:,:),gvy(:,:),gvz(:,:)
real, allocatable :: gnp(:,:),g_r(:,:),g_g(:,:),g_i(:,:),g_J(:,:),g_Ks(:,:)
real, allocatable :: g_fakeext_g(:,:),g_fakeext_r(:,:),g_redshift2mass(:,:)

real, allocatable :: mock_ra(:),mock_dec(:),mock_z(:),mock_mg(:),mock_mr(:),mock_ext_g(:),mock_ext_r(:)
logical, allocatable :: mock_visible(:)

real, allocatable :: mock2_ra(:),mock2_dec(:),mock2_z(:),mock2_mJ(:),mock2_mKs(:)
logical, allocatable :: mock2_visible(:)

real, allocatable :: Kmock_g(:),Kmock_r(:),Kmock_J(:),Kmock_Ks(:)
real, allocatable :: mockabsmag_g(:),mockabsmag_r(:),mockabsmag_J(:), mockabsmag_Ks(:)
real, allocatable :: mockdist(:),mock2dist(:)
real, allocatable :: mockxpos(:),mockypos(:),mockzpos(:),mock2xpos(:),mock2ypos(:),mock2zpos(:)
integer(kind=8), allocatable :: mockid(:),mockid2(:)

integer(kind=8), allocatable :: mgalid(:,:),mgalid2(:,:)

integer(kind=8), allocatable :: mmockid(:),mmockid2(:)

integer(kind=8), allocatable :: galid(:,:),galid2(:,:)

integer(kind=8), allocatable ::g_fofId(:,:),h_fofId(:,:)
integer(kind=8), allocatable :: mockfofid(:),mockfofid2(:)
integer(kind=8), allocatable :: hfofid(:),hfofid2(:)
integer, allocatable :: h_nsg(:,:),hmock_nsg(:)

real, allocatable :: gal_dist_c(:,:),gal_dist_l(:,:),gal_redshift(:,:),galpos_ra(:,:),galpos_dec(:,:)
real, allocatable :: gal_appmag_g(:,:),gal_appmag_r(:,:),gal_cosmored(:,:)
real, allocatable :: gal_appmag_J(:,:),gal_appmag_Ks(:,:)
logical, allocatable :: gal_visible(:,:),gal2_visible(:,:)
real ::  nx,ny,nz,err_redshift,err_pos,K_corr
real :: av_galext_g,sigma_galext_g,err_photo_g,av_galext_r,sigma_galext_r,err_photo_r
real ::  rand_angle,rand_value,rand_dec,rand_ra
real, dimension(0:7) :: snapshot_redshift
integer, dimension(1:6) :: snapshot_num,helpcount

real, allocatable :: hx(:,:),hy(:,:),hz(:,:),hnpart(:,:),hm(:,:),hdist(:,:),hredshift(:,:)
real, allocatable ::  gdist(:,:),gredshift(:,:)
real, allocatable :: g_mag_g(:),g_mag_r(:),g_mag_J(:),g_mag_Ks(:)
logical, allocatable :: halo_part_mock(:,:),gal_part_mock(:,:)
real, allocatable :: halo_x(:),halo_y(:),halo_z(:),halo_mass(:),halo_ra(:),halo_dec(:),halo_red(:)


real :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,moveinside,cubesize

real, dimension(1:6) :: mtot,m_ratio,partsum

logical, allocatable :: twomass_active(:),sdss_active(:),all_active(:)
real :: v_cmb,l_cmb,b_cmb,z_cmb,gal_b,gal_l,z_cmb_x,z_cmb_y,z_cmb_z,z_gal_x,z_gal_y,z_gal_z,arcsec
real :: limitmag_g,limitmag_r,limitmag_i,limitmag_z,limitmag_J,limitmag_H,limitmag_Ks
real :: limitmag_SDSS_official,limitmag_2MRS_official,completeness_2mrs
real, allocatable ::  helplist_J(:),helplist_H(:),helplist_Ks(:)

!K correction coefficients
real, dimension(0:5,0:3) :: K_coeff_r,K_coeff_i,K_coeff_z,K_coeff_J,K_coeff_Ks
real, dimension(0:7,0:3) :: K_coeff_g,K_coeff_H

real, dimension(0:5,0:3) :: invK_coeff_r,invK_coeff_g,invK_coeff_J,invK_coeff_Ks

real :: dummy_coeffient,d_J,e_J,f_J,d_Ks,e_Ks,f_Ks
real, allocatable :: K_g(:),K_r(:)

 real, allocatable :: magapp_g(:),magapp_r(:)
 real, allocatable :: K_J(:),K_Ks(:),z_2mass(:)
 real, allocatable :: magapp_J(:),magapp_Ks(:)
 real, allocatable :: distance_sdss(:),distance_2mass(:)
  real, allocatable :: loglumdist_sdss(:),loglumdist_2mass(:)
 real, allocatable :: magabs_g(:),magabs_r(:)
 real, allocatable :: magabs_J(:),magabs_Ks(:)
 real, allocatable :: xpos_sdss(:),ypos_sdss(:),zpos_sdss(:)
  real, allocatable :: xpos_2mass(:),ypos_2mass(:),zpos_2mass(:)
  


 real, allocatable ::g_gDust(:,:),g_rDust(:,:),g_iDust(:,:)
 real, allocatable :: h_m_Crit200(:,:),h_r_crit200(:,:)

    
   real ::  fiber_coll_r,decollided_sampling,loss_coll,final_sampling,angular_sep,masspart,usedsidelength
  
integer, dimension(0:500,0:400) :: mapmap
real :: help_x2,help_y2
integer :: help_x,help_y

real, dimension(0:110) :: bin_sdss,bin_2mass
real, dimension(0:2050,0:2050) :: red_map

real:: saturation_u,saturation_g,saturation_r,saturation_i,saturation_z,red_coeff_g,red_coeff_r
real:: rot_anglez,rot_angley,rot_x,rot_y,rot_z,rot_dec,rot_ra,rx,ry,rz,proj_x,proj_y

real ::gasdev!1,gasdev2,gasdev3,gasdev4,gasdev5,gasdev6,gasdev7,gasdev8

integer, dimension(1:8) :: u_array
integer, dimension(1:1) :: u_local
integer :: u


! !parallization variables		
! integer ierr, rankrank, nbnodes, namelen
!  character (len=MPI_MAX_PROCESSOR_NAME) :: name
! !parallization procedures
! call MPI_INIT(ierr)
! call MPI_COMM_RANK(MPI_COMM_WORLD, rankrank, ierr)
! call MPI_COMM_SIZE(MPI_COMM_WORLD, nbnodes, ierr)
! call MPI_GET_PROCESSOR_NAME(name, namelen, ierr)	

uuu=1



OPEN(33,file='logfile_ts_sor.txt')


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme MOCK CATALOGUE CREATOR started'
WRITE(*,*) '============================================================'
WRITE(33,*) '============================================================'
WRITE(33,*) '    programme MOCK CATALOGUE CREATOR started'
WRITE(33,*) '============================================================'

! random seed
 CALL SYSTEM_CLOCK(hcount)
 hcount=hcount-INT(hcount/100000)*100000
 CALL srand(hcount) 


! define constants
PI=ACOS(-1.D0)
arcsec=2.D0/3600.D0

n_mockcat=8

moveinside=10.D0
 cubesize=400.D0
usedsidelength=cubesize-2.D0*moveinside

n_qso=1889
sdss_cover=9376.D0/41253.D0
twomass_cover=0.91D0
mock_cover=1.D0/8.D0

fiber_coll_r=55.D0/3600.D0
decollided_sampling=0.99D0
loss_coll=0.D0
final_sampling=0.92D0
 completeness_2mrs=0.976D0

limitmag_SDSS_official=17.77D0
limitmag_2MRS_official=11.75D0


Omega_m=0.25D0
Omega_l=0.75D0
 
 masspart=8.61D8
 

q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
H0=73.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)
 cV_mill=(cubesize*1.D6/h)**3
m_expected=cV_mill*rho_crit*Omega_m
part_exp=(2160.D0**3)*((cubesize/500.D0)**3)


DO i=0,2050
DO ii=0,2050
red_map(i,ii)=0.D0
END DO
END DO


OPEN(52,file='snapshot_redshift_list.txt')
DO i=1,6
READ(52,*) iii,snapshot_redshift(i)
snapshot_num(i)=iii
! WRITE(*,*) snapshot_num(i),snapshot_redshift(i)
END DO
CLOSE(52)

snapshot_redshift(0)=0.D0
snapshot_redshift(7)=1.D0


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
OPEN(50,file='SDSS_DR12.txt')
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

 WRITE(*,*) 'arrays allocated'
 WRITE(33,*) 'arrays allocated'


 
! read file
OPEN(50,file='SDSS_DR12.txt')
READ(50,*)
DO i=1,n_sdss

READ(50,*) dummy_id,sdss_ra(i),sdss_dec(i),sdss_gb(i),sdss_gl(i),sdss_z(i),&
 cModelMag_g(i),cModelMag_r(i),dummy_id,cModelMagErr_g(i),cModelMagErr_r(i),dummy_id,&
 extinction_g(i),extinction_r(i),dummy_id

sdss_active(i)=.TRUE.

END DO

 CLOSE(50)


! read file
OPEN(50,file='SDSS_DR12.txt')
READ(50,*)
DO i=1,n_sdss 
READ(50,*) read_string
sdss_objID(i)=read_string(1:19)
END DO

 CLOSE(50)


WRITE(*,*) n_sdss,'galaxies from SDSS DR12'
WRITE(33,*) n_sdss,'galaxies from SDSS DR12'


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






  WRITE(*,*) '---------------------------------------------------------------'
 WRITE(33,*) '---------------------------------------------------------------'

OPEN(50,file='map_north.txt')
DO i=1,2048 
READ(50,*) red_map(i,1:2048)
END DO
 CLOSE(50)



OPEN(50,file='map_north.txt')
READ(50,*) dummy_coeffient
READ(50,*) red_coeff_g
READ(50,*) red_coeff_r
READ(50,*) dummy_coeffient
READ(50,*) dummy_coeffient
 CLOSE(50)







!prepare for calculation of correction for CMB motion
! ! values for the motion of the sun relativ to the CMB
! v_cmb=369.D0
! l_cmb=263.99D0
! b_cmb=48.26D0
! milkyway system
  v_cmb=-220.D0
  l_cmb=270.D0
  b_cmb=0.D0

! calculate sun's motion in z space relativ to the CMB
z_cmb=v_cmb/light
gal_b=b_cmb*PI/180.D0
gal_l=l_cmb*PI/180.D0
z_cmb_x=z_cmb*COS(gal_b)*COS(gal_l)
z_cmb_y=z_cmb*COS(gal_b)*SIN(gal_l)
z_cmb_z=z_cmb*SIN(gal_b)


DO i=1,n_sdss
z_cor(i)=sdss_z(i)

gal_b=sdss_gb(i)*PI/180.D0
gal_l=sdss_gl(i)*PI/180.D0
z_gal_x=sdss_z(i)*COS(gal_b)*COS(gal_l)
z_gal_y=sdss_z(i)*COS(gal_b)*SIN(gal_l)
z_gal_z=sdss_z(i)*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
z_cor(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


DO i=1,n_2mass
z_2mass(i)=twomass_cz(i)/light

gal_b=twomass_gb(i)*PI/180.D0
gal_l=twomass_gl(i)*PI/180.D0
z_gal_x=twomass_cz(i)/light*COS(gal_b)*COS(gal_l)
z_gal_y=twomass_cz(i)/light*COS(gal_b)*SIN(gal_l)
z_gal_z=twomass_cz(i)/light*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
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
 
 
 OPEN(52,file='saturation.txt')
READ(52,*) saturation_u
READ(52,*) saturation_g
READ(52,*) saturation_r
READ(52,*) saturation_i
READ(52,*) saturation_z
 CLOSE(52)
 

 !get parameters for inverse K correction  
OPEN(52,file='invK_correction_g.txt')
DO i=0,5
READ(52,*) invK_coeff_g(i,0:3)
END DO
 CLOSE(52)
OPEN(52,file='invK_correction_r.txt')
DO i=0,5
READ(52,*) invK_coeff_r(i,0:3)
END DO
 CLOSE(52)


OPEN(52,file='invK_correction_J.txt')
DO i=0,5
READ(52,*) invK_coeff_J(i,0:3)
END DO
 CLOSE(52)

OPEN(52,file='invK_correction_Ks.txt')
DO i=0,5
READ(52,*) invK_coeff_Ks(i,0:3)
END DO
 CLOSE(52)
 
 
!  WRITE(*,*) invK_coeff_g(0:5,0:3)
!  WRITE(*,*) invK_coeff_r(0:5,0:3)
!   WRITE(*,*) invK_coeff_J(0:5,0:3)
!    WRITE(*,*) invK_coeff_Ks(0:5,0:3)

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







 hcount=0
 DO i=1,n_2mass
 IF (twomass_active(i)) THEN
 hcount=hcount+1
 helplist_J(hcount)=magapp_J(i)
 helplist_Ks(hcount)=magapp_Ks(i)
 END IF
 END DO


 DO ii=1,hcount
 DO iii=ii,hcount

 IF (helplist_J(ii)<helplist_J(iii)) THEN
 helpflip=helplist_J(ii)
 helplist_J(ii)=helplist_J(iii)
 helplist_J(iii)=helpflip
 END IF

 IF (helplist_Ks(ii)<helplist_Ks(iii)) THEN
 helpflip=helplist_Ks(ii)
 helplist_Ks(ii)=helplist_Ks(iii)
 helplist_Ks(iii)=helpflip
 END IF


 END DO
 END DO

nedge=CEILING(DBLE(hcount)/1000.D0*3D0) !99.7% limit =~3 sigma

limitmag_J=helplist_J(nedge)


limitmag_Ks=helplist_Ks(nedge)

! IF (uuu==1) THEN
! OPEN(50,file='limitmag_2mass.txt')
! 
! WRITE(50,*) limitmag_J,limitmag_Ks
! 
! CLOSE(50)
! 
! END IF

DO i=1,n_2mass



IF (magapp_Ks(i)>(limitmag_2MRS_official+0.5D0)) THEN
twomass_active(i)=.FALSE.
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

DO i=1,n_sdss
IF (sdss_active(i)) THEN
gal_b=sdss_dec(i)*PI/180.D0
gal_l=sdss_ra(i)*PI/180.D0
xpos_sdss(i)=distance_sdss(i)*COS(gal_b)*COS(gal_l)
ypos_sdss(i)=distance_sdss(i)*COS(gal_b)*SIN(gal_l)
zpos_sdss(i)=distance_sdss(i)*SIN(gal_b)
END IF
END DO

DO i=1,n_2mass
IF (twomass_active(i)) THEN
gal_b=twomass_dec(i)*PI/180.D0
gal_l=twomass_ra(i)*PI/180.D0
xpos_2mass(i)=distance_2mass(i)*COS(gal_b)*COS(gal_l)
ypos_2mass(i)=distance_2mass(i)*COS(gal_b)*SIN(gal_l)
zpos_2mass(i)=distance_2mass(i)*SIN(gal_b)
END IF
END DO


! read file
! OPEN(50,file='SDSS_real/SDSS_obs_galaxies_pos3D.txt')
! DO i=1,n_sdss
! IF (sdss_active(i)) THEN
! WRITE(50,*) xpos_sdss(i),ypos_sdss(i),zpos_sdss(i),magabs_g(i),magabs_r(i)
! END IF
! END DO
! ! read file
! OPEN(50,file='2MRS_real/2MRS_obs_galaxies_pos3D.txt')
! DO i=1,n_2mass
! IF (twomass_active(i)) THEN
! WRITE(50,*) xpos_2mass(i),ypos_2mass(i),zpos_2mass(i),magabs_J(i),magabs_Ks(i)
! END IF
! END DO
! 

! read file
OPEN(50,file='ts_sor/SDSS_real/SDSS_obs_galaxies_radecz.txt')
DO i=1,n_sdss
IF (sdss_active(i)) THEN
WRITE(50,*) sdss_ra(i),sdss_dec(i),z_cor(i),magabs_g(i),magabs_r(i)
END IF
END DO


! 
! OPEN(50,file='SDSS_real/Prajwal.txt')
! DO i=1,n_sdss
! IF (sdss_active(i)) THEN
! WRITE(50,*) sdss_objID(i),sdss_ra(i),sdss_dec(i),z_cor(i),
! END IF
! END DO
! 


! read file
OPEN(50,file='ts_sor/2MRS_real/2MRS_obs_galaxies_radecz.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) twomass_ra(i),twomass_dec(i),z_2mass(i),magabs_J(i),magabs_Ks(i)
END IF
END DO




! read file
OPEN(50,file='ts_sor/SDSS_real/SDSS_id_file.txt')
DO i=1,n_sdss
IF (sdss_active(i)) THEN
WRITE(50,*) sdss_objID(i)
END IF
END DO

! read file
OPEN(50,file='ts_sor/2MRS_real/2MRS_id_file.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) twomass_ID1(i),twomass_ID2(i)
END IF
END DO



WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'
WRITE(33,*) '============================================================'
WRITE(33,*) '    programme complete'
WRITE(33,*) '============================================================'

CLOSE(33)






END PROGRAM







! function for Gaussian error
REAL FUNCTION gasdev(idum)

INTEGER :: idum
INTEGER :: iset
REAL :: fac,gset,rsq,v1,v2,ran1

SAVE iset,gset
DATA iset/0/



	IF (idum.lt.0) 	THEN
		iset=0
	END IF


	IF (iset.eq.0) THEN
		rsq=1.
		DO While (rsq.ge.1..or.rsq.eq.0)
			v1=2.*rand(idum)-1.
			v2=2.*rand(idum)-1.
			rsq=v1**2+v2**2
		END DO
		fac=SQRT(-2.*LOG(rsq)/rsq)
		gset=v1*fac
		gasdev=v2*fac
		iset=1
	ELSE
		gasdev=gset
		iset=0
	END IF

RETURN
END FUNCTION gasdev
