PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
! include 'mpif.h'



double precision :: PI,H0,q0,light,G,tophat,rho_crit,area_cover,dummy_var,P_help
integer(kind=8) :: io_err,n,i,ii,iii
character(200) :: filename,nm,apx
integer :: uuu,n_max,pos_index

integer, allocatable :: n_mock(:)

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor,R_factor,basic_link
double precision :: d_red,d_ang_dist,angular_sep,delta_z,hreal,l_next,hl_z,hl_a,l_old,mabs_vollim
double precision ::mass_halo_fit,twomass_cover,mock_cover,sdss_cover,dmax,dmin

double precision, allocatable ::  ra(:,:),dec(:,:),z(:,:),mag_g(:,:),mag_r(:,:)


double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y
integer, dimension(1:102) :: binclustermock,binclusterfof
double precision, dimension(0:110,0:8) :: bin_z_n,bin_vol,binz




! define constants
PI=ACOS(-1.D0)



apx='SDSS'
 
 
 Omega_m=0.25D0
 Omega_l=0.75D0
 H0=73.D0
 
 
q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
sdss_cover=9274.D0/41253.D0
twomass_cover=0.91D0
mock_cover=1.D0/8.D0



apx=TRIM(adjustl(apx))


allocate(n_mock(0:8))
n_max=0

DO uuu=0,8

WRITE(nm,*) uuu
nm=adjustl(nm)
 
 

 
!   WRITE(*,*) 'now in CPU',uuu
 
 IF (uuu==0) THEN
  filename=TRIM(apx)//'_real/'//TRIM(apx)//'_obs_galaxies_radecz.txt'
 ELSE
filename=TRIM(apx)//'_mock'//TRIM(nm)//'/'//TRIM(apx)//'_galaxies_mock'//TRIM(nm)//'_final.txt'
END IF
! get length of file
OPEN(50,file=TRIM(filename))
io_err=0
n_mock(uuu)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock(uuu)=n_mock(uuu)+1
END DO
 CLOSE(50)
n_mock(uuu)=n_mock(uuu)-1

! WRITE(*,*) n,'galaxies will be used in',uuu
IF (n_mock(uuu)>n_max) THEN
n_max=n_mock(uuu)
END IF

END DO




allocate(ra(1:n_max,0:8))
allocate(dec(1:n_max,0:8))
allocate(z(1:n_max,0:8))
allocate(mag_g(1:n_max,0:8))
allocate(mag_r(1:n_max,0:8))



DO uuu=0,8

WRITE(nm,*) uuu
nm=adjustl(nm)
 

 
 
 IF (uuu==0) THEN
  filename=TRIM(apx)//'_real/'//TRIM(apx)//'_obs_galaxies_radecz.txt'
 ELSE
filename=TRIM(apx)//'_mock'//TRIM(nm)//'/'//TRIM(apx)//'_galaxies_mock'//TRIM(nm)//'_final.txt'
END IF


! read file
OPEN(50,file=TRIM(filename))
DO i=1,n_mock(uuu)
READ(50,*) ra(i,uuu),dec(i,uuu),z(i,uuu),mag_g(i,uuu),mag_r(i,uuu)
END DO
CLOSE(50)

DO i=0,110
bin_z_n(i,uuu)=0.D0
binz(i,uuu)=(DBLE(i)+0.5D0)/1000.D0


divisor=(SQRT(1.D0+2.D0*q0*(binz(i,uuu)-0.0005))+1.D0+q0*(binz(i,uuu)-0.0005))
dmin=light/H0*(binz(i,uuu)-0.0005)*(1.D0+(((binz(i,uuu)-0.0005)*(1.D0-q0))/divisor))
dmin=dmin*((1.D0+(binz(i,uuu)-0.0005))**(-1.D0))

divisor=(SQRT(1.D0+2.D0*q0*(binz(i,uuu)+0.0005))+1.D0+q0*(binz(i,uuu)+0.0005))
dmax=light/H0*(binz(i,uuu)+0.0005)*(1.D0+(((binz(i,uuu)+0.0005)*(1.D0-q0))/divisor))
dmax=dmax*((1.D0+(binz(i,uuu)+0.0005))**(-1.D0))


IF (uuu==0) THEN
IF (TRIM(apx)=='SDSS') THEN
area_cover=sdss_cover
ELSE
area_cover=twomass_cover
END IF
ELSE
area_cover=mock_cover
END IF

bin_vol(i,uuu)=4.D0/3.D0*PI*area_cover*((dmax**3)-(dmin**3))

END DO

DO i=1,n_mock(uuu)
pos_index=NINT((z(i,uuu)*1000.D0)-0.5D0) 
IF ((pos_index>-1).AND.(pos_index<111)) THEN
bin_z_n(pos_index,uuu)=bin_z_n(pos_index,uuu)+1.D0
END IF

END DO

DO i=0,110
bin_z_n(i,uuu)=bin_z_n(i,uuu)/bin_vol(i,uuu)
END DO



END DO



! read file
OPEN(50,file='bin_ndens_'//TRIM(apx)//'.txt')
DO i=1,109
WRITE(50,*) binz(i,0),bin_z_n(i,0:8)
END DO
CLOSE(50)





END PROGRAM


 