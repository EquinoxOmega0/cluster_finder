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



OPEN(33,file='logfile.txt')


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





DO ii=1,6
WRITE(helpstring,*) snapshot_num(ii)
helpstring=adjustl(helpstring)
! get length of file
OPEN(50,file='all_fof/fof'//TRIM(helpstring)//'.txt')
io_err=0
n_halo(ii)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_halo(ii)=n_halo(ii)+1
END DO
 CLOSE(50)
n_halo(ii)=n_halo(ii)-1
END DO




DO ii=1,6
WRITE(helpstring,*) snapshot_num(ii)
helpstring=adjustl(helpstring)
! get length of file
OPEN(50,file='all_galaxies/gal'//TRIM(helpstring)//'.txt')
io_err=0
n_gal(ii)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_gal(ii)=n_gal(ii)+1
END DO
 CLOSE(50)
n_gal(ii)=n_gal(ii)-1
END DO



n_gal_max=0
DO ii=1,6
IF (n_gal(ii)>n_gal_max) THEN
n_gal_max=n_gal(ii)
END IF
END DO

n_halo_max=0
DO ii=1,6
IF (n_halo(ii)>n_halo_max) THEN
n_halo_max=n_halo(ii)
END IF
END DO
 WRITE(*,*) 'all file lengths measured'
 WRITE(33,*) 'all file lengths measured'


allocate(gx(1:n_gal_max,1:6))
allocate(gy(1:n_gal_max,1:6))
allocate(gz(1:n_gal_max,1:6))
allocate(gvx(1:n_gal_max,1:6))
allocate(gvy(1:n_gal_max,1:6))
allocate(gvz(1:n_gal_max,1:6))
allocate(gnp(1:n_gal_max,1:6))
allocate(g_r(1:n_gal_max,1:6))
allocate(g_g(1:n_gal_max,1:6))
allocate(g_i(1:n_gal_max,1:6))
allocate(g_J(1:n_gal_max,1:6))
allocate(g_Ks(1:n_gal_max,1:6))

allocate(g_fofID(1:n_gal_max,1:6))
allocate(g_gDust(1:n_gal_max,1:6))
allocate(g_rDust(1:n_gal_max,1:6))
allocate(g_iDust(1:n_gal_max,1:6))
 
allocate(g_fakeext_g(1:n_gal_max,1:6))
allocate(g_fakeext_r(1:n_gal_max,1:6))
allocate(g_redshift2mass(1:n_gal_max,1:6))
 
allocate(gal_dist_c(1:n_gal_max,1:6))
allocate(gal_dist_l(1:n_gal_max,1:6))
allocate(gal_redshift(1:n_gal_max,1:6))
allocate(gal_appmag_r(1:n_gal_max,1:6))
allocate(gal_appmag_g(1:n_gal_max,1:6))
allocate(galpos_ra(1:n_gal_max,1:6))
allocate(galpos_dec(1:n_gal_max,1:6))
allocate(gal_visible(1:n_gal_max,1:6))
allocate(gal2_visible(1:n_gal_max,1:6))
allocate(gal_cosmored(1:n_gal_max,1:6))
allocate(gal_part_mock(1:n_gal_max,1:6))

allocate(gal_appmag_J(1:n_gal_max,1:6))
allocate(gal_appmag_Ks(1:n_gal_max,1:6))

allocate(gdist(1:n_gal_max,1:6))
allocate(gredshift(1:n_gal_max,1:6))

allocate(galid(1:n_gal_max,1:6))
allocate(galid2(1:n_gal_max,1:6))

allocate(hx(1:n_halo_max,1:6))
allocate(hy(1:n_halo_max,1:6))
allocate(hz(1:n_halo_max,1:6))
allocate(hnpart(1:n_halo_max,1:6))
allocate(hm(1:n_halo_max,1:6))
allocate(hredshift(1:n_halo_max,1:6))
allocate(hdist(1:n_halo_max,1:6))
allocate(halo_part_mock(1:n_halo_max,1:6))

allocate(h_fofId(1:n_halo_max,1:6))
allocate(h_m_Crit200(1:n_halo_max,1:6))
allocate(h_r_crit200(1:n_halo_max,1:6))
allocate(h_nsg(1:n_halo_max,1:6))

allocate(mgalid(1:n_gal_max,1:6))
allocate(mgalid2(1:n_gal_max,1:6))


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


DO ii=1,6
DO i=1,n_gal_max
galid(i,ii)=i+(ii-1)*n_gal_max
galid2(i,ii)=i+(ii-1)*n_gal_max
END DO 
END DO


DO ii=1,6
WRITE(helpstring,*) snapshot_num(ii)
helpstring=adjustl(helpstring)
! get length of file
OPEN(50,file='all_galaxies/gal'//TRIM(helpstring)//'.txt')
! READ(50,*) 
DO i=1,n_gal(ii)
READ(50,*) mgalid(i,ii),g_fofID(i,ii),gx(i,ii),gy(i,ii),gz(i,ii),gvx(i,ii),gvy(i,ii),gvz(i,ii),&
gnp(i,ii),g_gDust(i,ii),g_rDust(i,ii),g_iDust(i,ii)
END DO

CLOSE(50)

END DO

DO ii=1,6
DO i=1,n_gal(ii)
g_g(i,ii)=g_gDust(i,ii)
g_r(i,ii)=g_rDust(i,ii)
g_i(i,ii)=g_iDust(i,ii)
END DO
END DO


hcount=0
DO ii=1,6
hcount=hcount+n_gal(ii)
END DO
 WRITE(*,*) hcount,'galaxies from the Millennium simulation'
 WRITE(33,*) hcount,'galaxies from the Millennium simulation'
 
 
DO ii=1,6
DO i=1,n_gal_max
mgalid2(i,ii)=mgalid(i,ii)
END DO 
END DO




DO ii=1,6
WRITE(helpstring,*) snapshot_num(ii)
helpstring=adjustl(helpstring)
! get length of file
OPEN(50,file='all_fof/fof'//TRIM(helpstring)//'.txt')
! READ(50,*) 
DO i=1,n_halo(ii)
READ(50,*) h_fofId(i,ii),hx(i,ii),hy(i,ii),hz(i,ii),h_m_Crit200(i,ii),h_r_crit200(i,ii),h_nsg(i,ii),hnpart(i,ii)
END DO
CLOSE(50)

END DO


hcount=0
DO ii=1,6
hcount=hcount+n_halo(ii)
END DO
 WRITE(*,*) hcount,'halos from the Millennium simulation'
 WRITE(33,*) hcount,'halos from the Millennium simulation'


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






  WRITE(*,*) '---------------------------------------------------------------'
 WRITE(33,*) '---------------------------------------------------------------'
 WRITE(*,*) 'mass completeness'
 WRITE(33,*) 'mass completeness'
 
 DO ii=1,6
  partsum(ii)=0.D0
 DO i=1,n_halo(ii)
  partsum(ii)=partsum(ii)+hnpart(i,ii)
 END DO
 partsum(ii)=partsum(ii)/part_exp*100.D0
 END DO
  DO ii=1,6
  
 WRITE(*,*) snapshot_num(ii),':',partsum(ii),'%'
 WRITE(33,*) snapshot_num(ii),':',partsum(ii),'%'
 END DO
 
 
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
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
z_cor(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


DO i=1,n_2mass
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

IF (uuu==1) THEN
OPEN(50,file='limitmag_2mass.txt')

WRITE(50,*) limitmag_J,limitmag_Ks

CLOSE(50)

END IF

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

IF (uuu==1) THEN

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
OPEN(50,file='SDSS_real/SDSS_obs_galaxies_radecz.txt')
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
OPEN(50,file='2MRS_real/2MRS_obs_galaxies_radecz.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) twomass_ra(i),twomass_dec(i),z_2mass(i),magabs_J(i),magabs_Ks(i)
END IF
END DO




! read file
OPEN(50,file='SDSS_real/SDSS_id_file.txt')
DO i=1,n_sdss
IF (sdss_active(i)) THEN
WRITE(50,*) sdss_objID(i)
END IF
END DO

! read file
OPEN(50,file='2MRS_real/2MRS_id_file.txt')
DO i=1,n_2mass
IF (twomass_active(i)) THEN
WRITE(50,*) twomass_ID1(i),twomass_ID2(i)
END IF
END DO






DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,n_sdss
IF (sdss_active(i)) THEN

help_x=NINT((loglumdist_sdss(i)-6.D0)/3.D0*500.D0)
help_y=NINT((-magabs_r(i)-10.D0)/20.D0*400.D0)
IF ((help_x>0).AND.(help_x<500)) THEN
IF ((help_y>0).AND.(help_y<400)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF

END IF
END DO


OPEN(61,file='sdss_absmag_vs_d.txt')
DO i=0,500
DO ii=0,400
help_x2=((DBLE(i))*3.D0/500.D0)+6.D0
help_y2=-((DBLE(ii))*20.D0/400.D0)-10.D0
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

help_x=NINT((loglumdist_2mass(i)-6.D0)/3.D0*500.D0)
help_y=NINT((-magabs_Ks(i)-10.D0)/20.D0*400.D0)
IF ((help_x>0).AND.(help_x<500)) THEN
IF ((help_y>0).AND.(help_y<400)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF

END IF
END DO


OPEN(61,file='2mass_absmag_vs_d.txt')
DO i=0,500
DO ii=0,400
help_x2=((DBLE(i))*3.D0/500.D0)+6.D0
help_y2=-((DBLE(ii))*20.D0/400.D0)-10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)



 


 
DO i=0,110
bin_sdss(i)=0
bin_2mass(i)=0
END DO

DO i=1,n_sdss

IF (sdss_active(i)) THEN
DO ii=0,110
binhelp=DBLE(ii)/1000.D0+0.0005
IF ((z_cor(i)>(binhelp-0.0005)).AND.(z_cor(i)<(binhelp+0.0005))) THEN
bin_sdss(ii)=bin_sdss(ii)+1  
END IF
END DO
END IF

END DO





DO i=1,n_2mass

IF (twomass_active(i)) THEN
DO ii=0,110
binhelp=DBLE(ii)/1000.D0+0.0005
IF ((z_2mass(i)>(binhelp-0.0005)).AND.(z_2mass(i)<(binhelp+0.0005))) THEN
bin_2mass(ii)=bin_2mass(ii)+1   
END IF
END DO
END IF

END DO


DO ii=0,110
zmax=DBLE(ii+1)/1000.D0
zmin=DBLE(ii)/1000.D0

divisor=(SQRT(1.D0+2.D0*q0*zmax)+1.D0+q0*zmax)

disthelp2=light/H0*zmax*(1.D0+((zmax*(1.D0-q0))/divisor))


disthelp2=disthelp2*((1.D0+zmax)**(-1.D0))

divisor=(SQRT(1.D0+2.D0*q0*zmin)+1.D0+q0*zmin)

disthelp1=light/H0*zmin*(1.D0+((zmin*(1.D0-q0))/divisor))

disthelp1=disthelp1*((1.D0+zmin)**(-1.D0))

 cmv=4.D0*PI/3.D0*((disthelp2**3)-(disthelp1**3))
 

bin_sdss(ii)=bin_sdss(ii)/(cmv*(sdss_cover))
bin_2mass(ii)=bin_2mass(ii)/(cmv*twomass_cover)


END DO



OPEN(70,file='redshift_distribution.txt')
DO i=0,110
binhelp=DBLE(i)/1000.D0+0.0005
WRITE(70,*) binhelp,bin_sdss(i),bin_2mass(i)
END DO
 CLOSE (70)


 
 
 
 
 
 END IF
 
 
 
 
 
 
 
 
 
 DO ii=1,6
 DO i=1,n_gal(ii)
 gx(i,ii)=(gx(i,ii)-moveinside)/h
 gy(i,ii)=(gy(i,ii)-moveinside)/h
 gz(i,ii)=(gz(i,ii)-moveinside)/h
 END DO
 END DO
 
    DO ii=1,6
  
 DO i=1,n_halo(ii)
 hx(i,ii)=(hx(i,ii)-moveinside)/h
 hy(i,ii)=(hy(i,ii)-moveinside)/h
 hz(i,ii)=(hz(i,ii)-moveinside)/h
 END DO
 
 END DO
 
!  WRITE(*,*) 'corners ok',uuu
 
err_photo_g=0.D0
err_photo_r=0.D0
av_galext_g=0.D0
av_galext_r=0.D0

av_red_2mass=0.D0
av_J_err=0.D0
av_Ks_err=0.D0

hcount=0
 DO i=1,n_sdss
 IF (sdss_active(i)) THEN
av_galext_g=av_galext_g+extinction_g(i)
av_galext_r=av_galext_r+extinction_r(i)
END IF
 END DO
  av_galext_g=av_galext_g/DBLE(n_sdss_used)
   av_galext_r=av_galext_r/DBLE(n_sdss_used)
 
DO i=1,n_2mass
IF (twomass_active(i)) THEN
av_red_2mass=av_red_2mass+twomass_ecz(i)/light
av_J_err=av_J_err+Jmagerr(i)
av_Ks_err=av_Ks_err+Kmagerr(i)
END IF
END DO
 
av_red_2mass=av_red_2mass/DBLE(n_2mass_used)
av_J_err=av_J_err/DBLE(n_2mass_used)
av_Ks_err=av_Ks_err/DBLE(n_2mass_used)

 
 
 DO i=1,n_sdss
IF (sdss_active(i)) THEN
err_photo_g=err_photo_g+cModelMagErr_g(i)
END IF
 END DO

 err_photo_g=err_photo_g/DBLE(n_sdss_used) 

  DO i=1,n_sdss
IF (sdss_active(i)) THEN
err_photo_r=err_photo_r+cModelMagErr_r(i)
END IF
 END DO

 err_photo_r=err_photo_r/DBLE(n_sdss_used) 


 sigma_galext_g=0.D0
  DO i=1,n_sdss
  IF (sdss_active(i)) THEN
sigma_galext_g=sigma_galext_g+((extinction_g(i)-av_galext_g)**2)
END IF
   END DO

    sigma_galext_g=sigma_galext_g/DBLE(n_sdss_used) 
    sigma_galext_g=SQRT(sigma_galext_g)
   
   
   
   
   sigma_galext_r=0.D0
  DO i=1,n_sdss
    IF (sdss_active(i)) THEN
sigma_galext_r=sigma_galext_r+((extinction_r(i)-av_galext_r)**2)
END IF
   END DO

    sigma_galext_r=sigma_galext_r/DBLE(n_sdss_used) 
    sigma_galext_r=SQRT(sigma_galext_r)
   
   
   
err_redshift=30.D0
err_pos=0.1D0/3600.D0/180.D0*PI

!  WRITE(*,*) 'averages ok',uuu




OPEN(77,file='error_2MRS.txt')
WRITE(77,*) 'z_err',av_red_2mass 
WRITE(77,*) 'J_err',av_J_err 
WRITE(77,*) 'Ks_err',av_Ks_err 

CLOSE(77)


OPEN(77,file='error_SDSS.txt')
WRITE(77,*) 'z_err',err_redshift 
WRITE(77,*) 'g_err',err_photo_g 
WRITE(77,*) 'r_err',err_photo_r
WRITE(77,*) 'extinction_g',av_galext_g,sigma_galext_g
WRITE(77,*) 'extinction_r',av_galext_r,sigma_galext_r
CLOSE(77)





 WRITE(*,*) 'start parallized part'

! DO u=1,8
! u_array(u)=u
! END DO
! 
! call MPI_Scatter(u_array,1,MPI_INTEGER,u_local,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! 
! uuu=u_local(1)

DO uuu=1,8

WRITE(number_mockset,*) uuu
number_mockset=adjustl(number_mockset)
 
 
 OPEN(34,file='logfile'//TRIM(number_mockset)//'.txt')
 
!   WRITE(*,*) 'now in CPU',uuu
 
 WRITE(*,*) 'now in mock',uuu
 


! 000
IF (uuu==1) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
  gx(i,ii)=gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=gz(i,ii)
 END DO
 END DO 
END IF

!100
IF (uuu==2) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
  gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=gz(i,ii)
 END DO
 END DO 
END IF

!010
IF (uuu==3) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
 
 gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=gz(i,ii)
 
  gx(i,ii)=gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=gz(i,ii)
 END DO
 END DO 
END IF

!001
IF (uuu==4) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
 
  gx(i,ii)=gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=gz(i,ii)
 
 
  gx(i,ii)=gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 END DO
 END DO 
END IF

!110
IF (uuu==5) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
 
 
gx(i,ii)=gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 
  gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=gz(i,ii)
 END DO
 END DO 
END IF

!011
IF (uuu==6) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
 
  gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=gz(i,ii)
 
 
  gx(i,ii)=gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 END DO
 END DO 
END IF

!101
IF (uuu==7) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
 
  gx(i,ii)=gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 
 
  gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 END DO
 END DO 
END IF

!111
IF (uuu==8) THEN
 DO ii=1,6
 DO i=1,n_gal(ii)
 
  gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 
 
  gx(i,ii)=(usedsidelength/h)-gx(i,ii)
  gy(i,ii)=(usedsidelength/h)-gy(i,ii)
  gz(i,ii)=(usedsidelength/h)-gz(i,ii)
 END DO
 END DO 
END IF


!  WRITE(*,*) 'corners ok',uuu

IF ((uuu>0).AND.(uuu<9)) THEN
 
 DO ii=1,6
 
 DO i=1,n_gal(ii)
   gal_visible(i,ii)=.TRUE.
   END DO
 
 
 DO i=1,n_gal(ii)
 !comoving distance
gal_dist_c(i,ii)=((gx(i,ii))**2)+((gy(i,ii))**2)+((gz(i,ii))**2)
gal_dist_c(i,ii)=SQRT(gal_dist_c(i,ii))

! END DO
! END DO
! 
!  WRITE(*,*) 'distances calculated',uuu
!  DO ii=1,6
!  
!  DO i=1,n_gal(ii)
!redshift
gal_redshift(i,ii)=(light**2)*((light**2)+(gal_dist_c(i,ii)**2)*(H0**2)*(1.D0-2.D0*q0))*(q0-1.D0)**2
gal_redshift(i,ii)=SQRT(gal_redshift(i,ii))
gal_redshift(i,ii)=gal_redshift(i,ii)+light*gal_dist_c(i,ii)*H0+(light**2)*(q0-1.D0)-(gal_dist_c(i,ii)**2)*(H0**2)*(q0**2)
gal_redshift(i,ii)=gal_redshift(i,ii)/((light-gal_dist_c(i,ii)*H0*q0)**2)

gal_cosmored(i,ii)=gal_redshift(i,ii)

 IF (0.13D0<gal_cosmored(i,ii)) THEN
 gal_visible(i,ii)=.FALSE.
 END IF

 IF (gal_visible(i,ii)) THEN
!luminosity distance
gal_dist_l(i,ii)=gal_dist_c(i,ii)*(1.D0+gal_redshift(i,ii))


nx=(gx(i,ii))/gal_dist_c(i,ii)
ny=(gy(i,ii))/gal_dist_c(i,ii)
nz=(gz(i,ii))/gal_dist_c(i,ii)
!peculiar motion

z_pec=(nx*gvx(i,ii)+ny*gvy(i,ii)+nz*gvz(i,ii))/light


gal_redshift(i,ii)=((1.D0+gal_redshift(i,ii))*(1.D0+z_pec))-1.D0

g_redshift2mass(i,ii)=gal_redshift(i,ii)

! END DO
! END DO
! WRITE(*,*) 'redshifts calculated',uuu
!  
!  
!  DO ii=1,6
!  
!  DO i=1,n_gal(ii)



g_J(i,ii)=g_g(i,ii)-(d_J*(g_g(i,ii)-g_r(i,ii))+e_J*(g_r(i,ii)-g_i(i,ii))+f_J)
g_Ks(i,ii)=g_g(i,ii)-(d_Ks*(g_g(i,ii)-g_r(i,ii))+e_Ks*(g_r(i,ii)-g_i(i,ii))+f_Ks)

! END DO
! END DO
! 
! WRITE(*,*) '2Mass transformed',uuu
!   
!  DO ii=1,6
!  
!  DO i=1,n_gal(ii)

IF (gal_dist_l(i,ii)>0.D0) THEN
gal_appmag_r(i,ii)=g_r(i,ii)+5.D0*LOG10(gal_dist_l(i,ii)*1.D6)-5.D0
gal_appmag_g(i,ii)=g_g(i,ii)+5.D0*LOG10(gal_dist_l(i,ii)*1.D6)-5.D0

gal_appmag_J(i,ii)=g_J(i,ii)+5.D0*LOG10(gal_dist_l(i,ii)*1.D6)-5.D0
gal_appmag_Ks(i,ii)=g_Ks(i,ii)+5.D0*LOG10(gal_dist_l(i,ii)*1.D6)-5.D0
ELSE
WRITE(*,*) 'error: negative distance:',uuu,i,ii
END IF

END IF
END DO
END DO
WRITE(*,*) 'appmag calculated',uuu
 
 DO ii=1,6
 
 DO i=1,n_gal(ii)
  IF (gal_visible(i,ii)) THEN
!perform Anti-K-correction
K_corr=0.D0
DO iiii=0,5
DO iii=0,3
K_corr=K_corr+invK_coeff_r(iiii,iii)*(gal_redshift(i,ii)**iiii)*((g_g(i,ii)-g_r(i,ii))**iii)
END DO
END DO
!  WRITE(*,*) 'Kcorr r ok',uuu,i,ii
! IF ((K_corr>1.D0).OR.(K_corr<-1.D0)) THEN
! WRITE(*,*) 'Error',K_corr,'at',uuu,i,ii
! END IF
! invK_coeff_r


!ANTI-K-Correction 
gal_appmag_r(i,ii)=gal_appmag_r(i,ii)+K_corr




K_corr=0.D0
DO iiii=0,5
DO iii=0,3
K_corr=K_corr+invK_coeff_g(iiii,iii)*(gal_redshift(i,ii)**iiii)*((g_g(i,ii)-g_r(i,ii))**iii)
END DO
END DO


!ANTI-K-Correction 
gal_appmag_g(i,ii)=gal_appmag_g(i,ii)+K_corr
!  WRITE(*,*) 'Kcorr g ok',uuu,i,ii
! IF ((K_corr>1.D0).OR.(K_corr<-1.D0)) THEN
! WRITE(*,*) 'Error',K_corr,'at',uuu,i,ii
! END IF


K_corr=0.D0
DO iiii=0,5
DO iii=0,3
K_corr=K_corr+invK_coeff_J(iiii,iii)*(g_redshift2mass(i,ii)**iiii)*((g_J(i,ii)-g_Ks(i,ii))**iii)
END DO
END DO

!ANTI-K-Correction
gal_appmag_J(i,ii)=gal_appmag_J(i,ii)+K_corr

! IF ((K_corr>1.D0).OR.(K_corr<-1.D0)) THEN
! WRITE(*,*) 'Error',K_corr,'at',uuu,i,ii
! END IF



K_corr=0.D0
DO iiii=0,5
DO iii=0,3
K_corr=K_corr+invK_coeff_Ks(iiii,iii)*(g_redshift2mass(i,ii)**iiii)*((g_J(i,ii)-g_Ks(i,ii))**iii)
END DO
END DO

!ANTI-K-Correction
gal_appmag_Ks(i,ii)=gal_appmag_Ks(i,ii)+K_corr

! IF ((K_corr>1.D0).OR.(K_corr<-1.D0)) THEN
! WRITE(*,*) 'Error',K_corr,'at',uuu,i,ii
 END IF
END DO
END DO
WRITE(*,*) 'Kcordone',uuu

! 
rot_anglez=-45.0*PI/180.0
rot_angley=-20.0*PI/180.0


  DO ii=1,6
 
 DO i=1,n_gal(ii)
  IF (gal_visible(i,ii)) THEN
  !ANTI-extinction

  
rot_x=gx(i,ii)
rot_y=gy(i,ii)
rot_z=gz(i,ii)

rx=COS(rot_anglez)*rot_x-SIN(rot_anglez)*rot_y
ry=SIN(rot_anglez)*rot_x+COS(rot_anglez)*rot_y
rz=rot_z

rot_x=rx
rot_y=ry
rot_z=rz

rx=COS(rot_angley)*rot_x+SIN(rot_angley)*rot_z
ry=rot_y
rz=-SIN(rot_angley)*rot_x+COS(rot_angley)*rot_z

rot_x=rx
rot_y=ry
rot_z=rz

rot_dec=ASIN(rot_z/gal_dist_c(i,ii))*180.0/PI
rot_ra=ATAN2(rot_y,rot_x)*180.0/PI

rot_dec=90.0-rot_dec


proj_x=rot_dec*COS(rot_ra*PI/180.0)/90.0
proj_y=rot_dec*SIN(rot_ra*PI/180.0)/90.0

proj_x=(proj_x+1.0)*1024.0
proj_y=(proj_y+1.0)*1024.0

  

IF ((proj_x>0.D0).AND.(proj_x<2050.D0).AND.(proj_y>0.D0).AND.(proj_y<2050.D0)) THEN
g_fakeext_r(i,ii)=red_coeff_r*red_map(NINT(proj_x),NINT(proj_y))
g_fakeext_g(i,ii)=red_coeff_g*red_map(NINT(proj_x),NINT(proj_y))
ELSE
g_fakeext_r(i,ii)=0.D0
g_fakeext_g(i,ii)=0.D0
WRITE(*,*) 'outside area'
END IF



!av_galext_r+sigma_galext_r*gasdev(0)
gal_appmag_r(i,ii)=gal_appmag_r(i,ii)+g_fakeext_r(i,ii)
 
 !av_galext_g+sigma_galext_g*gasdev(0)
gal_appmag_g(i,ii)=gal_appmag_g(i,ii)+g_fakeext_g(i,ii)

 END IF
 END DO
END DO

WRITE(*,*) 'extinction done',uuu

 

 DO ii=1,6
 
 DO i=1,n_gal(ii)
 IF (gal_visible(i,ii)) THEN
!scatter redshift measurment with error
gal_redshift(i,ii)=gal_redshift(i,ii)+err_redshift/light*gasdev(0)
g_redshift2mass(i,ii)=g_redshift2mass(i,ii)+av_red_2mass*gasdev(0)
END IF
END DO
END DO

 

 
 WRITE(*,*) 'redshift error done',uuu
 
 
 DO ii=1,6
  DO i=1,n_gal(ii)
   IF (gal_visible(i,ii)) THEN
!photometric error
gal_appmag_r(i,ii)=gal_appmag_r(i,ii)+err_photo_r*gasdev(0)
gal_appmag_g(i,ii)=gal_appmag_g(i,ii)+err_photo_g*gasdev(0)

gal_appmag_J(i,ii)=gal_appmag_J(i,ii)+av_J_err*gasdev(0)
gal_appmag_Ks(i,ii)=gal_appmag_Ks(i,ii)+av_Ks_err*gasdev(0)
END IF
END DO
END DO




 DO ii=1,6
 
 DO i=1,n_gal(ii)
IF (gal_dist_c(i,ii)<0.1D0) THEN
gal_visible(i,ii)=.false.
END IF
END DO 
END DO

WRITE(*,*) 'photometric error done',uuu
 
 

 DO ii=1,6
 
 DO i=1,n_gal(ii)
 IF (gal_visible(i,ii)) THEN
galpos_ra(i,ii)=ATAN2(gy(i,ii),gx(i,ii))
galpos_dec(i,ii)=ASIN(gz(i,ii)/gal_dist_c(i,ii))
END IF
END DO
END DO
WRITE(*,*) 'pos done',uuu
 DO ii=1,6
 
 DO i=1,n_gal(ii)
 IF (gal_visible(i,ii)) THEN
rand_angle=2.D0*PI*RAND()
rand_value=err_pos*gasdev(0)


rand_dec=rand_value*COS(rand_angle)
rand_ra=rand_value*SIN(rand_angle)


galpos_ra(i,ii)=galpos_ra(i,ii)+rand_ra*COS(galpos_dec(i,ii))
galpos_dec(i,ii)=galpos_dec(i,ii)+rand_dec
END IF
END DO
END DO

WRITE(*,*) 'before filter',uuu
 DO ii=1,6
 
 DO i=1,n_gal(ii)
 IF (gal_visible(i,ii)) THEN

IF (galpos_dec(i,ii)>(PI/2.D0)) THEN
galpos_dec(i,ii)=(PI/2.D0)-(galpos_dec(i,ii)-(PI/2.D0))
galpos_ra(i,ii)=galpos_ra(i,ii)+PI
END IF
IF (galpos_ra(i,ii)>(2.D0*PI)) THEN
galpos_ra(i,ii)=galpos_ra(i,ii)-(2.D0*PI)
END IF
IF (galpos_ra(i,ii)<(0.D0)) THEN
galpos_ra(i,ii)=galpos_ra(i,ii)+(2.D0*PI)
END IF

galpos_ra(i,ii)=galpos_ra(i,ii)*180.D0/PI
galpos_dec(i,ii)=galpos_dec(i,ii)*180.D0/PI
END IF
  END DO
   
 END DO
 

  WRITE(*,*) 'observed values calculated',uuu
 
 ! use the right cube for the evolution
 DO ii=1,6
 

 
 zmax=(snapshot_redshift(ii)+snapshot_redshift(ii+1))/2.D0
 zmin=(snapshot_redshift(ii)+snapshot_redshift(ii-1))/2.D0

 DO i=1,n_gal(ii)
 IF (gal_visible(i,ii)) THEN
  
 IF (gal_cosmored(i,ii)>zmax) THEN
 gal_visible(i,ii)=.FALSE.
 END IF 
 
 IF (gal_cosmored(i,ii)<zmin) THEN
 gal_visible(i,ii)=.FALSE.
 END IF   
 
 END IF
 gal2_visible(i,ii)=gal_visible(i,ii)
 
  END DO
 
 END DO
 
 
  WRITE(*,*) 'evolution',uuu
 
  DO ii=1,6
 
 helpcount(ii)=0
DO i=1,n_gal(ii)
IF (gal_visible(i,ii)) THEN
 helpcount(ii)= helpcount(ii)+1
END IF
END DO
END DO 

 
 DO ii=1,6
 DO i=1,n_gal(ii)
IF (limitmag_SDSS_official<gal_appmag_r(i,ii)) THEN
 gal_visible(i,ii)=.FALSE.
 END IF
 IF (0.11D0<gal_redshift(i,ii)) THEN
 gal_visible(i,ii)=.FALSE.
 END IF
 
 IF (galpos_ra(i,ii)<0.D0) THEN
 gal_visible(i,ii)=.FALSE. 
 END IF
 IF (galpos_ra(i,ii)>90.D0) THEN
 gal_visible(i,ii)=.FALSE. 
 END IF
 IF (galpos_dec(i,ii)<0.D0) THEN
 gal_visible(i,ii)=.FALSE. 
 END IF

 
 END DO
 END DO
 
 
 
  DO ii=1,6
 DO i=1,n_gal(ii)
 
  IF (galpos_ra(i,ii)<0.D0) THEN
 gal2_visible(i,ii)=.FALSE. 
 END IF
 IF (galpos_ra(i,ii)>90.D0) THEN
 gal2_visible(i,ii)=.FALSE. 
 END IF
 IF (galpos_dec(i,ii)<0.D0) THEN
 gal2_visible(i,ii)=.FALSE. 
 END IF
 
 
 IF (limitmag_2MRS_official<gal_appmag_Ks(i,ii)) THEN
 gal2_visible(i,ii)=.FALSE.
 END IF
 
  END DO
 END DO
 
 
 
  WRITE(*,*) 'biases ok',uuu
 
 
 
 
 DO ii=1,6


DO i=1,n_gal(ii)
IF (gal2_visible(i,ii)) THEN
IF (RAND()>completeness_2mrs) THEN
 gal2_visible(i,ii)=.FALSE. 
END IF
END IF
END DO


END DO
 
 
!  
!   DO ii=1,6
! DO i=1,n_gal(ii)
! IF (gx(i,ii)>340.D0) THEN
!  gal_visible(i,ii)=.FALSE. 
!  gal2_visible(i,ii)=.FALSE. 
! END IF 
! IF (gy(i,ii)>340.D0) THEN
!  gal_visible(i,ii)=.FALSE. 
!  gal2_visible(i,ii)=.FALSE. 
! END IF 
! IF (gz(i,ii)>340.D0) THEN
!  gal_visible(i,ii)=.FALSE. 
!  gal2_visible(i,ii)=.FALSE. 
! END IF 
! END DO
!  END DO 
 
 
 
 
 DO ii=1,6
 
 helpcount(ii)=0
DO i=1,n_gal(ii)
IF (gal_visible(i,ii)) THEN
 helpcount(ii)= helpcount(ii)+1
END IF
END DO
END DO 



DO ii=1,6


DO i=1,n_gal(ii)
IF (gal_visible(i,ii)) THEN
IF (RAND()>decollided_sampling) THEN
 gal_visible(i,ii)=.FALSE. 
END IF
END IF
END DO


END DO


  n_unified=0
  DO ii=1,6
DO i=1,n_gal(ii)
IF (gal_visible(i,ii)) THEN
n_unified=n_unified+1
END IF
END DO
END DO



  n_unified2=0
  DO ii=1,6
DO i=1,n_gal(ii)
IF (gal2_visible(i,ii)) THEN
n_unified2=n_unified2+1
END IF
END DO
END DO


allocate(mockid2(1:n_unified2))
allocate(mmockid2(1:n_unified2))
allocate(mockfofid2(1:n_unified2))
allocate(mock2_ra(1:n_unified2))
allocate(mock2_dec(1:n_unified2))
allocate(mock2_z(1:n_unified2))
allocate(mock2_mJ(1:n_unified2))
allocate(mock2_mKs(1:n_unified2))
allocate(mock2_visible(1:n_unified2))


allocate(mockid(1:n_unified))
allocate(mmockid(1:n_unified))
allocate(mockfofid(1:n_unified))
allocate(mock_ra(1:n_unified))
allocate(mock_dec(1:n_unified))
allocate(mock_z(1:n_unified))
allocate(mock_mg(1:n_unified))
allocate(mock_mr(1:n_unified))
allocate(mock_visible(1:n_unified))
allocate(mock_ext_g(1:n_unified))
allocate(mock_ext_r(1:n_unified))


 
 
 
 
allocate(Kmock_g(1:n_unified))
allocate(Kmock_r(1:n_unified))
allocate(Kmock_J(1:n_unified2))
allocate(Kmock_Ks(1:n_unified2))
allocate(mockabsmag_g(1:n_unified))
allocate(mockabsmag_r(1:n_unified))
allocate(mockabsmag_J(1:n_unified2))
allocate(mockabsmag_Ks(1:n_unified2))
allocate(mockdist(1:n_unified))
allocate(mock2dist(1:n_unified2))
 
 allocate(mockxpos(1:n_unified))
allocate(mockypos(1:n_unified))
allocate(mockzpos(1:n_unified))

allocate(mock2xpos(1:n_unified2))
allocate(mock2ypos(1:n_unified2))
allocate(mock2zpos(1:n_unified2))



hcount=0
DO ii=1,6
DO i=1,n_gal(ii)
IF (gal_visible(i,ii)) THEN

hcount=hcount+1

mockfofid(hcount)=g_fofID(i,ii)
mockid(hcount)=galid(i,ii)
mmockid(hcount)=mgalid(i,ii)
mock_ra(hcount)=galpos_ra(i,ii)
mock_dec(hcount)=galpos_dec(i,ii)
mock_z(hcount)=gal_redshift(i,ii)
mock_mg(hcount)=gal_appmag_g(i,ii)
mock_mr(hcount)=gal_appmag_r(i,ii)
mock_visible(hcount)=.TRUE.
mock_ext_r(hcount)=g_fakeext_r(i,ii)
mock_ext_g(hcount)=g_fakeext_g(i,ii)
END IF
END DO
END DO


hcount=0
DO ii=1,6
DO i=1,n_gal(ii)
IF (gal2_visible(i,ii)) THEN

hcount=hcount+1

mockfofid2(hcount)=g_fofID(i,ii)
mockid2(hcount)=galid2(i,ii)
mmockid2(hcount)=mgalid2(i,ii)
mock2_ra(hcount)=galpos_ra(i,ii)
mock2_dec(hcount)=galpos_dec(i,ii)
mock2_z(hcount)=g_redshift2mass(i,ii)
mock2_mJ(hcount)=gal_appmag_J(i,ii)
mock2_mKs(hcount)=gal_appmag_Ks(i,ii)
mock2_visible(hcount)=.TRUE.

END IF
END DO
END DO





DO i=1,(n_unified-1)
IF (mock_visible(i)) THEN

DO ii=(i+1),n_unified
IF (mock_visible(ii)) THEN

angular_sep=COS(mock_dec(i)*PI/180.D0)*COS(mock_dec(ii)*PI/180.D0)*COS((mock_ra(ii)-mock_ra(i))*PI/180.D0)
angular_sep=angular_sep+SIN(mock_dec(i)*PI/180.D0)*SIN(mock_dec(ii)*PI/180.D0)
angular_sep=ACOS(angular_sep)*180.D0/PI

IF (angular_sep<0.D0) THEN
angular_sep=-angular_sep
END IF



IF (angular_sep<fiber_coll_r)THEN
 mock_visible(ii)=.FALSE. 
END IF


END IF
END DO
END IF
END DO

 
 
 
 
 
 
 
 
hcount=0 
  DO ii=1,6
 hcount=hcount+helpcount(ii)
END DO 

 
  active_count=0
DO i=1,n_unified
IF (mock_visible(i)) THEN
active_count=active_count+1
END IF
END DO

 loss_coll=DBLE(active_count)/DBLE(hcount)

  
  ! Correction for galaxies, which are in the sample dispite fiber collision due to overlapping plates
  help_ratio=(final_sampling-loss_coll)/(decollided_sampling-loss_coll)

DO i=1,n_unified
IF (mock_visible(i).EQV..FALSE.) THEN
IF (RAND()<help_ratio) THEN
mock_visible(i)=.TRUE. 
END IF
END IF
END DO

    active_count=0
DO i=1,n_unified
IF (mock_visible(i)) THEN
active_count=active_count+1
END IF
END DO

 loss_coll=DBLE(active_count)/DBLE(hcount)


  
  ! klick out QSO
  qso_ratio=(DBLE(n_qso)*mock_cover/sdss_cover)/DBLE(active_count)
  
  
DO i=1,n_unified
IF (mock_visible(i)) THEN
IF (RAND()<qso_ratio) THEN
mock_visible(i)=.FALSE. 
END IF
END IF
END DO
  
  
  
  
  
    active_count=0
DO i=1,n_unified
IF (mock_visible(i)) THEN
active_count=active_count+1
END IF
END DO

  
  n_mock_final=active_count
  
      active_count=0
DO i=1,n_unified2
IF (mock2_visible(i)) THEN
active_count=active_count+1
END IF
END DO

  
  n_mock_final2=active_count
  
 WRITE(*,*) '-------------------------------------------------------------------' 
 WRITE(*,*) n_mock_final,'galaxies in the SDSS prototype mock catalogue',uuu
 WRITE(34,*) '-------------------------------------------------------------------' 
 WRITE(34,*) n_mock_final,'galaxies in the SDSS prototype mock catalogue',uuu
 
ratio=(DBLE(n_sdss_used)/sdss_cover)/(DBLE(n_mock_final)/mock_cover)
 
 WRITE(*,*) 'completeness:',(1.D0/ratio)
 WRITE(*,*) '-------------------------------------------------------------------'
  WRITE(34,*) 'completeness:',(1.D0/ratio)
 WRITE(34,*) '-------------------------------------------------------------------'
 
 
  WRITE(*,*) n_mock_final2,'galaxies in the 2MRS prototype mock catalogue',uuu
    WRITE(34,*) n_mock_final2,'galaxies in the 2MRS prototype mock catalogue',uuu
  ratio=(DBLE(n_2mass_used)/twomass_cover)/(DBLE(n_mock_final2)/mock_cover)
 
 WRITE(*,*) 'completeness:',(1.D0/ratio)
 WRITE(*,*) '-------------------------------------------------------------------'
  WRITE(34,*) 'completeness:',(1.D0/ratio)
 WRITE(34,*) '-------------------------------------------------------------------'
 

 
 
 
 
 
 
 
 
 
 ! saturation limit
 DO i=1,n_unified
IF (mock_visible(i)) THEN

IF (mock_mg(i)<saturation_g) THEN
mock_visible(i)=.FALSE.
END IF

IF (mock_mr(i)<saturation_r) THEN
mock_visible(i)=.FALSE.
END IF

END IF
END DO
 
 
 
 
 
!correct for galactic extinction
DO i=1,n_unified
IF (mock_visible(i)) THEN
mock_mg(i)=mock_mg(i)-mock_ext_g(i)
mock_mr(i)=mock_mr(i)-mock_ext_r(i)
END IF
END DO


!K-correction
DO i=1,n_unified
IF (mock_visible(i)) THEN
Kmock_g(i)=0.D0
DO ii=0,7
DO iii=0,3
Kmock_g(i)=Kmock_g(i)+K_coeff_g(ii,iii)*(mock_z(i)**ii)*((mock_mg(i)-mock_mr(i))**iii)
END DO
END DO

Kmock_r(i)=0.D0
DO ii=0,5
DO iii=0,3
Kmock_r(i)=Kmock_r(i)+K_coeff_r(ii,iii)*(mock_z(i)**ii)*((mock_mg(i)-mock_mr(i))**iii)
END DO
END DO


END IF
END DO 



!K-correction
DO i=1,n_unified2
IF (mock2_visible(i)) THEN
Kmock_J(i)=0.D0
DO ii=0,5
DO iii=0,3
Kmock_J(i)=Kmock_J(i)+K_coeff_J(ii,iii)*(mock2_z(i)**ii)*((mock2_mJ(i)-mock2_mKs(i))**iii)
END DO
END DO


Kmock_Ks(i)=0.D0
DO ii=0,5
DO iii=0,3
Kmock_Ks(i)=Kmock_Ks(i)+K_coeff_Ks(ii,iii)*(mock2_z(i)**ii)*((mock2_mJ(i)-mock2_mKs(i))**iii)
END DO
END DO

END IF
END DO


! corrected apparent magnitudes
DO i=1,n_unified
IF (mock_visible(i)) THEN
mock_mg(i)=mock_mg(i)-Kmock_g(i) 
mock_mr(i)=mock_mr(i)-Kmock_r(i) 

END IF
END DO



! corrected apparent magnitudes
DO i=1,n_unified2
IF (mock2_visible(i)) THEN

mock2_mJ(i)=mock2_mJ(i)-Kmock_J(i) 

mock2_mKs(i)=mock2_mKs(i)-Kmock_Ks(i) 

END IF
END DO


DO i=1,n_unified
IF (mock_visible(i)) THEN
IF (mock_z(i)>0.D0) THEN

divisor=(SQRT(1.D0+2.D0*q0*mock_z(i))+1.D0+q0*mock_z(i))
!luminosity distance
mockdist(i)=light/H0*mock_z(i)*(1.D0+((mock_z(i)*(1.D0-q0))/divisor))
!distance in kpc
mockdist(i)=mockdist(i)*1000.D0

ELSE
mock_visible(i)=.FALSE.
mock_z(i)=10000.D0
mockdist(i)=1000000.D0
END IF
END IF
END DO


DO i=1,n_unified
IF (mock_visible(i)) THEN

mockabsmag_g(i)=mock_mg(i)-5.D0*LOG10(mockdist(i)*1000.D0)+5.D0
mockabsmag_r(i)=mock_mr(i)-5.D0*LOG10(mockdist(i)*1000.D0)+5.D0

END IF
END DO


DO i=1,n_unified
IF (mock_visible(i)) THEN
IF (mockabsmag_r(i)>-15.D0) THEN
mock_visible(i)=.FALSE.
END IF


IF (mock_z(i)>0.11D0) THEN
mock_visible(i)=.FALSE.
END IF
END IF
END DO



DO i=1,n_unified2
IF (mock2_visible(i)) THEN
IF (mock2_z(i)>0.D0) THEN

divisor=(SQRT(1.D0+2.D0*q0*mock2_z(i))+1.D0+q0*mock2_z(i))
!luminosity distance
mock2dist(i)=light/H0*mock2_z(i)*(1.D0+((mock2_z(i)*(1.D0-q0))/divisor))
!distance in kpc
mock2dist(i)=mock2dist(i)*1000.D0

ELSE
mock2_visible(i)=.FALSE.
mock2_z(i)=10000.D0
mock2dist(i)=1000000.D0
END IF
END IF
END DO

DO i=1,n_unified2
IF (mock2_visible(i)) THEN

mockabsmag_J(i)=mock2_mJ(i)-5.D0*LOG10(mock2dist(i)*1000.D0)+5.D0
mockabsmag_Ks(i)=mock2_mKs(i)-5.D0*LOG10(mock2dist(i)*1000.D0)+5.D0

END IF
END DO


DO i=1,n_unified
IF (mock_visible(i)) THEN
mockdist(i)=mockdist(i)*((1.D0+mock_z(i))**(-1.D0))
END IF
END DO

DO i=1,n_unified2
IF (mock2_visible(i)) THEN
mock2dist(i)=mock2dist(i)*((1.D0+mock2_z(i))**(-1.D0))
END IF
END DO


DO i=1,n_unified
IF (mock_visible(i)) THEN

IF (mock_mr(i)>(limitmag_SDSS_official+0.5D0)) THEN
mock_visible(i)=.FALSE.
END IF
END IF
END DO

DO i=1,n_unified2
IF (mock2_visible(i)) THEN

IF (mock2_mKs(i)>(limitmag_2MRS_official+0.5D0)) THEN
mock2_visible(i)=.FALSE.
END IF


END IF
END DO






 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
  
    active_count=0
DO i=1,n_unified
IF (mock_visible(i)) THEN
active_count=active_count+1
END IF
END DO

  
  n_mock_final=active_count
  
      active_count=0
DO i=1,n_unified2
IF (mock2_visible(i)) THEN
active_count=active_count+1
END IF
END DO

  
  n_mock_final2=active_count
  
 
 WRITE(*,*) n_mock_final,'galaxies in the SDSS final mock catalogue',uuu
 WRITE(33,*) n_mock_final,'galaxies in the SDSS final mock catalogue',uuu
 WRITE(34,*) n_mock_final,'galaxies in the SDSS final mock catalogue',uuu
 
 ratio=(DBLE(n_sdss_used)/sdss_cover)/(DBLE(n_mock_final)/mock_cover)
 
 WRITE(*,*) 'completeness:',(1.D0/ratio)
 WRITE(*,*) '-------------------------------------------------------------------'
 WRITE(34,*) 'completeness:',(1.D0/ratio)
 WRITE(34,*) '-------------------------------------------------------------------'

 
  WRITE(*,*) n_mock_final2,'galaxies in the 2MRS final mock catalogue',uuu
    WRITE(33,*) n_mock_final2,'galaxies in the 2MRS final mock catalogue',uuu
        WRITE(34,*) n_mock_final2,'galaxies in the 2MRS final mock catalogue',uuu
  ratio=(DBLE(n_2mass_used)/twomass_cover)/(DBLE(n_mock_final2)/mock_cover)
 
 WRITE(*,*) 'completeness:',(1.D0/ratio)
 WRITE(*,*) '-------------------------------------------------------------------'
 WRITE(34,*) 'completeness:',(1.D0/ratio)
 WRITE(34,*) '-------------------------------------------------------------------'
  
 
 
 

 
 
 
 
 
 DO i=1,n_unified
IF (mock_visible(i)) THEN
gal_b=mock_dec(i)*PI/180.D0
gal_l=mock_ra(i)*PI/180.D0
mockxpos(i)=mockdist(i)*COS(gal_b)*COS(gal_l)
mockypos(i)=mockdist(i)*COS(gal_b)*SIN(gal_l)
mockzpos(i)=mockdist(i)*SIN(gal_b)
END IF
END DO

DO i=1,n_unified2
IF (mock2_visible(i)) THEN
gal_b=mock2_dec(i)*PI/180.D0
gal_l=mock2_ra(i)*PI/180.D0
mock2xpos(i)=mock2dist(i)*COS(gal_b)*COS(gal_l)
mock2ypos(i)=mock2dist(i)*COS(gal_b)*SIN(gal_l)
mock2zpos(i)=mock2dist(i)*SIN(gal_b)
END IF
END DO

 
 
 
 

OPEN(50,file='SDSS_mock'//TRIM(number_mockset)//'/SDSS_galaxies_mock'//TRIM(number_mockset)//'_final.txt')
DO i=1,n_unified
IF (mock_visible(i)) THEN
WRITE(50,*) mock_ra(i),mock_dec(i),mock_z(i),mockabsmag_g(i),mockabsmag_r(i)
END IF
END DO
 CLOSE(50)
 
 
 
OPEN(50,file='2MRS_mock'//TRIM(number_mockset)//'/2MRS_galaxies_mock'//TRIM(number_mockset)//'_final.txt')
DO i=1,n_unified2
IF (mock2_visible(i)) THEN
WRITE(50,*) mock2_ra(i),mock2_dec(i),mock2_z(i),mockabsmag_J(i),mockabsmag_Ks(i)
END IF
END DO
 CLOSE(50)
 

 
  OPEN(50,file='SDSS_mock'//TRIM(number_mockset)//'/SDSS_galaxies_mock'//TRIM(number_mockset)//'_final_simid.txt')
DO i=1,n_unified
IF (mock_visible(i)) THEN
WRITE(50,*) mmockid(i),mockid(i),mockfofid(i)
END IF
END DO
 CLOSE(50)
 
! read file
OPEN(50,file='2MRS_mock'//TRIM(number_mockset)//'/2MRS_galaxies_mock'//TRIM(number_mockset)//'_final_simid.txt')
DO i=1,n_unified2
IF (mock2_visible(i)) THEN
WRITE(50,*) mmockid2(i),mockid2(i),mockfofid2(i)
END IF
END DO
 CLOSE(50)
 

deallocate(mock2_ra)
deallocate(mock2_dec)
deallocate(mock2_z)
deallocate(mock2_mJ)
deallocate(mock2_mKs)
deallocate(mock2_visible)
deallocate(mockfofid2)

deallocate(mock_ra)
deallocate(mock_dec)
deallocate(mock_z)
deallocate(mock_mg)
deallocate(mock_mr)
deallocate(mock_visible)
deallocate(mock_ext_g)
deallocate(mock_ext_r)

deallocate(mockid2)
deallocate(mockid)
deallocate(mockfofid)
deallocate(mmockid2)
deallocate(mmockid)
 
deallocate(Kmock_g)
deallocate(Kmock_r)
deallocate(Kmock_J)
deallocate(Kmock_Ks)
deallocate(mockabsmag_g)
deallocate(mockabsmag_r)
deallocate(mockabsmag_J)
deallocate(mockabsmag_Ks)
deallocate(mockdist)
deallocate(mock2dist)
 
deallocate(mockxpos)
deallocate(mockypos)
deallocate(mockzpos)

deallocate(mock2xpos)
deallocate(mock2ypos)
deallocate(mock2zpos)

 
 
 
    

 
 
 
 ! Halos

 

! 000
IF (uuu==1) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=hz(i,ii)
 END DO
 END DO 
END IF

!100
IF (uuu==2) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=hz(i,ii)
 END DO
 END DO 
END IF

!010
IF (uuu==3) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
 
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=hz(i,ii)
 
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=hz(i,ii)
 END DO
 END DO 
END IF

!001
IF (uuu==4) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
 
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=hz(i,ii)
 
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii)
 END DO
 END DO 
END IF

!110
IF (uuu==5) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
 
 
 
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii)
 
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=hz(i,ii)
 END DO
 END DO 
END IF

!011
IF (uuu==6) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
 
 
 
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=hz(i,ii)
 
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii)
 END DO
 END DO 
END IF

!101
IF (uuu==7) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)
 
 
  hx(i,ii)=hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii)
 
 
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii)
 END DO
 END DO 
END IF

!111
IF (uuu==8) THEN
 DO ii=1,6
 DO i=1,n_halo(ii)

 
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii) 
 
 
  hx(i,ii)=(usedsidelength/h)-hx(i,ii)
  hy(i,ii)=(usedsidelength/h)-hy(i,ii)
  hz(i,ii)=(usedsidelength/h)-hz(i,ii)
 END DO
 END DO 
END IF


 
 
  DO ii=1,6
 mtot(ii)=0.D0
 partsum(ii)=0.D0
 
 DO i=1,n_halo(ii)
 partsum(ii)=partsum(ii)+hnpart(i,ii)
 hm(i,ii)=hnpart(i,ii)*masspart/h
 mtot(ii)=mtot(ii)+hm(i,ii)
 END DO

 m_ratio(ii)=mtot(ii)/m_expected
 
   END DO
  
  
  
    DO ii=1,6
   DO i=1,n_halo(ii)
   

  
  hdist(i,ii)=((hx(i,ii))**2)+((hy(i,ii))**2)+((hz(i,ii))**2)
hdist(i,ii)=SQRT(hdist(i,ii))

 
 !redshift
hredshift(i,ii)=(light**2)*((light**2)+(hdist(i,ii)**2)*(H0**2)*(1.D0-2.D0*q0))*(q0-1.D0)**2
hredshift(i,ii)=SQRT(hredshift(i,ii))
hredshift(i,ii)=hredshift(i,ii)+light*hdist(i,ii)*H0+(light**2)*(q0-1.D0)-(hdist(i,ii)**2)*(H0**2)*(q0**2)
hredshift(i,ii)=hredshift(i,ii)/((light-hdist(i,ii)*H0*q0)**2)

  
    END DO
    END DO
  
  ! use the right cube for the evolution
 DO ii=1,6
 

 
 zmax=(snapshot_redshift(ii)+snapshot_redshift(ii+1))/2.D0
 zmin=(snapshot_redshift(ii)+snapshot_redshift(ii-1))/2.D0

  DO i=1,n_halo(ii)
  halo_part_mock(i,ii)=.TRUE.
  

 IF (hredshift(i,ii)>zmax) THEN
 halo_part_mock(i,ii)=.FALSE.
 END IF 
 
 IF (hredshift(i,ii)<zmin) THEN
 halo_part_mock(i,ii)=.FALSE.
 END IF   
 
 
!  IF (hredshift(i,ii)>0.11D0) THEN
!  halo_part_mock(i,ii)=.FALSE.
!  END IF    
END DO
 
 END DO 


! DO ii=1,6
! DO i=1,n_halo(ii)
! IF (hx(i,ii)>340.D0) THEN
!  halo_part_mock(i,ii)=.FALSE.
! END IF 
! IF (hy(i,ii)>340.D0) THEN
!  halo_part_mock(i,ii)=.FALSE.
! END IF 
! IF (hz(i,ii)>340.D0) THEN
!  halo_part_mock(i,ii)=.FALSE.
! END IF 
! END DO
! END DO
!  
 
   n_hunified=0
  DO ii=1,6
DO i=1,n_halo(ii)
IF (halo_part_mock(i,ii)) THEN
n_hunified=n_hunified+1
END IF
END DO
END DO
 
 
 
 allocate(halo_x(1:n_hunified))
 allocate(halo_y(1:n_hunified))
 allocate(halo_z(1:n_hunified))
 allocate(halo_mass(1:n_hunified))
 allocate(halo_ra(1:n_hunified))
 allocate(halo_dec(1:n_hunified))
 allocate(halo_red(1:n_hunified))
  allocate(hfofid(1:n_hunified))
  allocate(hmock_nsg(1:n_hunified))
 
hcount=0
DO ii=1,6
DO i=1,n_halo(ii)
IF (halo_part_mock(i,ii)) THEN

hcount=hcount+1

halo_x(hcount)=hx(i,ii)
halo_y(hcount)=hy(i,ii)
halo_z(hcount)=hz(i,ii)
halo_mass(hcount)=hm(i,ii)
halo_red(hcount)=hredshift(i,ii)
hfofid(hcount)=h_fofId(i,ii)
hmock_nsg(hcount)=h_nsg(i,ii)

END IF
END DO
END DO

 
 
 
 DO i=1,n_hunified
 
divisor=((halo_x(i))**2)+((halo_y(i))**2)+((halo_z(i))**2)
divisor=SQRT(divisor)
 
 halo_ra(i)=ATAN2(halo_y(i),halo_x(i))
 halo_dec(i)=ASIN(halo_z(i)/divisor)



IF (halo_ra(i)>(2.D0*PI)) THEN
halo_ra(i)=halo_ra(i)-(2.D0*PI)
END IF
IF (halo_ra(i)<(0.D0)) THEN
halo_ra(i)=halo_ra(i)+(2.D0*PI)
END IF

halo_ra(i)=halo_ra(i)*180.D0/PI
halo_dec(i)=halo_dec(i)*180.D0/PI
 
 
 END DO
 
 
 
 
  WRITE(*,*) n_hunified,'halos in the halo mock catalogue',uuu

 WRITE(*,*) '-------------------------------------------------------------------'

  WRITE(34,*) n_hunified,'halos in the halo mock catalogue',uuu

 WRITE(34,*) '-------------------------------------------------------------------'
 
 
  ! read file
OPEN(50,file='halos'//TRIM(number_mockset)//'/halos_mock'//TRIM(number_mockset)//'_pos3D.txt')
DO i=1,n_hunified
WRITE(50,*) hfofid(i),halo_x(i),halo_y(i),halo_z(i),halo_mass(i),hmock_nsg(i)
END DO
 CLOSE(50)
 
 
!    ! read file
! OPEN(50,file='halos'//TRIM(number_mockset)//'/halos_mock'//TRIM(number_mockset)//'_radecz.txt')
! DO i=1,n_hunified
! WRITE(50,*) hfofid(i),halo_ra(i),halo_dec(i),halo_red(i),halo_mass(i),hmock_nsg(i)
! END DO
!  CLOSE(50)
!  

 
 
deallocate(halo_x)
deallocate(halo_y)
deallocate(halo_z)
deallocate(halo_mass)
deallocate(halo_ra)
deallocate(halo_dec)
deallocate(halo_red)
deallocate(hfofid)
deallocate(hmock_nsg)
 
 
 
 
 
 
 

  
    DO ii=1,6
   DO i=1,n_gal(ii)
   

  
  gdist(i,ii)=((gx(i,ii))**2)+((gy(i,ii))**2)+((gz(i,ii))**2)
gdist(i,ii)=SQRT(gdist(i,ii))

 
 !redshift
gredshift(i,ii)=(light**2)*((light**2)+(gdist(i,ii)**2)*(H0**2)*(1.D0-2.D0*q0))*(q0-1.D0)**2
gredshift(i,ii)=SQRT(gredshift(i,ii))
gredshift(i,ii)=gredshift(i,ii)+light*gdist(i,ii)*H0+(light**2)*(q0-1.D0)-(gdist(i,ii)**2)*(H0**2)*(q0**2)
gredshift(i,ii)=gredshift(i,ii)/((light-gdist(i,ii)*H0*q0)**2)

  
    END DO
    END DO
  
  ! use the right cube for the evolution
 DO ii=1,6
 

 
 zmax=(snapshot_redshift(ii)+snapshot_redshift(ii+1))/2.D0
 zmin=(snapshot_redshift(ii)+snapshot_redshift(ii-1))/2.D0

  DO i=1,n_gal(ii)
  gal_part_mock(i,ii)=.TRUE.
  

 IF (gredshift(i,ii)>zmax) THEN
 gal_part_mock(i,ii)=.FALSE.
 END IF 
 
 IF (gredshift(i,ii)<zmin) THEN
 gal_part_mock(i,ii)=.FALSE.
 END IF   
 
!  IF (gredshift(i,ii)>0.11D0) THEN
!  gal_part_mock(i,ii)=.FALSE.
!  END IF    
!  
  END DO
 
 END DO 

!    DO ii=1,6
! DO i=1,n_gal(ii)
! IF (gx(i,ii)>340.D0) THEN
!  gal_part_mock(i,ii)=.FALSE. 
! END IF 
! IF (gy(i,ii)>340.D0) THEN
!  gal_part_mock(i,ii)=.FALSE. 
! END IF 
! IF (gz(i,ii)>340.D0) THEN
!  gal_part_mock(i,ii)=.FALSE. 
! END IF 
! END DO
!  END DO 
 
   n_hunified=0
  DO ii=1,6
DO i=1,n_gal(ii)
IF (gal_part_mock(i,ii)) THEN
n_hunified=n_hunified+1
END IF
END DO
END DO
 
 
 
 allocate(halo_x(1:n_hunified))
 allocate(halo_y(1:n_hunified))
 allocate(halo_z(1:n_hunified))

 allocate(halo_ra(1:n_hunified))
 allocate(halo_dec(1:n_hunified))
 allocate(halo_red(1:n_hunified))
 
  allocate(g_mag_g(1:n_hunified))
  allocate(g_mag_r(1:n_hunified))
  allocate(g_mag_J(1:n_hunified))
  allocate(g_mag_Ks(1:n_hunified))
allocate(mockid(1:n_hunified))
allocate(mmockid(1:n_hunified))
  allocate(hfofid(1:n_hunified))
  
hcount=0
DO ii=1,6
DO i=1,n_gal(ii)
IF (gal_part_mock(i,ii)) THEN

hcount=hcount+1

mockid(hcount)=galid(i,ii)
mmockid(hcount)=mgalid(i,ii)
halo_x(hcount)=gx(i,ii)
halo_y(hcount)=gy(i,ii)
halo_z(hcount)=gz(i,ii)
g_mag_g(hcount)=g_g(i,ii)
g_mag_r(hcount)=g_r(i,ii)
g_mag_J(hcount)=g_J(i,ii)
g_mag_Ks(hcount)=g_Ks(i,ii)
halo_red(hcount)=hredshift(i,ii)
hfofid(hcount)=g_fofID(i,ii)
END IF
END DO
END DO

 
 
 
 DO i=1,n_hunified
 
divisor=((halo_x(i))**2)+((halo_y(i))**2)+((halo_z(i))**2)
divisor=SQRT(divisor)
 
 halo_ra(i)=ATAN2(halo_y(i),halo_x(i))
 halo_dec(i)=ASIN(halo_z(i)/divisor)



IF (halo_ra(i)>(2.D0*PI)) THEN
halo_ra(i)=halo_ra(i)-(2.D0*PI)
END IF
IF (halo_ra(i)<(0.D0)) THEN
halo_ra(i)=halo_ra(i)+(2.D0*PI)
END IF

halo_ra(i)=halo_ra(i)*180.D0/PI
halo_dec(i)=halo_dec(i)*180.D0/PI
 
 
 END DO
 
 
 
 
  WRITE(*,*) n_hunified,'galaxies in the all galaxies mock catalogue',uuu

 WRITE(*,*) '-------------------------------------------------------------------'

  WRITE(34,*) n_hunified,'galaxies in the all galaxies mock catalogue',uuu

 WRITE(34,*) '-------------------------------------------------------------------'
 
 
  ! read file
OPEN(50,file='2MRS_mock'//TRIM(number_mockset)//'/2MRS_allgalaxies_mock'//TRIM(number_mockset)//'_pos3D.txt')
DO i=1,n_hunified
WRITE(50,*) halo_x(i),halo_y(i),halo_z(i),g_mag_J(i),g_mag_Ks(i)
END DO
 CLOSE(50)
 
 
!    ! read file
! OPEN(50,file='2MRS_mock'//TRIM(number_mockset)//'/2MRS_allgalaxies_mock'//TRIM(number_mockset)//'_radecz.txt')
! DO i=1,n_hunified
! WRITE(50,*) halo_ra(i),halo_dec(i),halo_red(i),g_mag_J(i),g_mag_Ks(i)
! END DO
!  CLOSE(50)
 
OPEN(50,file='SDSS_mock'//TRIM(number_mockset)//'/SDSS_allgalaxies_mock'//TRIM(number_mockset)//'_pos3D.txt')
DO i=1,n_hunified
WRITE(50,*) halo_x(i),halo_y(i),halo_z(i),g_mag_g(i),g_mag_r(i)
END DO
 CLOSE(50)
 
 
!    ! read file
! OPEN(50,file='SDSS_mock'//TRIM(number_mockset)//'/SDSS_allgalaxies_mock'//TRIM(number_mockset)//'_radecz.txt')
! DO i=1,n_hunified
! WRITE(50,*) halo_ra(i),halo_dec(i),halo_red(i),g_mag_g(i),g_mag_r(i)
! END DO
!  CLOSE(50)
!  
 

OPEN(50,file='SDSS_mock'//TRIM(number_mockset)//'/SDSS_allgalaxies_mock'//TRIM(number_mockset)//'_id_file.txt')
DO i=1,n_hunified
WRITE(50,*) mmockid(i),mockid(i),hfofid(i)
END DO
 CLOSE(50)

! read file
OPEN(50,file='2MRS_mock'//TRIM(number_mockset)//'/2MRS_allgalaxies_mock'//TRIM(number_mockset)//'_id_file.txt')
DO i=1,n_hunified
WRITE(50,*) mmockid(i),mockid(i),hfofid(i)
END DO
 CLOSE(50)
 
deallocate(halo_x)
deallocate(halo_y)
deallocate(halo_z)

deallocate(halo_ra)
deallocate(halo_dec)
deallocate(halo_red)
 
  deallocate(g_mag_g)
  deallocate(g_mag_r)
  deallocate(g_mag_J)
  deallocate(g_mag_Ks)
 
 deallocate(mockid)
  deallocate(mmockid)
 
   deallocate(hfofid)
!  
! END DO


 END IF 
! call MPI_Scatter(positions,(3*npara),MPI_DOUBLE_PRECISION,positions_intern,(3*npara),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! call MPI_AllGather(inside_fi_intern,(npara),MPI_DOUBLE_PRECISION,help_fi_intern,(npara),MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

 CLOSE(34)


! call MPI_AllGather(u_local,1,MPI_INTEGER,u_array,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
! 
!  
! !finalize parallization
! call MPI_FINALIZE(ierr)
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
