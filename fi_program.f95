PROGRAM fi_regions
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer(kind=8) :: io_err,i,ii,iii,n_sdss,n_2mrs,help_i,n_all,n_overlap,n_overlap2,hc1,hc2
integer(kind=8) :: n_gal_sdss,n_gal_2mrs,limiting_funct,iiii,members_help,n_old,n_new,n_rand,n_rand_used,n_rand_fi
integer :: r_seed,n_mfi_final_max,n_mfi_first_max,evolution_seg,n_rand_multi
double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,dummy_var,av_rad_pec_vel
double precision :: help_mass_err,help_mass_err_weight
real :: out_px,out_py,out_pz,out_m,out_merr,out_r,out_rerr

double precision, allocatable :: ra(:),dec(:),red(:),mass(:),dist_c(:),px(:),py(:),pz(:),rad(:),mass_groups(:)
double precision, allocatable :: px_help,py_help,pz_help,mass_help,weights_help,weight_sum
double precision, allocatable :: ra_2mrs(:),dec_2mrs(:),z_2mrs(:),mass_2mrs(:),dist_c_2mrs(:)
double precision, allocatable :: ra_sdss(:),dec_sdss(:),z_sdss(:),mass_sdss(:),dist_c_sdss(:)
double precision, allocatable :: angrad_sdss(:),angrad_2mrs(:),angrad(:),sigma_sdss(:),sigma_2mrs(:),sigma(:)
double precision, allocatable :: m_halo_first(:,:),m_fi_first(:,:),m_halo_final(:,:),m_fi_final(:,:)
logical, allocatable :: active(:),sdss_in_2mrs(:),twomrs_in_sdss(:),gal_active_sdss(:)
logical, allocatable :: gal_active_2mrs(:),rand_act(:),rand_infi(:),rand_infi_multi(:)
double precision, allocatable :: gal_z_sdss(:),gal_mag_r(:),volume_weight_sdss(:),log_lum_dist_2mrs(:)
double precision, allocatable :: vol_int_sdss(:),log_lum_dist_sdss(:),v_cor_sdss(:),vol_sat_sdss(:)
double precision, allocatable :: gal_z_2mrs(:),gal_mag_K(:),volume_weight_2mrs(:),v_cor_2mrs(:)
double precision, allocatable :: b_eff(:),basic_link_a(:),R_eff(:),basic_link_R(:)
integer(kind=8), allocatable :: members_sdss(:),members_2mrs(:),members(:)
double precision, allocatable :: rand_x(:),rand_y(:),rand_z(:)
double precision, allocatable :: mass_sdss_err(:),mass_2mrs_err(:),mass_err(:),rad_err(:)
double precision, allocatable :: ra_fi(:),dec_fi(:),cdist_fi(:)


double precision :: sdss_cover,twomrs_cover,mock_cover,masspart,basic_link,d_help,d_limit_2mrs,Vol_fi
double precision :: mag_limit_2mrs,mabs_volsat,ts_modificator,d,mnew,logmh,logmfi,d_maxdepth,final_dens
double precision :: v_limit,v_sum,t1,t2,t3,D_limit,mag_limit_sdss,mag_saturation_sdss,mabs_vollim,z_limit
double precision :: b_factor,R_factor,lum_power,delta_z,angular_sep,sum_mass,V_inner,V_outer

double precision, dimension(1:3) :: Yvector
double precision, dimension(1:3,1:3) :: Amatrix
double precision, dimension(1:3,1:3,1:3) :: Bmatrices
double precision, dimension(1:3,1:6) :: solution_first,errorbar_first,solution_final,errorbar_final
double precision, allocatable :: datalist_first(:,:),datalist_final(:,:)
double precision, allocatable :: deviation_first(:),deviation_final(:)
double precision :: rms,det_A,otherside,gal_redshift,gal_dist_c,zmin,zmax
double precision, dimension(1:3) :: det_B
double precision, dimension(1:3,1:3) :: help_matrix
double precision, dimension(1:2,1:2) :: help_submatrix

integer, dimension(1:6) :: n_mfi_final,n_mfi_first
double precision, dimension (0:7) :: redshift_snap
double precision, dimension (1:6) :: rms_final,rms_first
character(1) :: numbersnap

double precision :: det2d,det3d

double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y
double precision, dimension(0:50) :: bin_mass_2mrs,bin_mass_sdss,bin_density_2mrs,bin_density_sdss,V_c
double precision, dimension(0:50) :: bin_mass_combi,bin_density_combi,bin_mass_fi,bin_density_fi


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme FINITE INFINITY REGIONS started'
WRITE(*,*) '============================================================'

 CALL SYSTEM_CLOCK(r_seed)
 r_seed=r_seed-INT(DBLE(r_seed)/100000.D0)*100000
 CALL srand(r_seed) 

  OPEN(50,file='globalparameters_SDSS.txt')
READ(50,*) basic_link
READ(50,*)
READ(50,*) av_rad_pec_vel
 CLOSE(50)
 
 
 OPEN(50,file='input_cf_2MRS.txt')

READ(50,*) 
READ(50,*) b_factor,R_factor,lum_power
READ(50,*)  
 CLOSE(50)

! define constants
PI=ACOS(-1.D0)

Omega_m=0.25D0
Omega_l=0.75D0
H0=73.D0

sdss_cover=9376.D0/41253.D0
twomrs_cover=0.91D0

q0=Omega_m/2.D0-Omega_l
light=3.D5
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)

ts_modificator=(50.1D0/61.7D0)**2
part_exp=(2160.D0**3)*((330.D0/500.D0)**3)
! rho_crit=3.D0*((H0)**2)/(8.D0*PI*(G*((1.D-6)**3)))
 cV_mill=(400.D6/h)**3
m_expected=cV_mill*rho_crit*Omega_m

! WRITE(*,*) m_expected
limiting_funct=1
mock_cover=1.D0/8.D0

WRITE(*,*) 'variables initalized'

! get length of file
OPEN(50,file='ts_sor/catalogues/cluster_list_all_SDSS.txt')
io_err=0
n_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss=n_sdss+1
END DO
 CLOSE(50)
n_sdss=n_sdss-1

! get length of file
OPEN(50,file='ts_sor/catalogues/cluster_list_all_2MRS.txt')
io_err=0
n_2mrs=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mrs=n_2mrs+1
END DO
 CLOSE(50)
n_2mrs=n_2mrs-1

! get length of file
OPEN(50,file='millimil/data/63/final_masses_fi.txt')
io_err=0
n_mfi_final(1)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_final(1)=n_mfi_final(1)+1
END DO
 CLOSE(50)
n_mfi_final(1)=n_mfi_final(1)-1

! get length of file
OPEN(50,file='millimil/data/63/first_masses_fi.txt')
io_err=0
n_mfi_first(1)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_first(1)=n_mfi_first(1)+1
END DO
 CLOSE(50)
n_mfi_first(1)=n_mfi_first(1)-1


! get length of file
OPEN(50,file='millimil/data/62/final_masses_fi.txt')
io_err=0
n_mfi_final(2)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_final(2)=n_mfi_final(2)+1
END DO
 CLOSE(50)
n_mfi_final(2)=n_mfi_final(2)-1

! get length of file
OPEN(50,file='millimil/data/62/first_masses_fi.txt')
io_err=0
n_mfi_first(2)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_first(2)=n_mfi_first(2)+1
END DO
 CLOSE(50)
n_mfi_first(2)=n_mfi_first(2)-1


! get length of file
OPEN(50,file='millimil/data/61/final_masses_fi.txt')
io_err=0
n_mfi_final(3)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_final(3)=n_mfi_final(3)+1
END DO
 CLOSE(50)
n_mfi_final(3)=n_mfi_final(3)-1

! get length of file
OPEN(50,file='millimil/data/61/first_masses_fi.txt')
io_err=0
n_mfi_first(3)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_first(3)=n_mfi_first(3)+1
END DO
 CLOSE(50)
n_mfi_first(3)=n_mfi_first(3)-1


! get length of file
OPEN(50,file='millimil/data/60/final_masses_fi.txt')
io_err=0
n_mfi_final(4)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_final(4)=n_mfi_final(4)+1
END DO
 CLOSE(50)
n_mfi_final(4)=n_mfi_final(4)-1

! get length of file
OPEN(50,file='millimil/data/60/first_masses_fi.txt')
io_err=0
n_mfi_first(4)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_first(4)=n_mfi_first(4)+1
END DO
 CLOSE(50)
n_mfi_first(4)=n_mfi_first(4)-1


! get length of file
OPEN(50,file='millimil/data/59/final_masses_fi.txt')
io_err=0
n_mfi_final(5)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_final(5)=n_mfi_final(5)+1
END DO
 CLOSE(50)
n_mfi_final(5)=n_mfi_final(5)-1

! get length of file
OPEN(50,file='millimil/data/59/first_masses_fi.txt')
io_err=0
n_mfi_first(5)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_first(5)=n_mfi_first(5)+1
END DO
 CLOSE(50)
n_mfi_first(5)=n_mfi_first(5)-1


! get length of file
OPEN(50,file='millimil/data/58/final_masses_fi.txt')
io_err=0
n_mfi_final(6)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_final(6)=n_mfi_final(6)+1
END DO
 CLOSE(50)
n_mfi_final(6)=n_mfi_final(6)-1


! get length of file
OPEN(50,file='millimil/data/58/first_masses_fi.txt')
io_err=0
n_mfi_first(6)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mfi_first(6)=n_mfi_first(6)+1
END DO
 CLOSE(50)
n_mfi_first(6)=n_mfi_first(6)-1




n_all=n_2mrs+n_sdss





WRITE(*,*) 'file sizes measured'



allocate(ra(1:n_all))
allocate(dec(1:n_all))
allocate(red(1:n_all))
allocate(mass(1:n_all))
allocate(dist_c(1:n_all))
allocate(mass_groups(1:n_all))
allocate(px(1:n_all))
allocate(py(1:n_all))
allocate(pz(1:n_all))
allocate(rad(1:n_all))
allocate(active(1:n_all))
allocate(members(1:n_all))
allocate(mass_err(1:n_all))
allocate(rad_err(1:n_all))
allocate(mass_sdss_err(1:n_sdss))
allocate(mass_2mrs_err(1:n_2mrs))
allocate(angrad_sdss(1:n_sdss))
allocate(angrad_2mrs(1:n_2mrs))
allocate(angrad(1:n_all))
allocate(sigma_sdss(1:n_sdss))
allocate(sigma_2mrs(1:n_2mrs))
allocate(sigma(1:n_all))
allocate(members_sdss(1:n_sdss))
allocate(members_2mrs(1:n_2mrs))
allocate(ra_2mrs(1:n_2mrs))
allocate(dec_2mrs(1:n_2mrs))
allocate(z_2mrs(1:n_2mrs))
allocate(mass_2mrs(1:n_2mrs))
allocate(ra_sdss(1:n_sdss))
allocate(dec_sdss(1:n_sdss))
allocate(z_sdss(1:n_sdss))
allocate(mass_sdss(1:n_sdss))
allocate(dist_c_2mrs(1:n_2mrs))
allocate(dist_c_sdss(1:n_sdss))

allocate(ra_fi(1:n_all))
allocate(dec_fi(1:n_all))
allocate(cdist_fi(1:n_all))



n_mfi_final_max=0
n_mfi_first_max=0
DO i=1,6
IF (n_mfi_final(i)>n_mfi_final_max) THEN
n_mfi_final_max=n_mfi_final(i)
END IF

IF (n_mfi_first(i)>n_mfi_first_max) THEN
n_mfi_first_max=n_mfi_first(i)
END IF
END DO




allocate(m_halo_first(1:n_mfi_first_max,1:6))
allocate(m_fi_first(1:n_mfi_first_max,1:6))
allocate(m_halo_final(1:n_mfi_final_max,1:6))
allocate(m_fi_final(1:n_mfi_final_max,1:6))

allocate(datalist_final(1:4,1:n_mfi_final_max))
allocate(deviation_final(1:n_mfi_final_max))
allocate(datalist_first(1:4,1:n_mfi_first_max))
allocate(deviation_first(1:n_mfi_first_max))

allocate(vol_int_sdss(1:n_sdss))
allocate(log_lum_dist_sdss(1:n_sdss))
allocate(v_cor_sdss(1:n_sdss))
allocate(vol_sat_sdss(1:n_sdss))

allocate(log_lum_dist_2mrs(1:n_2mrs))
allocate(v_cor_2mrs(1:n_2mrs))

allocate(b_eff(1:n_2mrs))
allocate(basic_link_a(1:n_2mrs))
allocate(R_eff(1:n_2mrs))
allocate(basic_link_R(1:n_2mrs))
allocate(sdss_in_2mrs(1:n_sdss))

WRITE(*,*) 'arrays allocated'




OPEN(50,file='snapshot_redshift_list.txt')
DO i=1,6
READ(50,*) ii,redshift_snap(i)
END DO
 CLOSE(50)
redshift_snap(0)=0.D0
redshift_snap(7)=1.D0

OPEN(50,file='ts_sor/catalogues/cluster_list_all_2MRS.txt')
DO i=1,n_2mrs
READ(50,*) ii,ra_2mrs(i),dec_2mrs(i),z_2mrs(i),&
dummy_var,dummy_var,dummy_var,mass_2mrs(i),mass_2mrs_err(i),dummy_var,&
sigma_2mrs(i),dummy_var,angrad_2mrs(i),dummy_var,members_2mrs(i)
END DO
 CLOSE(50)
 
 
 WRITE(*,*) n_2mrs,'clusters read in'



OPEN(50,file='ts_sor/catalogues/cluster_list_all_SDSS.txt')
DO i=1,n_sdss
READ(50,*) ii,ra_sdss(i),dec_sdss(i),z_sdss(i),&
dummy_var,dummy_var,dummy_var,mass_sdss(i),mass_sdss_err(i),dummy_var,&
sigma_sdss(i),dummy_var,angrad_sdss(i),dummy_var,members_sdss(i)
END DO
 CLOSE(50)

 WRITE(*,*) n_sdss,'clusters read in'
 
 OPEN(50,file='millimil/data/63/first_masses_fi.txt')
DO i=1,n_mfi_first(1)
READ(50,*)  m_halo_first(i,1),m_fi_first(i,1)
END DO
 CLOSE(50)
 WRITE(*,*) n_mfi_first(1),'masses read in'
 
 OPEN(50,file='millimil/data/63/final_masses_fi.txt')
DO i=1,n_mfi_final(1)
READ(50,*) m_halo_final(i,1),m_fi_final(i,1)
END DO
 CLOSE(50)
  WRITE(*,*) n_mfi_final(1),'masses read in'
 
 
  OPEN(50,file='millimil/data/62/first_masses_fi.txt')
DO i=1,n_mfi_first(2)
READ(50,*)  m_halo_first(i,2),m_fi_first(i,2)
END DO
 CLOSE(50)
 WRITE(*,*) n_mfi_first(2),'masses read in'
 
 OPEN(50,file='millimil/data/62/final_masses_fi.txt')
DO i=1,n_mfi_final(2)
READ(50,*) m_halo_final(i,2),m_fi_final(i,2)
END DO
 CLOSE(50)
  WRITE(*,*) n_mfi_final(2),'masses read in'
 
 
  OPEN(50,file='millimil/data/61/first_masses_fi.txt')
DO i=1,n_mfi_first(3)
READ(50,*)  m_halo_first(i,3),m_fi_first(i,3)
END DO
 CLOSE(50)
 WRITE(*,*) n_mfi_first(3),'masses read in'
 
 OPEN(50,file='millimil/data/61/final_masses_fi.txt')
DO i=1,n_mfi_final(3)
READ(50,*) m_halo_final(i,3),m_fi_final(i,3)
END DO
 CLOSE(50)
  WRITE(*,*) n_mfi_final(3),'masses read in'
 
 
  OPEN(50,file='millimil/data/60/first_masses_fi.txt')
DO i=1,n_mfi_first(4)
READ(50,*)  m_halo_first(i,4),m_fi_first(i,4)
END DO
 CLOSE(50)
 WRITE(*,*) n_mfi_first(4),'masses read in'
 
 OPEN(50,file='millimil/data/60/final_masses_fi.txt')
DO i=1,n_mfi_final(4)
READ(50,*) m_halo_final(i,4),m_fi_final(i,4)
END DO
 CLOSE(50)
  WRITE(*,*) n_mfi_final(4),'masses read in'
 
 
  OPEN(50,file='millimil/data/59/first_masses_fi.txt')
DO i=1,n_mfi_first(5)
READ(50,*)  m_halo_first(i,5),m_fi_first(i,5)
END DO
 CLOSE(50)
 WRITE(*,*) n_mfi_first(5),'masses read in'
 
 OPEN(50,file='millimil/data/59/final_masses_fi.txt')
DO i=1,n_mfi_final(5)
READ(50,*) m_halo_final(i,5),m_fi_final(i,5)
END DO
 CLOSE(50)
  WRITE(*,*) n_mfi_final(5),'masses read in'
 
 
  OPEN(50,file='millimil/data/58/first_masses_fi.txt')
DO i=1,n_mfi_first(6)
READ(50,*)  m_halo_first(i,6),m_fi_first(i,6)
END DO
 CLOSE(50)
 WRITE(*,*) n_mfi_first(6),'masses read in'
 
 OPEN(50,file='millimil/data/58/final_masses_fi.txt')
DO i=1,n_mfi_final(6)
READ(50,*) m_halo_final(i,6),m_fi_final(i,6)
END DO
 CLOSE(50)
  WRITE(*,*) n_mfi_final(6),'masses read in'
 
 
 
 
 
 
 
 
 WRITE(*,*) 'all data read in'
 
 
 
 
 DO i=1,n_2mrs
divisor=(SQRT(1.D0+2.D0*q0*z_2mrs(i))+1.D0+q0*z_2mrs(i))
dist_c_2mrs(i)=light/H0*z_2mrs(i)*(1.D0+((z_2mrs(i)*(1.D0-q0))/divisor))
log_lum_dist_2mrs(i)=LOG10(dist_c_2mrs(i)*1.D6)
dist_c_2mrs(i)=dist_c_2mrs(i)/(1.D0+z_2mrs(i))*1.D6
mass_2mrs(i)=(10.D0**mass_2mrs(i))
END DO
 
  DO i=1,n_sdss
divisor=(SQRT(1.D0+2.D0*q0*z_sdss(i))+1.D0+q0*z_sdss(i))
dist_c_sdss(i)=light/H0*z_sdss(i)*(1.D0+((z_sdss(i)*(1.D0-q0))/divisor))
log_lum_dist_sdss(i)=LOG10(dist_c_sdss(i)*1.D6)
dist_c_sdss(i)=dist_c_sdss(i)/(1.D0+z_sdss(i))*1.D6
mass_sdss(i)=(10.D0**mass_sdss(i))
END DO
 
 
 DO i=0,50
 bin_mass_2mrs(i)=0.D0
 bin_mass_sdss(i)=0.D0
 END DO
 
 DO i=1,n_2mrs
help_i=NINT((dist_c_2mrs(i)-5.D6)/10.D6)
IF (help_i<51) THEN
bin_mass_2mrs(help_i)=bin_mass_2mrs(help_i)+mass_2mrs(i)
END IF
END DO
 
 DO i=1,n_sdss
help_i=NINT((dist_c_sdss(i)-5.D6)/10.D6)
IF (help_i<51) THEN
bin_mass_sdss(help_i)=bin_mass_sdss(help_i)+mass_sdss(i)
END IF
END DO 
 
DO i=0,50
V_c(i)=4.D0/3.D0*PI*(((DBLE(i+1)*10.D6)**3)-((DBLE(i)*10.D6)**3))
bin_density_2mrs(i)=bin_mass_2mrs(i)/(V_c(i)*twomrs_cover)
bin_density_sdss(i)=bin_mass_sdss(i)/(V_c(i)*sdss_cover)
bin_density_2mrs(i)=bin_density_2mrs(i)/(rho_crit*Omega_m)
bin_density_sdss(i)=bin_density_sdss(i)/(rho_crit*Omega_m)
END DO

 

 OPEN(50,file='bin_dens_2mrs.txt')
DO i=0,50
dummy_var=(DBLE(i)*10.D6)+5.D6
WRITE(50,*) (dummy_var/1.D6),bin_density_2mrs(i)
END DO
 CLOSE(50)
 
 OPEN(50,file='bin_dens_sdss.txt')
DO i=0,50
dummy_var=(DBLE(i)*10.D6)+5.D6
WRITE(50,*) (dummy_var/1.D6),bin_density_sdss(i)
END DO
 CLOSE(50)
 
dummy_var=0.D0
 DO i=0,50
dummy_var=dummy_var+bin_mass_sdss(i)
END DO
! WRITE(*,*) dummy_var/sdss_cover
dummy_var=0.D0
 DO i=0,50
dummy_var=dummy_var+bin_mass_2mrs(i)
END DO
!  WRITE(*,*) dummy_var/twomrs_cover
!  
  WRITE(*,*) 'densities binned'
  
  DO ii=1,6
 DO i=1,n_mfi_first(ii)
m_halo_first(i,ii)=LOG10(m_halo_first(i,ii))
m_fi_first(i,ii)=LOG10(m_fi_first(i,ii))
END DO 
  
  DO i=1,n_mfi_final(ii)
m_halo_final(i,ii)=LOG10(m_halo_final(i,ii))
m_fi_final(i,ii)=LOG10(m_fi_final(i,ii))
END DO
END DO
 
 
DO iiii=1,6
 
 IF (iiii==1) THEN
 numbersnap='1'
 END IF
 IF (iiii==2) THEN
 numbersnap='2'
 END IF
 IF (iiii==3) THEN
 numbersnap='3'
 END IF
 IF (iiii==4) THEN
 numbersnap='4'
 END IF
 IF (iiii==5) THEN
 numbersnap='5'
 END IF
 IF (iiii==6) THEN
 numbersnap='6'
 END IF
 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_mfi_first(iiii)
help_x=NINT((m_halo_first(i,iiii)-10.D0)*100.D0)
help_y=NINT((m_fi_first(i,iiii)-10.D0)*100.D0)
IF ((help_x>0).AND.(help_x<500)) THEN
IF ((help_y>0).AND.(help_y<500)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='map_fi_first_'//TRIM(numbersnap)//'.txt')
DO i=0,500
DO ii=0,500
help_x2=(DBLE(i)/100.D0)+10.D0
help_y2=(DBLE(ii)/100.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_mfi_final(iiii)
help_x=NINT((m_halo_final(i,iiii)-10.D0)*100.D0)
help_y=NINT((m_fi_final(i,iiii)-10.D0)*100.D0)
IF ((help_x>0).AND.(help_x<500)) THEN
IF ((help_y>0).AND.(help_y<500)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='map_fi_final_'//TRIM(numbersnap)//'.txt')
DO i=0,500
DO ii=0,500
help_x2=(DBLE(i)/100.D0)+10.D0
help_y2=(DBLE(ii)/100.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 

DO i=1,n_mfi_first(iiii)
datalist_first(1,i)=m_fi_first(i,iiii)
datalist_first(2,i)=m_halo_first(i,iiii)**2
datalist_first(3,i)=m_halo_first(i,iiii)
datalist_first(4,i)=1.D0
END DO
 
  DO i=1,3
 Yvector(i)=0.D0
 DO ii=1,3
 Amatrix(i,ii)=0.D0
 DO iii=1,3
 Bmatrices(i,ii,iii)=0.D0 
 END DO
 END DO
 END DO
 
 DO i=1,n_mfi_first(iiii)
 
 DO ii=1,3
 DO iii=1,3
 Amatrix(ii,iii)=Amatrix(ii,iii)+datalist_first(ii+1,i)*datalist_first(iii+1,i)
 END DO
 END DO
 
 DO ii=1,3
 Yvector(ii)=Yvector(ii)+datalist_first(1,i)*datalist_first(ii+1,i)
 END DO
 
 END DO
 
 det_A=det3d(Amatrix)       !                    determinate

DO iii=1,3
DO i=1,3
DO ii=1,3
Bmatrices(i,ii,iii)=Amatrix(i,ii)
END DO
END DO
END DO

DO iii=1,3
DO i=1,3
Bmatrices(iii,i,iii)=Yvector(i)
END DO
END DO

DO iii=1,3
DO i=1,3
DO ii=1,3
help_matrix(i,ii)=Bmatrices(i,ii,iii)
END DO
END DO
det_B(iii)=det3d(help_matrix)                           ! determinate
END DO


DO iii=1,3
solution_first(iii,iiii)=det_B(iii)/det_A
END DO
! get solutions


! propagation of error
!root mean square
rms=0.D0


DO i=1,n_mfi_first(iiii)

otherside=0.D0
DO ii=1,3
otherside=otherside+solution_first(ii,iiii)*datalist_first(ii+1,i)
END DO 

deviation_first(i)=otherside-datalist_first(1,i)
rms=rms+(deviation_first(i)**2)

END DO
rms_first(iiii)=SQRT(rms/DBLE(n_mfi_first(iiii)))


! error of fitting parameters

DO i=1,3

hc1=0
DO ii=1,3
IF (ii.NE.i) THEN
hc1=hc1+1
hc2=0
DO iii=1,3
IF (iii.NE.i) THEN
hc2=hc2+1
help_submatrix(hc1,hc2)=Amatrix(ii,iii)
END IF
END DO
END IF
END DO


errorbar_first(i,iiii)=det2d(help_submatrix)/det_A            !subdeterminate

errorbar_first(i,iiii)=SQRT(errorbar_first(i,iiii))*rms_first(iiii)
END DO

 
 
 
 
 
 
 
 
 

DO i=1,n_mfi_final(iiii)
datalist_final(1,i)=m_fi_final(i,iiii)
datalist_final(2,i)=m_halo_final(i,iiii)**2
datalist_final(3,i)=m_halo_final(i,iiii)
datalist_final(4,i)=1.D0
END DO
 
  DO i=1,3
 Yvector(i)=0.D0
 DO ii=1,3
 Amatrix(i,ii)=0.D0
 DO iii=1,3
 Bmatrices(i,ii,iii)=0.D0 
 END DO
 END DO
 END DO
 
 DO i=1,n_mfi_final(iiii)
 
 DO ii=1,3
 DO iii=1,3
 Amatrix(ii,iii)=Amatrix(ii,iii)+datalist_final(ii+1,i)*datalist_final(iii+1,i)
 END DO
 END DO
 
 DO ii=1,3
 Yvector(ii)=Yvector(ii)+datalist_final(1,i)*datalist_final(ii+1,i)
 END DO
 
 END DO
 
 det_A=det3d(Amatrix)       !                    determinate

DO iii=1,3
DO i=1,3
DO ii=1,3
Bmatrices(i,ii,iii)=Amatrix(i,ii)
END DO
END DO
END DO

DO iii=1,3
DO i=1,3
Bmatrices(iii,i,iii)=Yvector(i)
END DO
END DO

DO iii=1,3
DO i=1,3
DO ii=1,3
help_matrix(i,ii)=Bmatrices(i,ii,iii)
END DO
END DO
det_B(iii)=det3d(help_matrix)                           ! determinate
END DO


DO iii=1,3
solution_final(iii,iiii)=det_B(iii)/det_A
END DO
! get solutions


! propagation of error
!root mean square
rms=0.D0


DO i=1,n_mfi_final(iiii)

otherside=0.D0
DO ii=1,3
otherside=otherside+solution_final(ii,iiii)*datalist_final(ii+1,i)
END DO 

deviation_final(i)=otherside-datalist_final(1,i)
rms=rms+(deviation_final(i)**2)

END DO
rms_final(iiii)=SQRT(rms/DBLE(n_mfi_final(iiii)))


! error of fitting parameters

DO i=1,3

hc1=0
DO ii=1,3
IF (ii.NE.i) THEN
hc1=hc1+1
hc2=0
DO iii=1,3
IF (iii.NE.i) THEN
hc2=hc2+1
help_submatrix(hc1,hc2)=Amatrix(ii,iii)
END IF
END DO
END IF
END DO


errorbar_final(i,iiii)=det2d(help_submatrix)/det_A            !subdeterminate

errorbar_final(i,iiii)=SQRT(errorbar_final(i,iiii))*rms_final(iiii)
END DO

 
 
  OPEN(61,file='fi_fit_para_all_'//TRIM(numbersnap)//'.txt')
WRITE(61,*) solution_first(1:3,iiii)
WRITE(61,*) errorbar_first(1:3,iiii)
WRITE(61,*) rms_first(iiii)
WRITE(61,*)
WRITE(61,*) solution_final(1:3,iiii)
WRITE(61,*) errorbar_final(1:3,iiii)
WRITE(61,*) rms_final(iiii)
 CLOSE(61)
 
 
 

 
 OPEN(61,file='fi_fit_para_'//TRIM(numbersnap)//'.txt')
WRITE(61,*) 'a_first_1=',solution_first(1,iiii)
WRITE(61,*) 'a_first_2=',solution_first(2,iiii)
WRITE(61,*) 'a_first_3=',solution_first(3,iiii)
WRITE(61,*) 'a_final_1=',solution_final(1,iiii)
WRITE(61,*) 'a_final_2=',solution_final(2,iiii)
WRITE(61,*) 'a_final_3=',solution_final(3,iiii)
 CLOSE(61)
 
 END DO
 
 
 WRITE(*,*) 'mass rescale functions fitted'
 
 
DO i=1,n_all
active(i)=.TRUE.
END DO
 
DO i=1,n_2mrs
ra(i)=ra_2mrs(i)
dec(i)=dec_2mrs(i)
red(i)=z_2mrs(i)
mass(i)=mass_2mrs(i)
mass_err(i)=mass_2mrs_err(i)
dist_c(i)=dist_c_2mrs(i)
members(i)=members_2mrs(i)
angrad(i)=angrad_2mrs(i)
sigma(i)=sigma_2mrs(i)
END DO

DO i=1,n_sdss
ra(i+n_2mrs)=ra_sdss(i)
dec(i+n_2mrs)=dec_sdss(i)
red(i+n_2mrs)=z_sdss(i)
mass(i+n_2mrs)=mass_sdss(i)
mass_err(i+n_2mrs)=mass_sdss_err(i)
dist_c(i+n_2mrs)=dist_c_sdss(i)
members(i+n_2mrs)=members_sdss(i)
angrad(i+n_2mrs)=angrad_sdss(i)
sigma(i+n_2mrs)=sigma_sdss(i)
END DO

 
DO i=1,n_all
ra(i)=ra(i)*PI/180.D0
dec(i)=dec(i)*PI/180.D0
angrad(i)=angrad(i)*PI/180.D0
sigma(i)=sigma(i)/light
dist_c(i)=dist_c(i)/1.D6
px(i)=dist_c(i)*COS(dec(i))*COS(ra(i))
py(i)=dist_c(i)*COS(dec(i))*SIN(ra(i))
pz(i)=dist_c(i)*SIN(dec(i))
END DO

 WRITE(*,*) 'parameters for combined catalogue calculated'
! 



IF (limiting_funct==1) THEN
! get length of file
OPEN(50,file='ts_sor/SDSS_real/SDSS_obs_galaxies_radecz.txt')
io_err=0
n_gal_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_gal_sdss=n_gal_sdss+1
END DO
 CLOSE(50)
n_gal_sdss=n_gal_sdss-1


allocate(gal_z_sdss(1:n_gal_sdss))
allocate(gal_mag_r(1:n_gal_sdss))
allocate(gal_active_sdss(1:n_gal_sdss))
allocate(volume_weight_sdss(1:n_gal_sdss))



! read file
OPEN(50,file='ts_sor/SDSS_real/SDSS_obs_galaxies_radecz.txt')
DO i=1,n_gal_sdss
READ(50,*) dummy_var,dummy_var,gal_z_sdss(i),dummy_var,gal_mag_r(i)
END DO
CLOSE(50)

 WRITE(*,*) n_gal_sdss,'SDSS galaxies read in'
! 
! 
mag_limit_sdss=17.77D0
mag_saturation_sdss=14.D0
! calculate volume weights to corret for Malmquist bias
V_sum=0.D0
DO i=1,n_gal_sdss

gal_active_sdss(i)=.TRUE.
volume_weight_sdss(i)=0.D0

IF (gal_mag_r(i)<-30.D0) THEN
gal_active_sdss(i)=.FALSE.
END IF

IF (gal_mag_r(i)>-15.D0) THEN
gal_active_sdss(i)=.FALSE.
END IF


IF (gal_active_sdss(i)) THEN
D_limit=-0.2D0*gal_mag_r(i)+((mag_limit_sdss+5.D0)/5.D0)
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


V_limit=V_limit*sdss_cover


volume_weight_sdss(i)=1.D0/V_limit

V_sum=V_sum+volume_weight_sdss(i)

END IF

END DO 

WRITE(*,*) 'SDSS weights calculated'
 

! rescale linking length according to the completeness of the luminosity function 
DO i=1,n_sdss
vol_sat_sdss(i)=1.D-15 !about zero, but to avoid division by zero later we make it a little larger
mabs_volsat=mag_saturation_sdss+5.D0-5.D0*log_lum_dist_sdss(i)
DO ii=1,n_gal_sdss
IF (gal_active_sdss(ii)) THEN
IF (gal_mag_r(ii)>mabs_volsat) THEN
vol_sat_sdss(i)=vol_sat_sdss(i)+volume_weight_sdss(ii)
END IF
END IF
END DO
END DO
WRITE(*,*) 'Saturation calculated'




! rescale linking length according to the completeness of the luminosity function 
DO i=1,n_sdss
vol_int_sdss(i)=1.D-15 !about zero, but to avoid division by zero later we make it a little larger
mabs_vollim=mag_limit_sdss+5.D0-5.D0*log_lum_dist_sdss(i)
mabs_volsat=mag_saturation_sdss+5.D0-5.D0*log_lum_dist_sdss(i)
DO ii=1,n_gal_sdss
IF (gal_active_sdss(ii)) THEN
IF (gal_mag_r(ii)<mabs_vollim) THEN
IF (gal_mag_r(ii)>mabs_volsat) THEN
vol_int_sdss(i)=vol_int_sdss(i)+volume_weight_sdss(ii)
END IF
END IF
END IF
END DO

END DO

WRITE(*,*) 'Saturation+Malmquist calculated'

! rescale linking length according to the completeness of the luminosity function 
DO i=1,n_sdss
v_cor_sdss(i)=1.D-15 !about zero, but to avoid division by zero later we make it a little larger
mabs_vollim=mag_limit_sdss+5.D0-5.D0*log_lum_dist_sdss(i)
DO ii=1,n_gal_sdss
IF (gal_active_sdss(ii)) THEN
IF (gal_mag_r(ii)<mabs_vollim) THEN
v_cor_sdss(i)=v_cor_sdss(i)+volume_weight_sdss(ii)
END IF
END IF
END DO

END DO
WRITE(*,*) 'Malmquist calculated'

DO i=1,n_sdss
vol_sat_sdss(i)=vol_sat_sdss(i)/V_sum
vol_int_sdss(i)=vol_int_sdss(i)/V_sum
v_cor_sdss(i)=v_cor_sdss(i)/V_sum
END DO



 D_limit=1000.D0*1.D6
DO i=1,n_sdss
IF (vol_sat_sdss(i)>0.95D0) THEN
IF (dist_c_sdss(i)<D_limit) THEN
D_limit=dist_c_sdss(i)
END IF
END IF
END DO
WRITE(*,*) D_limit










 OPEN(50,file='vol_cor_abs.txt')
DO i=1,n_sdss
WRITE(50,*) (dist_c_sdss(i)/1.D6),(vol_sat_sdss(i)*V_sum),(v_cor_sdss(i)*V_sum),(vol_int_sdss(i)*V_sum)
END DO
 CLOSE(50)
 

 OPEN(50,file='vol_cor.txt')
DO i=1,n_sdss
WRITE(50,*) (dist_c_sdss(i)/1.D6),vol_sat_sdss(i),v_cor_sdss(i),vol_int_sdss(i)
END DO
 CLOSE(50)
 
  OPEN(50,file='vol_cor_plot.txt')
DO i=1,n_sdss
IF ((dist_c_sdss(i)/1.D6)<50.D0) THEN
WRITE(50,*) (dist_c_sdss(i)/1.D6),vol_sat_sdss(i),v_cor_sdss(i),vol_int_sdss(i)
ELSE
IF (rand()>0.99D0) THEN
WRITE(50,*) (dist_c_sdss(i)/1.D6),vol_sat_sdss(i),v_cor_sdss(i),vol_int_sdss(i)
END IF
END IF
END DO
 CLOSE(50)
 
 

WRITE(*,*) 'SDSS completeness corrections calculated'
 
 
 
 
 mag_limit_2mrs=11.75D0
 
 
 
 
! get length of file
OPEN(50,file='ts_sor/2MRS_real/2MRS_obs_galaxies_radecz.txt')
io_err=0
n_gal_2mrs=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_gal_2mrs=n_gal_2mrs+1
END DO
 CLOSE(50)
n_gal_2mrs=n_gal_2mrs-1


allocate(gal_z_2mrs(1:n_gal_2mrs))
allocate(gal_mag_K(1:n_gal_2mrs))
allocate(gal_active_2mrs(1:n_gal_2mrs))
allocate(volume_weight_2mrs(1:n_gal_2mrs))



! read file
OPEN(50,file='ts_sor/2MRS_real/2MRS_obs_galaxies_radecz.txt')
DO i=1,n_gal_2mrs
READ(50,*) dummy_var,dummy_var,gal_z_2mrs(i),dummy_var,gal_mag_K(i)
END DO
CLOSE(50)

 WRITE(*,*) n_gal_2mrs,'2MRS galaxies read in'
 
 

 
! calculate volume weights to corret for Malmquist bias
V_sum=0.D0
DO i=1,n_gal_2mrs

gal_active_2mrs(i)=.TRUE.
volume_weight_2mrs(i)=0.D0

IF (gal_mag_K(i)<-30.D0) THEN
gal_active_2mrs(i)=.FALSE.
END IF

IF (gal_mag_K(i)>-18.D0) THEN
gal_active_2mrs(i)=.FALSE.
END IF

IF (gal_z_2mrs(i)>0.11D0) THEN
gal_active_2mrs(i)=.FALSE.
END IF


IF (gal_active_2mrs(i)) THEN
D_limit=-0.2D0*gal_mag_K(i)+((mag_limit_2mrs+5.D0)/5.D0)
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


V_limit=V_limit*twomrs_cover


volume_weight_2mrs(i)=1.D0/V_limit

V_sum=V_sum+volume_weight_2mrs(i)

END IF

END DO 


 
 WRITE(*,*) '2MRS weights calculated'
 
! rescale linking length according to the completeness of the luminosity function 
DO i=1,n_2mrs
v_cor_2mrs(i)=1.D-15 !about zero, but to avoid division by zero later we make it a little larger
mabs_vollim=mag_limit_2mrs+5.D0-5.D0*log_lum_dist_2mrs(i)
DO ii=1,n_gal_2mrs
IF (gal_active_2mrs(ii)) THEN
IF (gal_mag_K(ii)<mabs_vollim) THEN
v_cor_2mrs(i)=v_cor_2mrs(i)+volume_weight_2mrs(ii)
END IF
END IF
END DO

END DO


DO i=1,n_2mrs
v_cor_2mrs(i)=v_cor_2mrs(i)/V_sum
END DO

 
 WRITE(*,*) 'Malmquist calculated'
 


 OPEN(50,file='vol_cor2_abs.txt')
DO i=1,n_2mrs
WRITE(50,*) (dist_c_2mrs(i)/1.D6),(v_cor_2mrs(i)*V_sum)
END DO
 CLOSE(50)
 


 OPEN(50,file='vol_cor2.txt')
DO i=1,n_2mrs
WRITE(50,*) (dist_c_2mrs(i)/1.D6),v_cor_2mrs(i)
END DO
 CLOSE(50)
 
  OPEN(50,file='vol_cor2_plot.txt')
DO i=1,n_2mrs
IF ((dist_c_2mrs(i)/1.D6)<50.D0) THEN
WRITE(50,*) (dist_c_2mrs(i)/1.D6),v_cor_2mrs(i)
ELSE
IF (rand()>0.99D0) THEN
WRITE(50,*) (dist_c_2mrs(i)/1.D6),v_cor_2mrs(i)
END IF
END IF
END DO
 CLOSE(50)


WRITE(*,*) '2MRS completeness corrections calculated'
 
 
 ELSE
 
  OPEN(50,file='vol_cor.txt')
DO i=1,n_sdss
READ(50,*) dummy_var,vol_sat_sdss(i),v_cor_sdss(i),vol_int_sdss(i)
END DO
 CLOSE(50)
 
   OPEN(50,file='vol_cor2.txt')
DO i=1,n_2mrs
READ(50,*) dummy_var,v_cor_2mrs(i)
END DO
 CLOSE(50)
 

 
END IF

 
 D_limit=1000.D0*1.D6
DO i=1,n_sdss
IF (vol_sat_sdss(i)>0.95D0) THEN
IF (dist_c_sdss(i)<D_limit) THEN
D_limit=dist_c_sdss(i)
END IF
END IF
END DO



d_limit_2mrs=D_limit/1.D6
WRITE(*,*) d_limit_2mrs
iii=0
DO i=1,n_2mrs
IF (dist_c(i)>d_limit_2mrs) THEN
active(i)=.FALSE.
iii=iii+1
END IF
END DO 
WRITE(*,*) iii,'2MRS clusters removed for being to far away'




DO i=1,n_2mrs
IF (active(i)) THEN
D_limit=dist_c(i)*((1.D0+z_2mrs(i))**(-1.D0))
basic_link_a(i)=ATAN(basic_link/D_limit)!*180.D0/PI
IF (basic_link_a(i)<0.D0) THEN
basic_link_a(i)=-basic_link_a(i)
END IF
basic_link_R(i)=(basic_link*H0+av_rad_pec_vel*2.D0)/light
END IF
END DO
WRITE(*,*) 'redshift distortions and projection effects considered'




! modify linking lenght according to the modification calculated before
DO i=1,n_2mrs
IF (active(i)) THEN
b_eff(i)=b_factor*basic_link_a(i)*((v_cor_2mrs(i))**(-lum_power/3.D0))
R_eff(i)=R_factor*basic_link_R(i)*((v_cor_2mrs(i))**(-lum_power/3.D0))
END IF
END DO
! adjust linking length if smaller than group size
DO i=1,n_2mrs
IF (active(i)) THEN
IF (b_eff(i)<angrad(i)) THEN
b_eff(i)=angrad(i)
END IF
IF (R_eff(i)<sigma(i)) THEN
R_eff(i)=sigma(i)
END IF
END IF
END DO

WRITE(*,*) 'linking lengths for crossmatch calibrated'

help_i=0
iii=0
iiii=0
DO ii=1,n_sdss
sdss_in_2mrs(ii)=.FALSE.
IF (dist_c(ii+n_2mrs)<(1.1D0*d_limit_2mrs)) THEN
iiii=iiii+1
END IF
END DO
WRITE(*,*) iiii,'SDSS groups in the overlapping area'

  OPEN(50,file='merged_masses.txt')

WRITE(50,*)



DO i=1,n_2mrs
IF (active(i)) THEN

iiii=0
DO ii=1,n_sdss

sdss_in_2mrs(ii)=.FALSE.
IF (active(ii+n_2mrs)) THEN
IF (dist_c(ii+n_2mrs)<(1.1D0*d_limit_2mrs)) THEN

delta_z=z_2mrs(i)-z_sdss(ii)
IF (delta_z<0.D0) THEN
delta_z=-delta_z
END IF
IF (delta_z<R_eff(i)) THEN


angular_sep=COS(dec(i))*COS(dec(ii+n_2mrs))*COS((ra(ii+n_2mrs)-ra(i)))
angular_sep=angular_sep+SIN(dec(i))*SIN(dec(ii+n_2mrs))
IF (angular_sep>1.D0) THEN
angular_sep=1.D0
END IF
IF (angular_sep<-1.D0) THEN
angular_sep=-1.D0
END IF
angular_sep=ACOS(angular_sep)

IF (angular_sep<0.D0) THEN
angular_sep=-angular_sep
END IF

! help_mass_err,help_mass_err_weight
IF (angular_sep<b_eff(i)) THEN

WRITE(50,*) i,ii,dist_c(i),LOG10(mass(i)),LOG10(mass(ii+n_2mrs)),&
(mass(i)/mass(ii+n_2mrs)),members(i),members(ii+n_2mrs)
sdss_in_2mrs(ii)=.TRUE.
iii=iii+1
iiii=iiii+1
END IF

END IF
END IF
END IF
END DO

IF (iiii>0) THEN
help_i=help_i+1
members_help=members(i)
weights_help=DBLE(members(i))*v_cor_2mrs(i)
px_help=px(i)*mass(i)*weights_help
py_help=py(i)*mass(i)*weights_help
pz_help=pz(i)*mass(i)*weights_help
mass_help=mass(i)*weights_help
weight_sum=weights_help
help_mass_err=weights_help*(mass_err(i)**2)
! ! help_mass_err_weight=

DO ii=1,n_sdss
IF (sdss_in_2mrs(ii)) THEN

members_help=members_help+members(ii+n_2mrs)
weights_help=DBLE(members(ii+n_2mrs))*vol_int_sdss(ii)
px_help=px_help+px(ii+n_2mrs)*mass(ii+n_2mrs)*weights_help
py_help=py_help+py(ii+n_2mrs)*mass(ii+n_2mrs)*weights_help
pz_help=pz_help+pz(ii+n_2mrs)*mass(ii+n_2mrs)*weights_help
mass_help=mass_help+mass(ii+n_2mrs)*weights_help
weight_sum=weight_sum+weights_help
help_mass_err=help_mass_err+weights_help*(mass_err(ii+n_2mrs)**2)

active(ii+n_2mrs)=.FALSE.
END IF
END DO

mass_err(i)=SQRT(help_mass_err/weight_sum)

px(i)=px_help/mass_help
py(i)=py_help/mass_help
pz(i)=pz_help/mass_help
mass(i)=mass_help/weight_sum
dist_c(i)=SQRT((px(i)**2)+(py(i)**2)+(pz(i)**2))

END IF


END IF
END DO

 CLOSE(50)


WRITE(*,*) iii,'SDSS groups in overlap merged with',help_i,'2MRS groups'


iiii=0
DO i=1,n_all
IF (active(i)) THEN
iiii=iiii+1
END IF
END DO

WRITE(*,*) iiii,'groups in basic catalogue'





 DO i=0,50
 bin_mass_combi(i)=0.D0
 END DO
 
 DO i=1,n_all
IF (active(i)) THEN
help_i=NINT((dist_c(i)-5.D0)/10.D0)
IF (help_i<51) THEN
bin_mass_combi(help_i)=bin_mass_combi(help_i)+mass(i)
END IF
END IF
END DO
 


DO i=0,50
V_c(i)=4.D0/3.D0*PI*(((DBLE(i+1)*10.D6)**3)-((DBLE(i)*10.D6)**3))

IF ((DBLE(i+1)*10.D0)<d_limit_2mrs) THEN
bin_density_combi(i)=bin_mass_combi(i)/(V_c(i)*twomrs_cover)
ELSE
IF ((DBLE(i)*10.D0)>d_limit_2mrs) THEN
bin_density_combi(i)=bin_mass_combi(i)/(V_c(i)*sdss_cover)
ELSE
V_limit=4.D0/3.D0*PI*(((d_limit_2mrs*1.D6)**3)-((DBLE(i)*10.D6)**3))*twomrs_cover
V_limit=V_limit+4.D0/3.D0*PI*(((DBLE(i+1)*10.D6)**3)-((d_limit_2mrs*1.D6)**3))*sdss_cover
bin_density_combi(i)=bin_mass_combi(i)/V_limit
END IF
END IF

bin_density_combi(i)=bin_density_combi(i)/(rho_crit*Omega_m)
END DO
! 
 

 OPEN(50,file='bin_dens_combi.txt')
DO i=0,50
dummy_var=(DBLE(i)*10.D6)+5.D6
WRITE(50,*) (dummy_var/1.D6),bin_density_combi(i)
END DO
 CLOSE(50)

 
 
 
 z_limit=0.11D0
divisor=(SQRT(1.D0+2.D0*q0*z_limit)+1.D0+q0*z_limit)
D_limit=light/H0*z_limit*(1.D0+((z_limit*(1.D0-q0))/divisor))
D_limit=D_limit*((1.D0+z_limit)**(-1.D0))
 d_maxdepth=D_limit
 OPEN(50,file='bin_plot_help.txt')

WRITE(50,*) 'survey_limit=',d_maxdepth
WRITE(50,*) 'merge_limit=',d_limit_2mrs
 CLOSE(50)

 
 
WRITE(*,*) 'combined densities binned'


 sum_mass=0.D0
  DO i=1,n_all
IF (active(i)) THEN
 sum_mass=sum_mass+mass(i)
 END IF
END DO

V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)
V_inner=4.D0/3.D0*PI*((d_limit_2mrs*1.D6)**3)*twomrs_cover
V_outer=4.D0/3.D0*PI*(((d_maxdepth*1.D6)**3)-((d_limit_2mrs*1.D6)**3))*sdss_cover
V_sum=V_inner+V_outer

! 
WRITE(*,*) sum_mass
WRITE(*,*) V_inner,V_outer,V_sum
final_dens=sum_mass/V_sum
final_dens=final_dens/(rho_crit*Omega_m)
WRITE(*,*) final_dens,'of the total matter density in catalogue after combination'


DO i=1,n_all
IF (active(i)) THEN
mass_groups(i)=mass(i)
rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6

END IF
END DO
 
 
DO i=1,n_all
IF (active(i)) THEN



ii=0
DO WHILE (ii<n_all)
ii=ii+1
IF (active(ii)) THEN


IF (i/=ii) THEN

d=(px(i)-px(ii))**2+(py(i)-py(ii))**2+(pz(i)-pz(ii))**2
d=SQRT(d)
IF ((d+rad(ii))<rad(i)) THEN


active(ii)=.FALSE.

mass_groups(i)=mass_groups(i)+mass_groups(ii)

mass_err(i)=SQRT(((mass_err(i)**2)+(mass_err(ii)**2))/2.D0)
mnew=mass(i)+mass(ii)
px(i)=(px(i)*mass(i)+px(ii)*mass(ii))/mnew
py(i)=(py(i)*mass(i)+py(ii)*mass(ii))/mnew
pz(i)=(pz(i)*mass(i)+pz(ii)*mass(ii))/mnew
mass(i)=mnew



rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6

ii=0
END IF
END IF
END IF

END DO


END IF
END DO

iii=0
DO i=1,n_all
IF (active(i)) THEN
iii=iii+1
END IF
END DO
n_new=iii
 WRITE(*,*) iii,'groups remain after first fi filtering'
 
 
DO i=1,n_all
IF (active(i)) THEN
! 

gal_dist_c=((px(i)**2)+(py(i)**2)+(pz(i)**2))
gal_dist_c=SQRT(gal_dist_c)
gal_redshift=(light**2)*((light**2)+(gal_dist_c**2)*(H0**2)*(1.D0-2.D0*q0))*(q0-1.D0)**2
gal_redshift=SQRT(gal_redshift)
gal_redshift=gal_redshift+light*gal_dist_c*H0+(light**2)*(q0-1.D0)-(gal_dist_c**2)*(H0**2)*(q0**2)
gal_redshift=gal_redshift/((light-gal_dist_c*H0*q0)**2)

evolution_seg=1
 DO ii=1,6
 zmax=(redshift_snap(ii)+redshift_snap(ii+1))/2.D0
 zmin=(redshift_snap(ii)+redshift_snap(ii-1))/2.D0

 IF ((gal_redshift<zmax).AND.(gal_redshift>zmin)) THEN
  evolution_seg=ii
END IF
 
 END DO
 
 

logmh=LOG10(mass(i))

logmfi=solution_first(1,evolution_seg)*(logmh**2)+solution_first(2,evolution_seg)*logmh+solution_first(3,evolution_seg)

! ! GAUSS FEHLER
! help_mass_err=((errorbar_first(1,evolution_seg)*(logmh**2))**2)+((errorbar_first(2,evolution_seg)*logmh)**2)
! help_mass_err=help_mass_err+(errorbar_first(3,evolution_seg)**2)+&
! (((2.D0*solution_first(1,evolution_seg)*logmh+solution_first(2,evolution_seg))*mass_err(i))**2)
! mass_err(i)=SQRT(help_mass_err)
mass_err(i)=rms_first(evolution_seg)

mass(i)=10.D0**logmfi

 rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6

END IF
END DO

 WRITE(*,*) 'first fi rescaling'
 
 
 
 
  sum_mass=0.D0
  DO i=1,n_all
IF (active(i)) THEN
 sum_mass=sum_mass+mass(i)
 END IF
END DO

V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)
V_inner=4.D0/3.D0*PI*((d_limit_2mrs*1.D6)**3)*twomrs_cover
V_outer=4.D0/3.D0*PI*(((d_maxdepth*1.D6)**3)-((d_limit_2mrs*1.D6)**3))*sdss_cover
V_sum=V_inner+V_outer

! 
! WRITE(*,*) sum_mass
! WRITE(*,*) V_inner,V_outer,V_sum
final_dens=sum_mass/V_sum
final_dens=final_dens/(rho_crit*Omega_m)
WRITE(*,*) final_dens,'of the total matter density after first rescaling'





 
 
 
 
 
  n_old=0
 DO WHILE (n_old.NE.n_new)
  WRITE(*,*) 'status:',n_old,n_new
 n_old=n_new
 
 
 DO i=1,n_all
IF (active(i)) THEN



ii=0
DO WHILE (ii<n_all)
ii=ii+1
IF (active(ii)) THEN


IF (i/=ii) THEN

d=(px(i)-px(ii))**2+(py(i)-py(ii))**2+(pz(i)-pz(ii))**2
d=SQRT(d)
IF ((d+rad(ii))<rad(i)) THEN


active(ii)=.FALSE.
mass_groups(i)=mass_groups(i)+mass_groups(ii)
mass_err(i)=SQRT(((mass_err(i)**2)+(mass_err(ii)**2))/2.D0)
mnew=mass(i)+mass(ii)
px(i)=(px(i)*mass(i)+px(ii)*mass(ii))/mnew
py(i)=(py(i)*mass(i)+py(ii)*mass(ii))/mnew
pz(i)=(pz(i)*mass(i)+pz(ii)*mass(ii))/mnew
mass(i)=mnew

rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6

ii=0
END IF
END IF
END IF

END DO


END IF
END DO

iii=0
DO i=1,n_all
IF (active(i)) THEN
iii=iii+1
END IF
END DO
 n_new=iii
 
 END DO
 
 WRITE(*,*) n_new,'groups remain after interative fi filtering'
 
 
 
 DO i=1,n_all
IF (active(i)) THEN


gal_dist_c=((px(i)**2)+(py(i)**2)+(pz(i)**2))
gal_dist_c=SQRT(gal_dist_c)
gal_redshift=(light**2)*((light**2)+(gal_dist_c**2)*(H0**2)*(1.D0-2.D0*q0))*(q0-1.D0)**2
gal_redshift=SQRT(gal_redshift)
gal_redshift=gal_redshift+light*gal_dist_c*H0+(light**2)*(q0-1.D0)-(gal_dist_c**2)*(H0**2)*(q0**2)
gal_redshift=gal_redshift/((light-gal_dist_c*H0*q0)**2)

evolution_seg=1
 DO ii=1,6
 zmax=(redshift_snap(ii)+redshift_snap(ii+1))/2.D0
 zmin=(redshift_snap(ii)+redshift_snap(ii-1))/2.D0

 IF ((gal_redshift<zmax).AND.(gal_redshift>zmin)) THEN
  evolution_seg=ii
END IF
 
 END DO
 


logmh=LOG10(mass_groups(i))


logmfi=solution_final(1,evolution_seg)*(logmh**2)+solution_final(2,evolution_seg)*logmh+solution_final(3,evolution_seg)

help_mass_err=((errorbar_final(1,evolution_seg)*(logmh**2))**2)+((errorbar_final(2,evolution_seg)*logmh)**2)
help_mass_err=help_mass_err+(errorbar_final(3,evolution_seg)**2)+&
(((2.D0*solution_final(1,evolution_seg)*logmh+solution_final(2,evolution_seg))*mass_err(i))**2)
mass_err(i)=SQRT(help_mass_err)

mass_err(i)=SQRT((rms_first(evolution_seg)**2)+(rms_first(evolution_seg)**2))

mass(i)=10.D0**logmfi

 rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6

 rad_err(i)=((10.D0**mass_err(i))-1.D0)*rad(i)/SQRT(3.D0)
 
 
END IF
END DO

 sum_mass=0.D0
  DO i=1,n_all
IF (active(i)) THEN
 sum_mass=sum_mass+mass(i)
 END IF
END DO

V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)
V_inner=4.D0/3.D0*PI*((d_limit_2mrs*1.D6)**3)*twomrs_cover
V_outer=4.D0/3.D0*PI*(((d_maxdepth*1.D6)**3)-((d_limit_2mrs*1.D6)**3))*sdss_cover
V_sum=V_inner+V_outer

! 
! WRITE(*,*) sum_mass
! WRITE(*,*) V_inner,V_outer,V_sum
final_dens=sum_mass/V_sum
final_dens=final_dens/(rho_crit*Omega_m)
WRITE(*,*) final_dens,'of the total matter density in catalogue'














 DO i=0,50
 bin_mass_fi(i)=0.D0
 END DO
 
 DO i=1,n_all
IF (active(i)) THEN
dist_c(i)=SQRT((px(i)**2)+(py(i)**2)+(pz(i)**2))
help_i=NINT((dist_c(i)-5.D0)/10.D0)
IF (help_i<51) THEN
bin_mass_fi(help_i)=bin_mass_fi(help_i)+mass(i)
END IF
END IF
END DO
 


DO i=0,50
V_c(i)=4.D0/3.D0*PI*(((DBLE(i+1)*10.D6)**3)-((DBLE(i)*10.D6)**3))

IF ((DBLE(i+1)*10.D0)<d_limit_2mrs) THEN
bin_density_fi(i)=bin_mass_fi(i)/(V_c(i)*twomrs_cover)
ELSE
IF ((DBLE(i)*10.D0)>d_limit_2mrs) THEN
bin_density_fi(i)=bin_mass_fi(i)/(V_c(i)*sdss_cover)
ELSE
V_limit=4.D0/3.D0*PI*(((d_limit_2mrs*1.D6)**3)-((DBLE(i)*10.D6)**3))*twomrs_cover
V_limit=V_limit+4.D0/3.D0*PI*(((DBLE(i+1)*10.D6)**3)-((d_limit_2mrs*1.D6)**3))*sdss_cover
bin_density_fi(i)=bin_mass_fi(i)/V_limit
END IF
END IF

bin_density_fi(i)=bin_density_fi(i)/(rho_crit*Omega_m)
END DO
! 
 

 OPEN(50,file='bin_dens_fi.txt')
DO i=0,50
dummy_var=(DBLE(i)*10.D6)+5.D6
WRITE(50,*) (dummy_var/1.D6),bin_density_fi(i)
END DO
 CLOSE(50)

 
 







n_rand=100000


allocate(rand_x(1:n_rand))
allocate(rand_y(1:n_rand))
allocate(rand_z(1:n_rand))
allocate(rand_act(1:n_rand))
allocate(rand_infi(1:n_rand))
allocate(rand_infi_multi(1:n_rand))


DO i=1,n_rand
rand_x(i)=(RAND()*d_maxdepth*2.D0)-d_maxdepth
rand_y(i)=(RAND()*d_maxdepth*2.D0)-d_maxdepth
rand_z(i)=(RAND()*d_maxdepth*2.D0)-d_maxdepth
rand_act(i)=.TRUE.
rand_infi(i)=.FALSE.
rand_infi_multi(i)=.FALSE.
END DO

n_rand_used=n_rand
DO i=1,n_rand
d=(rand_x(i)**2)+(rand_y(i)**2)+(rand_z(i)**2)
IF (d>(d_maxdepth**2)) THEN
rand_act(i)=.FALSE.
n_rand_used=n_rand_used-1
END IF
END DO
WRITE(*,*) n_rand_used

DO i=1,n_all
! WRITE(*,*) i
IF (active(i)) THEN
DO ii=1,n_rand
IF (rand_act(ii)) THEN
d=((rand_x(ii)-px(i))**2)+((rand_y(ii)-py(i))**2)+((rand_z(ii)-pz(i))**2)
IF (d<(rad(i)**2)) THEN

IF (rand_infi(ii)) THEN
rand_infi_multi(i)=.TRUE.
rand_act(ii)=.FALSE.
END IF

rand_infi(ii)=.TRUE.
END IF
END IF
END DO
END IF
END DO

n_rand_fi=0
DO i=1,n_rand
IF (rand_infi(i)) THEN
n_rand_fi=n_rand_fi+1
END IF
END DO


n_rand_multi=0
DO i=1,n_rand
IF (rand_infi_multi(i)) THEN
n_rand_multi=n_rand_multi+1
END IF
END DO

V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)
Vol_fi=(DBLE(n_rand_multi)/V_sum)/(DBLE(n_rand_used)/V_limit)

WRITE(*,*) Vol_fi,'of total survey volume covered by overlapping fi regions'


V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)
WRITE(*,*) n_rand_fi
WRITE(*,*) V_sum,V_limit,(V_sum/V_limit)

Vol_fi=(DBLE(n_rand_fi)/V_sum)/(DBLE(n_rand_used)/V_limit)

WRITE(*,*) Vol_fi,'of total survey volume occupied by fi regions'



!  

DO i=1,n_all
IF (active(i)) THEN
 cdist_fi(i)=SQRT((px(i)**2)+(py(i)**2)+(pz(i)**2))
 IF (cdist_fi(i)>0.D0) THEN
 dec_fi(i)=ASIN(pz(i)/cdist_fi(i))
 ELSE
 dec_fi(i)=0.D0
 END IF
 ra_fi(i)=ATAN2(py(i),px(i))
 
 dec_fi(i)=dec_fi(i)*180.D0/PI
 ra_fi(i)=ra_fi(i)*180.D0/PI 
END IF
END DO
 
 
 

OPEN(50,file='ts_sor/catalogues/fi_info.txt')
WRITE(50,*) final_dens
WRITE(50,*) Vol_fi
 CLOSE(50)

OPEN(50,file='ts_sor/catalogues/fi_regions.txt')
 DO i=1,n_all
IF (active(i)) THEN
WRITE(50,*) ra_fi(i),dec_fi(i),cdist_fi(i),px(i),py(i),pz(i),LOG10(mass(i)),mass_err(i),rad(i),rad_err(i)
END IF
END DO
 CLOSE(50)
 
 
 OPEN(50,file='final_catalogues/fi_regions_basic.txt')
 DO i=1,n_all
IF (active(i)) THEN
WRITE(50,"(2F14.8,4F10.3,4F10.5)") ra_fi(i),dec_fi(i),cdist_fi(i),px(i),py(i),pz(i),&
LOG10(mass(i)),mass_err(i),rad(i),rad_err(i)
END IF
END DO
 CLOSE(50)
 
 
DO i=1,n_all
IF (active(i)) THEN
mass(i)=mass(i)/final_dens
mass_err(i)=mass_err(i) !relative Fehler bleibt gleich
 rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6
END IF
END DO




DO i=1,n_all
IF (active(i)) THEN



ii=0
DO WHILE (ii<n_all)
ii=ii+1
IF (active(ii)) THEN


IF (i/=ii) THEN

d=(px(i)-px(ii))**2+(py(i)-py(ii))**2+(pz(i)-pz(ii))**2
d=SQRT(d)
IF ((d+rad(ii))<rad(i)) THEN


active(ii)=.FALSE.

! mass_groups(i)=mass_groups(i)+mass_groups(ii)
mnew=mass(i)+mass(ii)
mass_err(i)=SQRT(((mass_err(i)**2)+(mass_err(ii)**2))/2.D0)
px(i)=(px(i)*mass(i)+px(ii)*mass(ii))/mnew
py(i)=(py(i)*mass(i)+py(ii)*mass(ii))/mnew
pz(i)=(pz(i)*mass(i)+pz(ii)*mass(ii))/mnew
mass(i)=mnew

rad(i)=((3.D0/(4.D0*PI)*mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0))/1.D6
 rad_err(i)=((10.D0**mass_err(i))-1.D0)*rad(i)/SQRT(3.D0)

ii=0
END IF
END IF
END IF

END DO


END IF
END DO

iii=0
DO i=1,n_all
IF (active(i)) THEN
iii=iii+1
END IF
END DO
n_new=iii
 WRITE(*,*) iii,'groups remain after density adjusted fi filtering'
 



DO i=1,n_rand
rand_x(i)=(RAND()*d_maxdepth*2.D0)-d_maxdepth
rand_y(i)=(RAND()*d_maxdepth*2.D0)-d_maxdepth
rand_z(i)=(RAND()*d_maxdepth*2.D0)-d_maxdepth
rand_act(i)=.TRUE.
rand_infi(i)=.FALSE.
rand_infi_multi(i)=.FALSE.
END DO

n_rand_used=n_rand
DO i=1,n_rand
d=(rand_x(i)**2)+(rand_y(i)**2)+(rand_z(i)**2)
IF (d>(d_maxdepth**2)) THEN
rand_act(i)=.FALSE.
n_rand_used=n_rand_used-1
END IF
END DO

DO i=1,n_all
IF (active(i)) THEN
DO ii=1,n_rand
IF (rand_act(ii)) THEN
d=((rand_x(ii)-px(i))**2)+((rand_y(ii)-py(i))**2)+((rand_z(ii)-pz(i))**2)
IF (d<(rad(i)**2)) THEN

IF (rand_infi(ii)) THEN
rand_infi_multi(i)=.TRUE.
rand_act(ii)=.FALSE.
END IF

rand_infi(ii)=.TRUE.
END IF
END IF
END DO
END IF
END DO



n_rand_multi=0
DO i=1,n_rand
IF (rand_infi_multi(i)) THEN
n_rand_multi=n_rand_multi+1
END IF
END DO

V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)
Vol_fi=(DBLE(n_rand_multi)/V_sum)/(DBLE(n_rand_used)/V_limit)

WRITE(*,*) Vol_fi,'of total survey volume covered by overlapping fi regions'



n_rand_fi=0
DO i=1,n_rand
IF (rand_infi(i)) THEN
n_rand_fi=n_rand_fi+1
END IF
END DO
V_limit=4.D0/3.D0*PI*((d_maxdepth*1.D6)**3)

Vol_fi=(DBLE(n_rand_fi)/V_sum)/(DBLE(n_rand_used)/V_limit)




DO i=1,n_all
IF (active(i)) THEN
 cdist_fi(i)=SQRT((px(i)**2)+(py(i)**2)+(pz(i)**2))
 IF (cdist_fi(i)>0.D0) THEN
 dec_fi(i)=ASIN(pz(i)/cdist_fi(i))
 ELSE
 dec_fi(i)=0.D0
 END IF
 ra_fi(i)=ATAN2(py(i),px(i))
 
 dec_fi(i)=dec_fi(i)*180.D0/PI
 ra_fi(i)=ra_fi(i)*180.D0/PI 
END IF
END DO
 
 
 
 



WRITE(*,*) Vol_fi,'of total volume occupied by fi regions after rescaling'

OPEN(50,file='ts_sor/catalogues/fi_info_fullmass.txt')
WRITE(50,*) 1.D0
WRITE(50,*) Vol_fi
 CLOSE(50)
 

 
 
 
 
 
OPEN(50,file='ts_sor/catalogues/fi_regions_fullmass.txt')
 DO i=1,n_all
IF (active(i)) THEN
WRITE(50,*) ra_fi(i),dec_fi(i),cdist_fi(i),px(i),py(i),pz(i),LOG10(mass(i)),mass_err(i),rad(i),rad_err(i)
END IF
END DO
 CLOSE(50)
 
  
 OPEN(50,file='final_catalogues/fi_regions_fullmass.txt')
 DO i=1,n_all
IF (active(i)) THEN
WRITE(50,"(2F14.8,4F10.3,4F10.5)") ra_fi(i),dec_fi(i),cdist_fi(i),px(i),py(i),pz(i),&
LOG10(mass(i)),mass_err(i),rad(i),rad_err(i)
END IF
END DO
 CLOSE(50)
 
 
 
 WRITE(*,*) 'all output written'
 
 
 
 
 
 
 
 
 
WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM













double precision FUNCTION det3d(matrix)

double precision, dimension(1:3,1:3) :: matrix
double precision :: helpsum

 helpsum=matrix(1,1)*matrix(2,2)*matrix(3,3)
helpsum=helpsum+matrix(2,1)*matrix(3,2)*matrix(1,3)
helpsum=helpsum+matrix(1,2)*matrix(2,3)*matrix(3,1)
helpsum=helpsum-matrix(1,3)*matrix(2,2)*matrix(3,1)
helpsum=helpsum-matrix(2,1)*matrix(1,2)*matrix(3,3)
helpsum=helpsum-matrix(2,3)*matrix(3,2)*matrix(1,1)

 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det3d'
 END IF 
 


 det3d=helpsum

! 
RETURN
END FUNCTION det3d

 
 
double precision FUNCTION det2d(matrix)

double precision, dimension(1:2,1:2) :: matrix
double precision :: helpsum

 helpsum=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)

 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det2d'
 END IF 
 
 det2d=helpsum

! 
RETURN
END FUNCTION det2d

 
