PROGRAM massfunction
! declaration of variables
IMPLICIT NONE

double precision :: PI,H0,q0,light,G,tophat,rho_crit,Omega_l,Omega_m
integer :: io_err,n,i,ii,iii,n_mock,n_halo,n_id,n_mock_all,n_id_all,n_help,help_index,help_count
character(200) :: filename_mock,filename_id,filename_halo,filename_mock_all,filename_id_all
double precision :: h,sdss_red_measure,weight,dummy_var,V_limit,v_sum_sdss
double precision :: angular_sep,divisor,dist_help,mag_sol_g,mhelpr,mhelpg,galred,m_expected
double precision :: t1,t2,t3,D_limit,z_limit,lum_sum_sdss,mabs_vollim,area_cover,vol,d_com,z_max

double precision, allocatable ::  ra(:),dec(:),z(:),mag_g(:),mag_r(:)
logical, allocatable :: hactive(:),mock_unused(:),mock_unused_all(:),iter_help(:),mock_act(:),active_sdss(:)
integer(kind=8), allocatable :: mmockid(:),mockid(:),mockfofid(:)
integer(kind=8), allocatable :: h_fofid(:),h_n(:)
double precision, allocatable ::  hx(:),hy(:),hz(:),h_mass(:),volume_weight_sdss(:),lum_gal_sdss(:)
! double precision :: sdss_cover,twomass_cover,mock_cover
double precision, allocatable ::  lumsum(:),mock_lum(:),comovdist(:),ang_distance(:),lum_dist_sdss(:)
double precision, allocatable ::  x_all(:),y_all(:),z_all(:),mag_g_all(:),mag_r_all(:),lum_int_sdss(:)
double precision, allocatable :: sigma_gap(:),sigma(:),v_order(:),vred(:),m_dyn(:),lumsum_cor(:)
double precision, allocatable :: radius(:),ang_rad(:),rad_order(:),dummy_core(:)
double precision, allocatable :: cluster_ra(:),cluster_dec(:),cluster_z(:),order_list(:)

double precision, allocatable ::  lumsum_all(:),mock_lum_all(:),comovdist_all(:)
integer(kind=8), allocatable :: mmockid_all(:),mockid_all(:),mockfofid_all(:)
integer(kind=8), allocatable :: n_visible(:),halo_id_of_mockgal(:),n_mock_group(:)
double precision, allocatable :: lumdist(:)
double precision :: magmax_sdss,magmin_sdss,magmax_2mrs,magmin_2mrs,mag_limit_sdss,mag_limit_2mrs


double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2,mag_sol_r
integer :: help_x,help_y,help_lum,uuu
integer, dimension(1:102) :: binclustermock,binclusterfof
integer, dimension(1:8) :: n_cores
double precision :: v1,v2,v3,v4,v5,v6,v7,v8

character(200) :: appendix
character(20) :: nm,ax
!parallization variables		
integer ierr, rankrank, nbnodes, namelen
!  character (len=MPI_MAX_PROCESSOR_NAME) :: name
integer, dimension(1:8) :: u_array
integer, dimension(1:1) :: u_local
double precision, dimension(1:8) :: S_array
double precision, dimension(1:1) :: S_local


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme MASS FUNCTION started'
WRITE(*,*) '============================================================'




! appendix='SDSS'
! mag_sol_r=4.71D0
! mag_sol_g=5.31D0     
! sdss_red_measure=30.D0
! mag_limit_sdss=17.77D0
! magmax_sdss=-30.D0
! magmin_sdss=-15.D0
z_max=0.11D0

appendix='2MRS'
mag_sol_r=3.28D0 !Ks
mag_sol_g=3.64D0 !J
sdss_red_measure=32.D0 !31.96D0
mag_limit_sdss=11.75D0
magmax_sdss=-30.D0
magmin_sdss=-18.D0


divisor=(SQRT(1.D0+2.D0*q0*z_max)+1.D0+q0*z_max)
d_com=light/H0*z_max*(1.D0+((z_max*(1.D0-q0))/divisor))
d_com=d_com/(1.D0+z_max)*1.D6

vol=4.D0*PI/3.D0*(d_com**3)

m_expected=vol*rho_crit*Omega_m

area_cover=1.D0


 Omega_m=0.25D0
 Omega_l=0.75D0
 H0=73.D0
 
 
! define constants
PI=ACOS(-1.D0)
 
q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)



! DO i=1,8
! u_array(i)=i
! END DO
! 
DO uuu=1,8
! call MPI_Scatter(u_array,1,MPI_INTEGER,u_local,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! 
! uuu=DBLE(u_local(1))
WRITE(nm,*) uuu
nm=adjustl(nm)
 
S_local=0.0D0

appendix=TRIM(adjustl(appendix))
ax=TRIM(adjustl(appendix))

filename_mock=TRIM(ax)//'_mock'//TRIM(nm)//'/'//TRIM(ax)//'_galaxies_mock'//TRIM(nm)//'_final.txt'
filename_id=TRIM(ax)//'_mock'//TRIM(nm)//'/'//TRIM(ax)//'_galaxies_mock'//TRIM(nm)//'_final_simid.txt'
filename_halo='halos'//TRIM(nm)//'/halos_mock'//TRIM(nm)//'_pos3D.txt'
filename_mock_all=TRIM(ax)//'_mock'//TRIM(nm)//'/'//TRIM(ax)//'_allgalaxies_mock'//TRIM(nm)//'_pos3D.txt'
filename_id_all=TRIM(ax)//'_mock'//TRIM(nm)//'/'//TRIM(ax)//'_allgalaxies_mock'//TRIM(nm)//'_id_file.txt'

! get length of file
OPEN(50,file=filename_mock)
io_err=0
n_mock=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock=n_mock+1
END DO
 CLOSE(50)
n_mock=n_mock-1

OPEN(50,file=filename_id)
io_err=0
n_id=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_id=n_id+1
END DO
 CLOSE(50)
n_id=n_id-1

IF (n_id.NE.n_mock) THEN
WRITE(*,*) 'WARNING: mock catalogue and ID do not match!!!'
END IF

WRITE(*,*) n_mock,'galaxies from mock catalogue will be used'

! get length of file
OPEN(50,file=filename_halo)
io_err=0
n_halo=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_halo=n_halo+1
END DO
 CLOSE(50)
n_halo=n_halo-1

WRITE(*,*) n_halo,'halos will be used'

OPEN(50,file=filename_mock_all)
io_err=0
n_mock_all=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock_all=n_mock_all+1
END DO
 CLOSE(50)
n_mock_all=n_mock_all-1

OPEN(50,file=filename_id_all)
io_err=0
n_id_all=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_id_all=n_id_all+1
END DO
 CLOSE(50)
n_id_all=n_id_all-1

IF (n_id_all.NE.n_mock_all) THEN
WRITE(*,*) 'WARNING: complete catalogue and ID do not match!!!'
END IF

WRITE(*,*) n_mock_all,'galaxies from the complete catalogue will be used'


allocate(ra(1:n_mock))
allocate(dec(1:n_mock))
allocate(z(1:n_mock))
allocate(mag_g(1:n_mock))
allocate(mag_r(1:n_mock))
allocate(volume_weight_sdss(1:n_mock))
allocate(lum_gal_sdss(1:n_mock))
allocate(mmockid(1:n_mock))
allocate(mockid(1:n_mock))
allocate(mockfofid(1:n_mock))
allocate(v_order(1:n_mock))
allocate(halo_id_of_mockgal(1:n_mock))
allocate(vred(1:n_mock))
allocate(iter_help(1:n_mock))
allocate(lum_dist_sdss(1:n_halo))
allocate(lum_int_sdss(1:n_halo))
allocate(radius(1:n_halo))
allocate(m_dyn(1:n_halo))
allocate(ang_rad(1:n_halo))
allocate(rad_order(1:n_mock))
allocate(cluster_ra(1:n_halo))
allocate(cluster_dec(1:n_halo))
allocate(cluster_z(1:n_halo))
allocate(n_mock_group(1:n_mock))
allocate(order_list(1:n_mock))
allocate(dummy_core(1:n_mock))
allocate(mock_act(1:n_mock))
allocate(ang_distance(1:n_halo))

allocate(x_all(1:n_mock_all))
allocate(y_all(1:n_mock_all))
allocate(z_all(1:n_mock_all))
allocate(mag_g_all(1:n_mock_all))
allocate(mag_r_all(1:n_mock_all))


allocate(mmockid_all(1:n_mock_all))
allocate(mockid_all(1:n_mock_all))
allocate(mockfofid_all(1:n_mock_all))

allocate(h_fofid(1:n_halo))
allocate(hx(1:n_halo))
allocate(hy(1:n_halo))
allocate(hz(1:n_halo))
allocate(h_mass(1:n_halo))
allocate(h_n(1:n_halo))

! read file
OPEN(50,file=filename_mock)
DO i=1,n_mock
READ(50,*) ra(i),dec(i),z(i),mag_g(i),mag_r(i)
END DO
CLOSE(50)


OPEN(50,file=filename_id)
DO i=1,n_mock
READ(50,*) mmockid(i),mockid(i),mockfofid(i)
END DO
CLOSE(50)
WRITE(*,*) 'all mock galaxies read in'

! read file
OPEN(50,file=filename_mock_all)
DO i=1,n_mock_all
READ(50,*) x_all(i),y_all(i),z_all(i),mag_g_all(i),mag_r_all(i)
END DO
CLOSE(50)

OPEN(50,file=filename_id_all)
DO i=1,n_mock_all
READ(50,*) mmockid_all(i),mockid_all(i),mockfofid_all(i)
END DO
CLOSE(50)
WRITE(*,*) 'all galaxies read in'

OPEN(50,file=filename_halo)
DO i=1,n_halo
READ(50,*) h_fofid(i),hx(i),hy(i),hz(i),h_mass(i),h_n(i)
END DO
CLOSE(50)
WRITE(*,*) 'all halos read in'
WRITE(*,*) 'all files read in'

allocate(hactive(1:n_halo))
allocate(mock_unused(1:n_mock))
allocate(lumsum(1:n_halo))
allocate(lumsum_cor(1:n_halo))
allocate(n_visible(1:n_halo))
allocate(mock_lum(1:n_mock))
allocate(comovdist(1:n_halo))
allocate(lumdist(1:n_halo))
allocate(sigma_gap(1:n_halo))
allocate(sigma(1:n_halo))
allocate(active_sdss(1:n_mock))
allocate(lumsum_all(1:n_halo))
allocate(mock_lum_all(1:n_mock_all))
allocate(mock_unused_all(1:n_mock_all))

DO i=1,n_mock
ra(i)=ra(i)*PI/180.D0
dec(i)=dec(i)*PI/180.D0
END DO

DO i=1,n_halo
hactive(i)=.FALSE.
lumsum(i)=0.D0
lumsum_all(i)=0.D0
n_visible(i)=0
 comovdist(i)=SQRT((hx(i)**2)+(hy(i)**2)+(hz(i)**2))
 
 
 !redshift
galred=(light**2)*((light**2)+(comovdist(i)**2)*(H0**2)*(1.D0-2.D0*q0))*(q0-1.D0)**2
galred=SQRT(galred)
galred=galred+light*comovdist(i)*H0+(light**2)*(q0-1.D0)-(comovdist(i)**2)*(H0**2)*(q0**2)
galred=galred/((light-comovdist(i)*H0*q0)**2)


!luminosity distance
lumdist(i)=comovdist(i)*(1.D0+galred)

 
!  
 sigma(i)=0.D0
 sigma_gap(i)=0.D0
radius(i)=0.D0
ang_rad(i)=0.D0
 cluster_ra(i)=0.D0
 cluster_dec(i)=0.D0
 cluster_z(i)=0.D0
 
END DO

DO i=1,n_mock
mock_lum(i)=10.D0**(-0.4D0*(mag_r(i)-mag_sol_r))
! lum_g(i)=10.D0**(-0.4D0*(mag_g(i)-mag_sol_g))
mock_unused(i)=.TRUE.
n_mock_group(i)=0
vred(i)=((((1.D0+z(i))**2)-1.D0)/(((1.D0+z(i))**2)+1.D0))*light
END DO

DO i=1,n_mock_all
mock_lum_all(i)=10.D0**(-0.4D0*(mag_r_all(i)-mag_sol_r))
mock_unused_all(i)=.TRUE.
END DO




DO i=1,n_halo
WRITE(*,*) i,n_halo,uuu
DO ii=1,n_mock
IF (mock_unused(ii)) THEN
IF (h_fofid(i)==mockfofid(ii)) THEN
hactive(i)=.TRUE.
mock_unused(ii)=.FALSE.
n_visible(i)=n_visible(i)+1
lumsum(i)=lumsum(i)+mock_lum(ii)
halo_id_of_mockgal(ii)=i
END IF
END IF
END DO

 
END DO
 

 
 
 WRITE(*,*) 'all mock galaxies and halos correlated'
 
 DO i=1,n_halo
IF (hactive(i)) THEN
DO ii=1,n_mock_all
IF (mock_unused_all(ii)) THEN
IF (h_fofid(i)==mockfofid_all(ii)) THEN

mock_unused_all(ii)=.FALSE.

lumsum_all(i)=lumsum_all(i)+mock_lum_all(ii)

END IF
END IF
END DO

 END IF
END DO
! 



! 




! calculate volume weights to corret for Malmquist bias
V_sum_sdss=0.D0
lum_sum_sdss=0.D0
DO i=1,n_mock
lum_gal_sdss(i)=mock_lum(i)
active_sdss(i)=.TRUE.
volume_weight_sdss(i)=0.D0

IF (mag_r(i)<magmax_sdss) THEN
active_sdss(i)=.FALSE.
END IF

IF (mag_r(i)>magmin_sdss) THEN
active_sdss(i)=.FALSE.
END IF

IF (mag_g(i)<(magmax_sdss-1.5D0)) THEN
active_sdss(i)=.FALSE.
END IF

IF (mag_g(i)>(magmin_sdss+1.5D0)) THEN
active_sdss(i)=.FALSE.
END IF


IF ((mag_g(i)-mag_r(i))>2.D0) THEN
active_sdss(i)=.FALSE.
END IF

IF ((mag_g(i)-mag_r(i))<-1.D0) THEN
active_sdss(i)=.FALSE.
END IF


IF (active_sdss(i)) THEN


D_limit=-0.2D0*mag_r(i)+((mag_limit_sdss+5.D0)/5.D0)
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
! WRITE(*,*) 'lumgal + weight',
END IF

END DO 

WRITE(*,*) lum_sum_sdss

DO i=1,n_mock
IF (active_sdss(i)) THEN
lum_gal_sdss(i)=lum_gal_sdss(i)/lum_sum_sdss
END IF
END DO



WRITE(*,*) 'weights calculated'
!  

DO i=1,n_halo
IF (lumdist(i).LE.0.D0) THEN
hactive(i)=.FALSE.
END IF
END DO

! 
DO i=1,n_halo
IF (hactive(i)) THEN
lum_dist_sdss(i)=LOG10(lumdist(i))+6.D0
lum_int_sdss(i)=0.D0
mabs_vollim=mag_limit_sdss+5.D0-5.D0*lum_dist_sdss(i)

DO ii=1,n_mock
IF (active_sdss(ii)) THEN
IF (mag_r(ii)<mabs_vollim) THEN
lum_int_sdss(i)=lum_int_sdss(i)+lum_gal_sdss(ii)
END IF
END IF
END DO

lumsum_cor(i)=lumsum(i)/lum_int_sdss(i)

END IF
END DO







 DO i=1,n_halo
IF (hactive(i)) THEN
IF (n_visible(i)>1) THEN
DO ii=1,n_mock
v_order(ii)=0.D0
END DO
n_help=n_visible(i)
iii=0
DO ii=1,n_mock
IF (halo_id_of_mockgal(ii)==i) THEN
iii=iii+1
v_order(iii)=vred(ii)
END IF
END DO

IF (iii.NE.n_help) THEN
WRITE(*,*) 'problem with halo ',i,'with',iii,n_help
END IF

DO ii=1,n_help
DO iii=ii,n_help
IF (v_order(ii)<v_order(iii)) THEN
dummy_var=v_order(ii)
v_order(ii)=v_order(iii)
v_order(iii)=dummy_var
END IF
END DO
END DO

sigma_gap(i)=0.D0

DO ii=1,(n_help-1)
weight=DBLE(ii*(n_help-ii))
sigma_gap(i)=sigma_gap(i)+weight*ABS(v_order(ii+1)-v_order(ii))
END DO

sigma_gap(i)=sigma_gap(i)*SQRT(PI)/DBLE(n_help*(n_help-1))

IF (sigma_gap(i)<sdss_red_measure) THEN
sigma_gap(i)=sdss_red_measure
END IF

sigma(i)=SQRT(DBLE(n_help)/DBLE(n_help-1)*(sigma_gap(i)**2)-(sdss_red_measure**2))



IF (sigma(i)<sdss_red_measure) THEN
sigma(i)=sdss_red_measure

END IF


END IF
END IF
END DO
!  



 DO i=1,n_halo
IF (hactive(i)) THEN
IF (n_visible(i)>1) THEN

DO ii=1,n_mock
order_list(ii)=0.D0
END DO


n_help=n_visible(i)
iii=0
DO ii=1,n_mock
iter_help(ii)=.FALSE.
mock_act(ii)=.FALSE.
IF (halo_id_of_mockgal(ii)==i) THEN
iii=iii+1
order_list(iii)=z(ii)
iter_help(ii)=.TRUE.
mock_act(ii)=.TRUE.
END IF
END DO
! WRITE(*,*) iii,n_visible(i)

DO ii=1,n_help
DO iii=ii,n_help
IF (order_list(ii)<order_list(iii)) THEN
dummy_var=order_list(ii)
order_list(ii)=order_list(iii)
order_list(iii)=dummy_var
END IF
END DO
END DO


IF (MOD(n_help,2)==0) THEN
help_index=n_help/2
 cluster_z(i)=(order_list(help_index)+order_list(help_index+1))/2.D0
ELSE
help_index=n_help/2
 cluster_z(i)=order_list(help_index+1)
END IF



help_count=n_help







DO WHILE (help_count>1) 
 cluster_ra(i)=0.D0
 cluster_dec(i)=0.D0
help_lum=0.D0

DO ii=1,n_mock
IF (iter_help(ii)) THEN
 cluster_ra(i)=cluster_ra(i)+(ra(ii)*mock_lum(ii))
 cluster_dec(i)=cluster_dec(i)+(dec(ii)*mock_lum(ii))
help_lum=help_lum+mock_lum(ii)
END IF
END DO
 cluster_ra(i)=cluster_ra(i)/help_lum
 cluster_dec(i)=cluster_dec(i)/help_lum



DO ii=1,n_mock
IF (iter_help(ii)) THEN
angular_sep=COS(cluster_dec(i))*COS(dec(ii))*COS((ra(ii)-cluster_ra(i)))
angular_sep=angular_sep+SIN(cluster_dec(i))*SIN(dec(ii))
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

dummy_core(ii)=angular_sep

END IF
END DO

angular_sep=0.D0
help_index=0
DO ii=1,n_mock
IF (iter_help(ii)) THEN
IF (angular_sep<dummy_core(ii)) THEN
angular_sep=dummy_core(ii)
help_index=ii
END IF
END IF
END DO
iter_help(help_index)=.FALSE.

help_count=0
DO ii=1,n_mock
IF (iter_help(ii)) THEN
help_count=help_count+1
END IF
END DO


END DO

DO ii=1,n_mock
IF (iter_help(ii)) THEN
 cluster_ra(i)=ra(ii)
 cluster_dec(i)=dec(ii)
END IF
END DO









help_count=0
DO ii=1,n_mock
IF (mock_act(ii)) THEN
help_count=help_count+1

angular_sep=COS(cluster_dec(i))*COS(dec(ii))*COS((ra(ii)-cluster_ra(i)))
angular_sep=angular_sep+SIN(cluster_dec(i))*SIN(dec(ii))
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

IF (ISNAN(angular_sep)) THEN
WRITE(*,*) cluster_ra(i),cluster_dec(i),ra(ii),dec(ii)
END IF

dummy_core(help_count)=angular_sep
END IF
END DO

IF (help_count.NE.n_help) THEN
WRITE(*,*) 'problem with halo ',i,'with',help_count,n_visible(i)
END IF



DO ii=1,help_count
DO iii=ii,help_count
IF (dummy_core(ii)<dummy_core(iii)) THEN
dummy_var=dummy_core(ii)
dummy_core(ii)=dummy_core(iii)
dummy_core(iii)=dummy_var
END IF
END DO
END DO



IF (MOD(help_count,2)==0) THEN
help_index=help_count/2
ang_rad(i)=(dummy_core(help_index)+dummy_core(help_index+1))/2.D0
ELSE
help_index=help_count/2
ang_rad(i)=dummy_core(help_index+1)
END IF

divisor=(SQRT(1.D0+2.D0*q0*cluster_z(i))+1.D0+q0*cluster_z(i))
!luminosity distance
ang_distance(i)=light/H0*cluster_z(i)*(1.D0+((cluster_z(i)*(1.D0-q0))/divisor))


ang_distance(i)=ang_distance(i)*((1.D0+cluster_z(i))**(-2.D0))*1000.D0  !angular diameter distance in kpc


radius(i)=TAN(ang_rad(i))*ang_distance(i) ! radius in kpc






END IF
END IF
END DO
!  
 DO i=1,n_halo
IF (hactive(i)) THEN
IF (n_visible(i)>1) THEN
m_dyn(i)=10.D0/G*(3.D0*(sigma(i)**2))*radius(i) !WURZEL(3)*sigma_rad=sigma_3d
END IF
END IF
END DO


  WRITE(*,*) 'all galaxies and halos correlated'
 
 
   
OPEN(51,file='mass/all_mock/log_dep_single_'//TRIM(ax)//TRIM(nm)//'.txt')
DO i=1,n_halo
IF (hactive(i)) THEN
IF (n_visible(i)==1) THEN
WRITE(51,*) LOG10(h_mass(i)),LOG10(lumsum_cor(i)),LOG10(lumdist(i)),LOG10(lumsum(i)),LOG10(lumsum_all(i))
END IF
END IF
END DO
CLOSE(51)

 
OPEN(51,file='mass/all_mock/log_dep_multi_high_'//TRIM(ax)//TRIM(nm)//'.txt')
DO i=1,n_halo
IF (hactive(i)) THEN
IF (n_visible(i)>4) THEN
WRITE(51,*) LOG10(h_mass(i)),LOG10(lumsum_cor(i)),LOG10(lumdist(i)),LOG10(m_dyn(i)),&
LOG10(sigma(i)),LOG10(radius(i)),LOG10(DBLE(n_visible(i))),LOG10(lumsum(i)),LOG10(lumsum_all(i))
END IF
END IF
END DO
CLOSE(51)
 
 
 OPEN(51,file='mass/all_mock/log_dep_multi_low_'//TRIM(ax)//TRIM(nm)//'.txt')
DO i=1,n_halo
IF (hactive(i)) THEN
IF ((n_visible(i)>1).AND.(n_visible(i)<5)) THEN
WRITE(51,*) LOG10(h_mass(i)),LOG10(lumsum_cor(i)),LOG10(lumdist(i)),LOG10(m_dyn(i)),&
LOG10(sigma(i)),LOG10(radius(i)),LOG10(DBLE(n_visible(i))),LOG10(lumsum(i)),LOG10(lumsum_all(i))
END IF
END IF
END DO
CLOSE(51)


deallocate(ra)
deallocate(dec)
deallocate(z)
deallocate(mag_g)
deallocate(mag_r)
deallocate(volume_weight_sdss)
deallocate(lum_gal_sdss)
deallocate(mmockid)
deallocate(mockid)
deallocate(mockfofid)
deallocate(v_order)
deallocate(halo_id_of_mockgal)
deallocate(vred)
deallocate(iter_help)
deallocate(lum_dist_sdss)
deallocate(lum_int_sdss)
deallocate(radius)
deallocate(m_dyn)
deallocate(ang_rad)
deallocate(rad_order)
deallocate(cluster_ra)
deallocate(cluster_dec)
deallocate(cluster_z)
deallocate(n_mock_group)
deallocate(order_list)
deallocate(dummy_core)
deallocate(mock_act)
deallocate(ang_distance)
deallocate(x_all)
deallocate(y_all)
deallocate(z_all)
deallocate(mag_g_all)
deallocate(mag_r_all)
deallocate(mmockid_all)
deallocate(mockid_all)
deallocate(mockfofid_all)
deallocate(h_fofid)
deallocate(hx)
deallocate(hy)
deallocate(hz)
deallocate(h_mass)
deallocate(h_n)
deallocate(hactive)
deallocate(mock_unused)
deallocate(lumsum)
deallocate(lumsum_cor)
deallocate(n_visible)
deallocate(mock_lum)
deallocate(comovdist)
deallocate(lumdist)
deallocate(sigma_gap)
deallocate(sigma)
deallocate(active_sdss)
deallocate(lumsum_all)
deallocate(mock_lum_all)
deallocate(mock_unused_all)

END DO
 
WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


 
 
 
 
 