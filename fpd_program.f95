PROGRAM fundamentalplanedistances
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit,H0_m,H0_s,H0_p
integer(kind=8) :: io_err,n,i,ii,n_clusters,n_galaxies,n_fp,n_fp_p,n_fp_s,n_fp_m,n_ell_used,help0,n_cluster_used

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,help1,help2,help3,divisor,av_ratio
double precision :: fp_rms_p,fp_rms_s,fp_rms_m,err_basis_p,err_basis_s,err_basis_m,z_err
logical, allocatable :: act_fp_m(:),act_fp_s(:),act_fp_p(:),act_fp(:),act_cluster(:)
double precision, allocatable :: ra_clustercore(:),dec_clustercore(:),z_clustercore(:),lumtot_clustercore(:)
double precision, allocatable :: sigma_clustercore(:),dummy_core(:),radius_clustercore(:),angradius_clustercore(:)
double precision, allocatable :: lumdistance_clustercore(:),mass_clustercore(:),colour_clustercore(:)
double precision, allocatable :: lumotherband(:),lum_g(:),sigma_gap(:),lumobs_clustercore(:) 
double precision, allocatable :: mass_err(:),lumtot_err(:) 
double precision, allocatable :: mdyn_clustercore(:) 
integer(kind=8), allocatable :: members_clustercore(:),clusterindex(:),clustergalnmember(:)
integer(kind=8), allocatable :: SDSS_id_gal(:),numberofgal(:),fp_clusterindex(:),objid(:)
double precision, allocatable :: fp_adist_p(:),fp_cdist_p(:),fp_ldist_p(:)
double precision, allocatable :: fp_adist_s(:),fp_cdist_s(:),fp_ldist_s(:)
double precision, allocatable :: fp_adist_m(:),fp_cdist_m(:),fp_ldist_m(:)
double precision, allocatable :: dist_err_m(:),dist_err_s(:),dist_err_p(:)


double precision, allocatable :: cluster_fp_adist_p(:),cluster_fp_cdist_p(:),cluster_fp_ldist_p(:)
double precision, allocatable :: cluster_fp_adist_s(:),cluster_fp_cdist_s(:),cluster_fp_ldist_s(:)
double precision, allocatable :: cluster_fp_adist_m(:),cluster_fp_cdist_m(:),cluster_fp_ldist_m(:)
double precision, allocatable :: cluster_z_adist_p(:),cluster_z_cdist_p(:),cluster_z_ldist_p(:)
double precision, allocatable :: cluster_z_adist_s(:),cluster_z_cdist_s(:),cluster_z_ldist_s(:)
double precision, allocatable :: cluster_z_adist_m(:),cluster_z_cdist_m(:),cluster_z_ldist_m(:)
double precision, allocatable :: err_fpd_m(:),err_fpd_s(:),err_fpd_p(:)
double precision, allocatable :: err_red_m(:),err_red_s(:),err_red_p(:)
double precision, allocatable :: err_z_m(:),err_z_s(:),err_z_p(:)

integer(kind=8), allocatable :: n_fp_in_cluster(:)

double precision, allocatable ::  angular_dist(:),ra_rad(:),dec_rad(:),lum(:),lum_dist(:)

double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y
integer, dimension(1:102) :: bincluster
double precision , dimension(1:100) :: av_m,av_s,av_p,sigma_m,sigma_s,sigma_p



WRITE(*,*) '============================================================'
WRITE(*,*) '    programme FUNDAMENTAL PLANE DISTANCES started'
WRITE(*,*) '============================================================'



! define constants
PI=ACOS(-1.D0)
G=4.302D-3 !in pc/Msol * (km/s)**2
light=3.D5
H0_m=73.D0
H0_s=70.D0
H0_p=67.3D0
z_err=30.D0/light

! get length of file
OPEN(50,file='catalogues/cluster_list_all_SDSS.txt')
io_err=0
n_clusters=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_clusters=n_clusters+1
END DO
 CLOSE(50)
n_clusters=n_clusters-1

WRITE(*,*) n_clusters,'clusters will be used'



! get length of file
OPEN(50,file='catalogues/cluster_members_SDSS.txt')
io_err=0
n_galaxies=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_galaxies=n_galaxies+1
END DO
 CLOSE(50)
n_galaxies=n_galaxies-1

WRITE(*,*) n_galaxies,'galaxies will be used'

OPEN(61,file='fp_distances/n_max.txt')
READ(61,*) n_fp
 CLOSE(61)


 OPEN(50,file='fp_distances/fp_z_planck.txt')
io_err=0
n_fp_p=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_fp_p=n_fp_p+1
END DO
 CLOSE(50)
n_fp_p=n_fp_p-2
 
  OPEN(50,file='fp_distances/fp_z_simple.txt')
io_err=0
n_fp_s=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_fp_s=n_fp_s+1
END DO
 CLOSE(50)
n_fp_s=n_fp_s-2
 
  OPEN(50,file='fp_distances/fp_z_millennium.txt')
io_err=0
n_fp_m=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_fp_m=n_fp_m+1
END DO
 CLOSE(50)
n_fp_m=n_fp_m-2
 

allocate(fp_adist_p(1:n_fp))
allocate(fp_cdist_p(1:n_fp))
allocate(fp_ldist_p(1:n_fp))

allocate(fp_adist_s(1:n_fp))
allocate(fp_cdist_s(1:n_fp))
allocate(fp_ldist_s(1:n_fp))

allocate(fp_adist_m(1:n_fp))
allocate(fp_cdist_m(1:n_fp))
allocate(fp_ldist_m(1:n_fp))
 
allocate(cluster_fp_adist_p(1:n_clusters))
allocate(cluster_fp_cdist_p(1:n_clusters))
allocate(cluster_fp_ldist_p(1:n_clusters))

allocate(cluster_fp_adist_s(1:n_clusters))
allocate(cluster_fp_cdist_s(1:n_clusters))
allocate(cluster_fp_ldist_s(1:n_clusters))

allocate(cluster_fp_adist_m(1:n_clusters))
allocate(cluster_fp_cdist_m(1:n_clusters))
allocate(cluster_fp_ldist_m(1:n_clusters))

allocate(cluster_z_adist_p(1:n_clusters))
allocate(cluster_z_cdist_p(1:n_clusters))
allocate(cluster_z_ldist_p(1:n_clusters))

allocate(cluster_z_adist_s(1:n_clusters))
allocate(cluster_z_cdist_s(1:n_clusters))
allocate(cluster_z_ldist_s(1:n_clusters))

allocate(cluster_z_adist_m(1:n_clusters))
allocate(cluster_z_cdist_m(1:n_clusters))
allocate(cluster_z_ldist_m(1:n_clusters))

allocate(n_fp_in_cluster(1:n_clusters))
allocate(act_cluster(1:n_clusters))
 
allocate(dist_err_m(1:n_clusters))
allocate(dist_err_s(1:n_clusters))
allocate(dist_err_p(1:n_clusters))

 
 allocate(act_fp_m(1:n_fp))
 allocate(act_fp_s(1:n_fp)) 
 allocate(act_fp_p(1:n_fp)) 
 
  
 allocate(act_fp(1:n_fp))
 allocate(objid(1:n_fp))

 allocate(fp_clusterindex(1:n_fp))


allocate(ra_clustercore(1:n_clusters)) 
allocate(dec_clustercore(1:n_clusters)) 
allocate(z_clustercore(1:n_clusters)) 
allocate(lumtot_clustercore(1:n_clusters)) 
allocate(mass_clustercore(1:n_clusters)) 
allocate(sigma_clustercore(1:n_clusters)) 
allocate(radius_clustercore(1:n_clusters)) 
allocate(angradius_clustercore(1:n_clusters)) 
allocate(lumdistance_clustercore(1:n_clusters)) 
allocate(colour_clustercore(1:n_clusters)) 
allocate(members_clustercore(1:n_clusters)) 
allocate(lumobs_clustercore(1:n_clusters)) 
allocate(mdyn_clustercore(1:n_clusters)) 
allocate(mass_err(1:n_clusters)) 
allocate(lumtot_err(1:n_clusters)) 
allocate(numberofgal(1:n_galaxies)) 
allocate(clusterindex(1:n_galaxies)) 
allocate(SDSS_id_gal(1:n_galaxies)) 

allocate(err_fpd_m(1:n_clusters)) 
allocate(err_fpd_s(1:n_clusters)) 
allocate(err_fpd_p(1:n_clusters)) 
allocate(err_red_m(1:n_clusters)) 
allocate(err_red_s(1:n_clusters)) 
allocate(err_red_p(1:n_clusters)) 
allocate(err_z_m(1:n_clusters)) 
allocate(err_z_s(1:n_clusters)) 
allocate(err_z_p(1:n_clusters)) 



DO i=1,n_fp
act_fp_m(i)=.FALSE.
act_fp_s(i)=.FALSE.
act_fp_p(i)=.FALSE.
END DO
 

OPEN(61,file='fp_distances/fp_z_planck.txt')
READ(61,*) 
DO i=1,n_fp_p
READ(61,*) ii,help0,help1,help2,help3
objid(ii)=help0
fp_adist_p(ii)=help1
fp_cdist_p(ii)=help2
fp_ldist_p(ii)=help3
act_fp_p(ii)=.TRUE.
END DO
 CLOSE(61)

OPEN(61,file='fp_distances/fp_z_simple.txt')
READ(61,*) 
DO i=1,n_fp_s
READ(61,*) ii,help0,help1,help2,help3
objid(ii)=help0
fp_adist_s(ii)=help1
fp_cdist_s(ii)=help2
fp_ldist_s(ii)=help3
act_fp_s(ii)=.TRUE.
END DO
 CLOSE(61)
 
OPEN(61,file='fp_distances/fp_z_millennium.txt')
READ(61,*) 
DO i=1,n_fp_m
READ(61,*) ii,help0,help1,help2,help3
objid(ii)=help0
fp_adist_m(ii)=help1
fp_cdist_m(ii)=help2
fp_ldist_m(ii)=help3
act_fp_m(ii)=.TRUE.
END DO
 CLOSE(61) 
 

 
OPEN(50,file='catalogues/cluster_list_all_SDSS.txt')
DO i=1,n_clusters
READ(50,*) help0,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
lumtot_clustercore(i),lumtot_err(i),lumobs_clustercore(i),&
mass_clustercore(i),mass_err(i),mdyn_clustercore(i),&
sigma_clustercore(i),radius_clustercore(i),angradius_clustercore(i),lumdistance_clustercore(i),members_clustercore(i)
END DO
 CLOSE(50)

 
OPEN(50,file='catalogues/cluster_members_SDSS.txt')
DO i=1,n_galaxies
   READ(50,*) numberofgal(i),clusterindex(i),SDSS_id_gal(i) 
END DO
 CLOSE(50)
 
 OPEN(50,file='fp_rms.txt')
   READ(50,*) fp_rms_p
   READ(50,*) fp_rms_s
   READ(50,*) fp_rms_m
 CLOSE(50)
 
 err_basis_p=(10.D0**fp_rms_p)-1.D0
 err_basis_s=(10.D0**fp_rms_s)-1.D0
 err_basis_m=(10.D0**fp_rms_m)-1.D0
 
 
 
 DO i=1,n_fp
 act_fp(i)=.TRUE.
 IF (act_fp_m(i).EQV..FALSE.) THEN
 act_fp(i)=.FALSE. 
 END IF
 IF (act_fp_s(i).EQV..FALSE.) THEN
 act_fp(i)=.FALSE. 
 END IF
 IF (act_fp_p(i).EQV..FALSE.) THEN
 act_fp(i)=.FALSE. 
 END IF
 END DO

  DO i=1,n_fp 
 IF (act_fp(i)) THEN
 DO ii=1,n_fp
 IF (act_fp(i)) THEN
 IF (i.NE.ii) THEN
 
 IF (objid(i)==objid(ii)) THEN
 act_fp(ii)=.FALSE.
 END IF
 
 END IF
 END IF
 END DO
 END IF
 END DO
 
 
 
 
 
 
 
 n_ell_used=0
 
  DO i=1,n_fp
 IF (act_fp(i)) THEN
n_ell_used=n_ell_used+1
 END IF
 END DO
WRITE(*,*) n_ell_used,'elliptical galaxies used'
 
 DO i=1,n_fp
 IF (act_fp(i)) THEN 
 fp_clusterindex(i)=0
 DO ii=1,n_galaxies
 IF (SDSS_id_gal(ii)==objid(i)) THEN 
 fp_clusterindex(i)=clusterindex(ii)
 END IF 
 END DO
 
 IF (fp_clusterindex(i)==0) THEN
 act_fp(i)=.FALSE. 
 END IF
 
 END IF  
 END DO
 
 n_ell_used=0
   DO i=1,n_fp
 IF (act_fp(i)) THEN
n_ell_used=n_ell_used+1
 END IF
 END DO
WRITE(*,*) n_ell_used,'elliptical galaxies in the cluster catalogue'
 
 
 
 
 DO i=1,n_clusters
 act_cluster(i)=.TRUE.
n_fp_in_cluster(i)=0 
 cluster_fp_adist_p(i)=0.D0
 cluster_fp_cdist_p(i)=0.D0
 cluster_fp_ldist_p(i)=0.D0
 cluster_fp_adist_s(i)=0.D0
 cluster_fp_cdist_s(i)=0.D0
 cluster_fp_ldist_s(i)=0.D0
 cluster_fp_adist_m(i)=0.D0
 cluster_fp_cdist_m(i)=0.D0
 cluster_fp_ldist_m(i)=0.D0

 DO ii=1,n_fp
 IF (act_fp(ii)) THEN
 IF (fp_clusterindex(ii)==i) THEN
n_fp_in_cluster(i)=n_fp_in_cluster(i)+1 
 
 cluster_fp_adist_p(i)=cluster_fp_adist_p(i)+fp_adist_p(ii)
 cluster_fp_cdist_p(i)=cluster_fp_cdist_p(i)+fp_cdist_p(ii)
 cluster_fp_ldist_p(i)=cluster_fp_ldist_p(i)+fp_ldist_p(ii)
 cluster_fp_adist_s(i)=cluster_fp_adist_s(i)+fp_adist_s(ii)
 cluster_fp_cdist_s(i)=cluster_fp_cdist_s(i)+fp_cdist_s(ii)
 cluster_fp_ldist_s(i)=cluster_fp_ldist_s(i)+fp_ldist_s(ii)
 cluster_fp_adist_m(i)=cluster_fp_adist_m(i)+fp_adist_m(ii)
 cluster_fp_cdist_m(i)=cluster_fp_cdist_m(i)+fp_cdist_m(ii)
 cluster_fp_ldist_m(i)=cluster_fp_ldist_m(i)+fp_ldist_m(ii)
 
 END IF
 END IF
 END DO
 
IF (n_fp_in_cluster(i)==0) THEN
 act_cluster(i)=.FALSE.
 
ELSE
 
 cluster_fp_adist_p(i)=cluster_fp_adist_p(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_cdist_p(i)=cluster_fp_cdist_p(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_ldist_p(i)=cluster_fp_ldist_p(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_adist_s(i)=cluster_fp_adist_s(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_cdist_s(i)=cluster_fp_cdist_s(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_ldist_s(i)=cluster_fp_ldist_s(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_adist_m(i)=cluster_fp_adist_m(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_cdist_m(i)=cluster_fp_cdist_m(i)/DBLE(n_fp_in_cluster(i))
 cluster_fp_ldist_m(i)=cluster_fp_ldist_m(i)/DBLE(n_fp_in_cluster(i))
 
END IF
 
 END DO

 
 ! Planck cosmology redshift distances
 Omega_m=0.315D0
 Omega_l=0.685D0
 H0=67.3D0
  q0=Omega_m/2.D0-Omega_l
light=3.D5

     DO i=1,n_clusters
 IF (act_cluster(i)) THEN
 divisor=(SQRT(1.D0+2.D0*q0*z_clustercore(i))+1.D0+q0*z_clustercore(i))
 cluster_z_ldist_p(i)=light/H0*z_clustercore(i)*(1.D0+((z_clustercore(i)*(1.D0-q0))/divisor))
 cluster_z_cdist_p(i)=cluster_z_ldist_p(i)*((1.D0+z_clustercore(i))**(-1.D0))
 cluster_z_adist_p(i)=cluster_z_ldist_p(i)*((1.D0+z_clustercore(i))**(-2.D0)) 
 END IF
 END DO
 
 
 ! simple cosmology redshift distances
 Omega_m=0.3D0
 Omega_l=0.7D0
 H0=70.D0
  q0=Omega_m/2.D0-Omega_l
light=3.D5

     DO i=1,n_clusters
 IF (act_cluster(i)) THEN
 divisor=(SQRT(1.D0+2.D0*q0*z_clustercore(i))+1.D0+q0*z_clustercore(i))
 cluster_z_ldist_s(i)=light/H0*z_clustercore(i)*(1.D0+((z_clustercore(i)*(1.D0-q0))/divisor))
 cluster_z_cdist_s(i)=cluster_z_ldist_s(i)*((1.D0+z_clustercore(i))**(-1.D0))
 cluster_z_adist_s(i)=cluster_z_ldist_s(i)*((1.D0+z_clustercore(i))**(-2.D0)) 
 END IF
 END DO 
 
 
 
  ! millennium cosmology redshift distances
 Omega_m=0.25D0
 Omega_l=0.75D0
 H0=73.D0
  q0=Omega_m/2.D0-Omega_l
light=3.D5

     DO i=1,n_clusters
 IF (act_cluster(i)) THEN
 divisor=(SQRT(1.D0+2.D0*q0*z_clustercore(i))+1.D0+q0*z_clustercore(i))
 cluster_z_ldist_m(i)=light/H0*z_clustercore(i)*(1.D0+((z_clustercore(i)*(1.D0-q0))/divisor))
 cluster_z_cdist_m(i)=cluster_z_ldist_m(i)*((1.D0+z_clustercore(i))**(-1.D0))
 cluster_z_adist_m(i)=cluster_z_ldist_m(i)*((1.D0+z_clustercore(i))**(-2.D0)) 
 END IF
 END DO 
 
 
 
 n_cluster_used=0
    DO i=1,n_clusters
 IF (act_cluster(i)) THEN
n_cluster_used=n_cluster_used+1
 END IF
 END DO
WRITE(*,*) n_cluster_used,'clusters with elliptical galaxies'
 
 
 
 
  n_cluster_used=0
    DO i=1,n_clusters
 IF (act_cluster(i)) THEN
 IF (n_fp_in_cluster(i)>1) THEN
n_cluster_used=n_cluster_used+1
END IF
 END IF
 END DO
WRITE(*,*) n_cluster_used,'clusters with 2 or more elliptical galaxies'
 
 
 
   n_cluster_used=0
    DO i=1,n_clusters
 IF (act_cluster(i)) THEN
 IF (n_fp_in_cluster(i)>2) THEN
n_cluster_used=n_cluster_used+1
END IF
 END IF
 END DO
WRITE(*,*) n_cluster_used,'clusters with 3 or more elliptical galaxies'
 
 
 
  
   n_cluster_used=0
    DO i=1,n_clusters
 IF (act_cluster(i)) THEN
 IF (n_fp_in_cluster(i)>4) THEN
n_cluster_used=n_cluster_used+1
END IF
 END IF
 END DO
WRITE(*,*) n_cluster_used,'clusters with 5 or more elliptical galaxies'
 
 
 
DO i=1,n_clusters
IF (act_cluster(i)) THEN
err_fpd_m(i)=err_basis_m/SQRT(DBLE(n_fp_in_cluster(i)))
err_fpd_s(i)=err_basis_s/SQRT(DBLE(n_fp_in_cluster(i)))
err_fpd_p(i)=err_basis_p/SQRT(DBLE(n_fp_in_cluster(i)))



err_z_m(i)=z_err/SQRT(DBLE(members_clustercore(i)))
err_z_s(i)=z_err/SQRT(DBLE(members_clustercore(i)))
err_z_p(i)=z_err/SQRT(DBLE(members_clustercore(i)))

err_red_m(i)=err_z_m(i)/z_clustercore(i)
err_red_s(i)=err_z_s(i)/z_clustercore(i)
err_red_p(i)=err_z_p(i)/z_clustercore(i)


END IF
END DO 
 
 
 
 
 
 
 
 
 
OPEN(50,file='catalogues/cluster_fp_distances_SDSS_p.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),err_z_p(i),&
 cluster_fp_adist_p(i),(err_fpd_p(i)*cluster_fp_adist_p(i)),&
 cluster_fp_cdist_p(i),(err_fpd_p(i)*cluster_fp_cdist_p(i)),&
 cluster_fp_ldist_p(i),(err_fpd_p(i)*cluster_fp_ldist_p(i)),&
 cluster_z_adist_p(i),(err_red_p(i)*cluster_z_adist_p(i)),&
 cluster_z_cdist_p(i),(err_red_p(i)*cluster_z_cdist_p(i)),&
 cluster_z_ldist_p(i),(err_red_p(i)*cluster_z_ldist_p(i)),&
n_fp_in_cluster(i),members_clustercore(i)
END IF
END DO
 CLOSE(50)

 
  
OPEN(50,file='catalogues/cluster_fp_distances_SDSS_s.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),err_z_s(i),&
 cluster_fp_adist_s(i),(err_fpd_s(i)*cluster_fp_adist_s(i)),&
 cluster_fp_cdist_s(i),(err_fpd_s(i)*cluster_fp_cdist_s(i)),&
 cluster_fp_ldist_s(i),(err_fpd_s(i)*cluster_fp_ldist_s(i)),&
 cluster_z_adist_s(i),(err_red_s(i)*cluster_z_adist_s(i)),&
 cluster_z_cdist_s(i),(err_red_s(i)*cluster_z_cdist_s(i)),&
 cluster_z_ldist_s(i),(err_red_s(i)*cluster_z_ldist_s(i)),&
n_fp_in_cluster(i),members_clustercore(i)
END IF
END DO
 CLOSE(50)

 
  
OPEN(50,file='catalogues/cluster_fp_distances_SDSS_m.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),err_z_m(i),&
 cluster_fp_adist_m(i),(err_fpd_m(i)*cluster_fp_adist_m(i)),&
 cluster_fp_cdist_m(i),(err_fpd_m(i)*cluster_fp_cdist_m(i)),&
 cluster_fp_ldist_m(i),(err_fpd_m(i)*cluster_fp_ldist_m(i)),&
 cluster_z_adist_m(i),(err_red_m(i)*cluster_z_adist_m(i)),&
 cluster_z_cdist_m(i),(err_red_m(i)*cluster_z_cdist_m(i)),&
 cluster_z_ldist_m(i),(err_red_m(i)*cluster_z_ldist_m(i)),&
n_fp_in_cluster(i),members_clustercore(i)
END IF
END DO
 CLOSE(50)

 
 
 
 
 
 
 OPEN(50,file='final_catalogues/cluster_fp_distances_SDSS_p.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(50,"(I8,2F14.8,2F10.7,12F10.4,2I6)") i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),err_z_p(i),&
 cluster_fp_adist_p(i),(err_fpd_p(i)*cluster_fp_adist_p(i)),&
 cluster_fp_cdist_p(i),(err_fpd_p(i)*cluster_fp_cdist_p(i)),&
 cluster_fp_ldist_p(i),(err_fpd_p(i)*cluster_fp_ldist_p(i)),&
 cluster_z_adist_p(i),(err_red_p(i)*cluster_z_adist_p(i)),&
 cluster_z_cdist_p(i),(err_red_p(i)*cluster_z_cdist_p(i)),&
 cluster_z_ldist_p(i),(err_red_p(i)*cluster_z_ldist_p(i)),&
n_fp_in_cluster(i),members_clustercore(i)
END IF
END DO
 CLOSE(50)


  
OPEN(50,file='final_catalogues/cluster_fp_distances_SDSS_s.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(50,"(I8,2F14.8,2F10.7,12F10.4,2I6)") i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),err_z_s(i),&
 cluster_fp_adist_s(i),(err_fpd_s(i)*cluster_fp_adist_s(i)),&
 cluster_fp_cdist_s(i),(err_fpd_s(i)*cluster_fp_cdist_s(i)),&
 cluster_fp_ldist_s(i),(err_fpd_s(i)*cluster_fp_ldist_s(i)),&
 cluster_z_adist_s(i),(err_red_s(i)*cluster_z_adist_s(i)),&
 cluster_z_cdist_s(i),(err_red_s(i)*cluster_z_cdist_s(i)),&
 cluster_z_ldist_s(i),(err_red_s(i)*cluster_z_ldist_s(i)),&
n_fp_in_cluster(i),members_clustercore(i)
END IF
END DO
 CLOSE(50)

 
  
OPEN(50,file='final_catalogues/cluster_fp_distances_SDSS_m.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(50,"(I8,2F14.8,2F10.7,12F10.4,2I6)") i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),err_z_m(i),&
 cluster_fp_adist_m(i),(err_fpd_m(i)*cluster_fp_adist_m(i)),&
 cluster_fp_cdist_m(i),(err_fpd_m(i)*cluster_fp_cdist_m(i)),&
 cluster_fp_ldist_m(i),(err_fpd_m(i)*cluster_fp_ldist_m(i)),&
 cluster_z_adist_m(i),(err_red_m(i)*cluster_z_adist_m(i)),&
 cluster_z_cdist_m(i),(err_red_m(i)*cluster_z_cdist_m(i)),&
 cluster_z_ldist_m(i),(err_red_m(i)*cluster_z_ldist_m(i)),&
n_fp_in_cluster(i),members_clustercore(i)
END IF
END DO
 CLOSE(50)

 
 
 
 
 
 
 
  
DO i=1,102
bincluster(i)=0
END DO

DO i=1,100
DO ii=1,n_clusters
IF (act_cluster(ii)) THEN
IF (n_fp_in_cluster(ii)==i) THEN
bincluster(i)=bincluster(i)+1
END IF
END IF
END DO
END DO


DO ii=1,n_clusters
IF (act_cluster(ii)) THEN
IF ((n_fp_in_cluster(ii)>100).AND.(n_fp_in_cluster(ii)<1000)) THEN
bincluster(101)=bincluster(101)+1
END IF
END IF
END DO

DO ii=1,n_clusters
IF (act_cluster(ii)) THEN
IF (n_fp_in_cluster(ii)>1000) THEN
bincluster(102)=bincluster(102)+1
END IF
END IF
END DO


OPEN(50,file='catalogues/bin_earlytypes_percluster.txt')
DO i=1,102
WRITE(50,*) i,bincluster(i)
END DO
 CLOSE(50)
 

DO i=1,n_clusters
IF (act_cluster(i)) THEN
! dist_err_m(i)=(((cluster_fp_cdist_m(i)-cluster_z_cdist_m(i))/cluster_z_cdist_m(i))*100.D0)
! dist_err_s(i)=(((cluster_fp_cdist_s(i)-cluster_z_cdist_s(i))/cluster_z_cdist_s(i))*100.D0)
! dist_err_p(i)=(((cluster_fp_cdist_p(i)-cluster_z_cdist_p(i))/cluster_z_cdist_p(i))*100.D0)
dist_err_m(i)=LOG10(cluster_fp_cdist_m(i)/cluster_z_cdist_m(i))
dist_err_s(i)=LOG10(cluster_fp_cdist_s(i)/cluster_z_cdist_s(i))
dist_err_p(i)=LOG10(cluster_fp_cdist_p(i)/cluster_z_cdist_p(i))
END IF
END DO
 
OPEN(61,file='fp_distance_variations_1.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>0) THEN
WRITE(61,*) cluster_fp_cdist_m(i),cluster_z_cdist_m(i),dist_err_m(i)
END IF
END IF
END DO
 CLOSE(61)
 
 OPEN(61,file='fp_distance_variations_10.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>10) THEN
WRITE(61,*) cluster_fp_cdist_m(i),cluster_z_cdist_m(i),dist_err_m(i)
END IF
END IF
END DO
 CLOSE(61)
 
 OPEN(61,file='fp_distance_variations_2.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>1) THEN
WRITE(61,*) cluster_fp_cdist_m(i),cluster_z_cdist_m(i),dist_err_m(i)
END IF
END IF
END DO
 CLOSE(61)
 
 OPEN(61,file='fp_distance_variations_3.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>2) THEN
WRITE(61,*) cluster_fp_cdist_m(i),cluster_z_cdist_m(i),dist_err_m(i)
END IF
END IF
END DO
 CLOSE(61)
 
 OPEN(61,file='fp_distance_variations_5.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>4) THEN
WRITE(61,*) cluster_fp_cdist_m(i),cluster_z_cdist_m(i),dist_err_m(i)
END IF
END IF
END DO
 CLOSE(61)
 

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
help_x=NINT(cluster_z_cdist_m(i)/4.D0)
! help_y=NINT(dist_err_m(i)+100.D0)
help_y=NINT((dist_err_m(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_map_1.txt')
DO i=0,200
DO ii=0,200
help_x2=DBLE(i)*4.D0
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61) 

 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>1) THEN
help_x=NINT(cluster_z_cdist_m(i)/4.D0)
! help_y=NINT(dist_err_m(i)+100.D0)
help_y=NINT((dist_err_m(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_map_2.txt')
DO i=0,200
DO ii=0,200
help_x2=DBLE(i)*4.D0
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61) 
 
 
 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>2) THEN
help_x=NINT(cluster_z_cdist_m(i)/4.D0)
! help_y=NINT(dist_err_m(i)+100.D0)
help_y=NINT((dist_err_m(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_map_3.txt')
DO i=0,200
DO ii=0,200
help_x2=DBLE(i)*4.D0
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61) 
 
 
 
 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>4) THEN
help_x=NINT(cluster_z_cdist_m(i)/4.D0)
! help_y=NINT(dist_err_m(i)+100.D0)
help_y=NINT((dist_err_m(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_map_5.txt')
DO i=0,200
DO ii=0,200
help_x2=DBLE(i)*4.D0
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61) 
 
 
 
 
 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i)>9) THEN
help_x=NINT(cluster_z_cdist_m(i)/4.D0)
! help_y=NINT(dist_err_m(i)+100.D0)
help_y=NINT((dist_err_m(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_map_10.txt')
DO i=0,200
DO ii=0,200
help_x2=DBLE(i)*4.D0
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61) 
 
 
 
 
 
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
help_x=n_fp_in_cluster(i)
! help_y=NINT(dist_err_m(i)+100.D0)
help_y=NINT((dist_err_m(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<160)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_err_m.txt')
DO i=0,100
DO ii=0,160
help_x2=DBLE(i)
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61)


 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
help_x=n_fp_in_cluster(i)
help_y=NINT((dist_err_s(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<160)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_err_s.txt')
DO i=0,100
DO ii=0,160
help_x2=DBLE(i)
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61)


 DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_clusters
IF (act_cluster(i)) THEN
help_x=n_fp_in_cluster(i)
help_y=NINT((dist_err_p(i)+0.4D0)*200.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<160)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END IF
END DO

OPEN(61,file='fp_distance_err_p.txt')
DO i=0,100
DO ii=0,160
help_x2=DBLE(i)
help_y2=(DBLE(ii)/200.D0)-0.4D0
WRITE(61,*) help_x2,help_y2,LOG10(mapmap(i,ii)+1.D0)
END DO
END DO
 CLOSE(61)

! WRITE(*,*) 'noch ok 1'
 
  av_m=0.D0
  av_s=0.D0
  av_p=0.D0
  
  sigma_m=0.D0
  sigma_s=0.D0
  sigma_p=0.D0
  
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF ((n_fp_in_cluster(i)>0).AND.(n_fp_in_cluster(i)<101)) THEN

av_p(n_fp_in_cluster(i))=av_p(n_fp_in_cluster(i))+dist_err_p(i)

av_m(n_fp_in_cluster(i))=av_m(n_fp_in_cluster(i))+dist_err_m(i)

av_s(n_fp_in_cluster(i))=av_s(n_fp_in_cluster(i))+dist_err_s(i)

END IF
END IF
END DO
!   WRITE(*,*) 'noch ok 2'
  
  DO i=1,100
  IF (bincluster(i)>0.D0) THEN
av_p(i)=av_p(i)/DBLE(bincluster(i))
av_s(i)=av_s(i)/DBLE(bincluster(i))
av_m(i)=av_m(i)/DBLE(bincluster(i))
  ELSE
av_p(i)=0.D0
av_s(i)=0.D0
av_m(i)=0.D0
  END IF
  END DO
  
!   WRITE(*,*) 'noch ok 3'
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF ((n_fp_in_cluster(i)>0).AND.(n_fp_in_cluster(i)<101)) THEN

sigma_s(n_fp_in_cluster(i))=sigma_s(n_fp_in_cluster(i))+(dist_err_s(i)-av_s(n_fp_in_cluster(i)))**2

sigma_m(n_fp_in_cluster(i))=sigma_m(n_fp_in_cluster(i))+(dist_err_m(i)-av_m(n_fp_in_cluster(i)))**2

sigma_p(n_fp_in_cluster(i))=sigma_p(n_fp_in_cluster(i))+(dist_err_p(i)-av_p(n_fp_in_cluster(i)))**2

END IF
END IF
END DO 

!   WRITE(*,*) 'noch ok 4'
  DO i=1,100
  IF (bincluster(i)>1.D0) THEN
sigma_p(i)=SQRT(sigma_p(i)/DBLE(bincluster(i)-1))
sigma_s(i)=SQRT(sigma_s(i)/DBLE(bincluster(i)-1))
sigma_m(i)=SQRT(sigma_m(i)/DBLE(bincluster(i)-1))
  ELSE
sigma_p(i)=-1.D0
sigma_s(i)=-1.D0
sigma_m(i)=-1.D0
  END IF
  END DO 
 
 
!  WRITE(*,*) 'noch ok 5'
 
 
 OPEN(61,file='fp_distance_statistics_m.txt')
DO i=1,100
IF (sigma_m(i)>0.D0) THEN
WRITE(61,*) i,av_m(i),sigma_m(i),(av_m(i)+sigma_m(i)),(av_m(i)-sigma_m(i))
END IF
END DO
 CLOSE(61)

  OPEN(61,file='fp_distance_statistics_p.txt')
DO i=1,100
IF (sigma_p(i)>0.D0) THEN
WRITE(61,*) i,av_p(i),sigma_p(i),(av_p(i)+sigma_p(i)),(av_p(i)-sigma_p(i))
END IF
END DO
 CLOSE(61)
 
  OPEN(61,file='fp_distance_statistics_s.txt')
DO i=1,100
IF (sigma_s(i)>0.D0) THEN
WRITE(61,*) i,av_s(i),sigma_s(i),(av_s(i)+sigma_s(i)),(av_s(i)-sigma_s(i))
END IF
END DO
 CLOSE(61)
 
 
   OPEN(61,file='fp_etg_ratios_all.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
WRITE(61,*) i,n_fp_in_cluster(i),members_clustercore(i),(DBLE(n_fp_in_cluster(i))/DBLE(members_clustercore(i)))
END IF
END DO
 CLOSE(61)
 
 

 help0=0
 av_ratio=0.D0
    OPEN(61,file='fp_etg_ratios_rich.txt')
DO i=1,n_clusters
IF (act_cluster(i)) THEN
IF (n_fp_in_cluster(i).GE.10) THEN
WRITE(61,*) i,n_fp_in_cluster(i),members_clustercore(i),(DBLE(n_fp_in_cluster(i))/DBLE(members_clustercore(i))),&
 cluster_fp_cdist_m(i),cluster_z_cdist_m(i)
av_ratio=av_ratio+(DBLE(n_fp_in_cluster(i))/DBLE(members_clustercore(i)))
help0=help0+1
END IF
END IF
END DO
av_ratio=av_ratio/DBLE(help0)
WRITE(61,*) av_ratio
 CLOSE(61)
 

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


 
 
 
 
 