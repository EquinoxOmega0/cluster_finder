PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer(kind=8) :: io_err,n,i,ii,iii,hint,counter,old_fof_n,fof_n,n_dreieck,help_count,help_index,central_index
character(200) :: filename,intial_values,listname,fulllistname
integer(kind=8) :: n_l,n_fl,calibrationmode,n_truegroups_mock,n_truegroups_all,n2plus_truemock,n2plus_trueall
logical :: flipvar
integer(kind=8) :: n_g2_fof,n_g2_mock,n_g2_bij,n_multi_cut,n_help,n_used

logical, allocatable :: group_bij(:),samecluster(:),oldact(:)
double precision, allocatable ::  P_mock(:),P_fof(:)

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor,R_factor,basic_link
double precision :: d_red,d_ang_dist,angular_sep,delta_z,hreal,l_next,hl_z,hl_a,l_old,mabs_vollim
double precision :: avlogmass,avmass,mag_sol_r,help_lum,s_power,mag_sol_g,weight,dummy_var,dummy_var2

double precision, allocatable ::  ra(:),dec(:),z(:),mag_g(:),mag_r(:)
logical, allocatable :: active(:),fof(:),fof_recursive(:),iter_help(:)
double precision, allocatable :: ra_clustercore(:),dec_clustercore(:),z_clustercore(:),lumtot_clustercore(:)
double precision, allocatable :: sigma_clustercore(:),dummy_core(:),radius_clustercore(:),angradius_clustercore(:)
double precision, allocatable :: lumdistance_clustercore(:),mass_clustercore(:),colour_clustercore(:),mstar_cluster(:)
double precision, allocatable :: lumotherband(:),lum_g(:),sigma_gap(:),lumobs_clustercore(:),mdyn_clustercore(:)
integer(kind=8), allocatable :: members_clustercore(:),clusterindex(:),clustergalnmember(:),n_stellar(:)
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

double precision :: sigma_lumtot_s,sigma_lumtot_ml,sigma_lumtot_mh


double precision :: E_fof,E_mock,E_tot,Q_fof,Q_mock,Q_tot,S_tot,P_help,rms_s,rms_ml,rms_mh

double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y
integer, dimension(1:102) :: binclustermock,binclusterfof

character(20) :: appendix


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme RESIDUAL PLOT started'
WRITE(*,*) '============================================================'



! define constants
PI=ACOS(-1.D0)
appendix='SDSS'


 appendix=TRIM(adjustl(appendix))

 filename='final_catalogues/cluster_list_'//TRIM(appendix)//'.txt'
 




! get length of file
OPEN(50,file=filename)
io_err=0
n=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n=n+1
END DO
 CLOSE(50)
n=n-1

WRITE(*,*) n,'galaxies will be used'





allocate(active(1:n))
allocate(oldact(1:n))
allocate(volume_weight(1:n))
allocate(vol_int(1:n))
allocate(lum_dist(1:n))
allocate(b_eff(1:n))
allocate(R_eff(1:n))
! allocate(mass_stretch(1:n))
allocate(angular_dist(1:n))
allocate(basic_link_a(1:n))
allocate(basic_link_R(1:n))
allocate(ra_rad(1:n))
allocate(dec_rad(1:n))
allocate(lum(1:n))
allocate(ra_clustercore(1:n))
allocate(dec_clustercore(1:n))
allocate(z_clustercore(1:n))
allocate(lumtot_clustercore(1:n))
allocate(sigma_clustercore(1:n))
allocate(members_clustercore(1:n))
allocate(radius_clustercore(1:n))
allocate(angradius_clustercore(1:n))
allocate(lumdistance_clustercore(1:n))
allocate(mass_clustercore(1:n))
allocate(lumobs_clustercore(1:n))
allocate(mdyn_clustercore(1:n))
allocate(fof(1:n))
allocate(fof_recursive(1:n))
allocate(iter_help(1:n))
allocate(dummy_core(1:n))
allocate(clusterindex(1:n))
allocate(colour_clustercore(1:n))
allocate(lumotherband(1:n))
allocate(lum_g(1:n))
allocate(vred(1:n))
allocate(sigma_gap(1:n))
allocate(lum_gal(1:n))
allocate(lum_int(1:n))
allocate(mass_err(1:n))
allocate(lumtot_err(1:n))
allocate(mstar_cluster(1:n))
allocate(n_stellar(1:n))




OPEN(50,file=filename)
DO i=1,n
!READ(50,*) ii,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
!LOG10(lumtot_clustercore(i)),lumtot_err(i),LOG10(lumobs_clustercore(i)),&
!mass_clustercore(i),mass_err(i),LOG10(mdyn_clustercore(i)),&
!sigma_clustercore(i),radius_clustercore(i),(angradius_clustercore(i)/PI*180.D0),lumdistance_clustercore(i),&
!members_clustercore(i)

READ(50,*) ii,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
lumtot_clustercore(i),lumtot_err(i),lumobs_clustercore(i),&
mass_clustercore(i),mass_err(i),mstar_cluster(i),n_stellar(i),mdyn_clustercore(i),&
sigma_clustercore(i),radius_clustercore(i),angradius_clustercore(i),lumdistance_clustercore(i),&
members_clustercore(i)
END DO
 CLOSE(50)

 
 

OPEN(50,file='distdep_catplot_all_'//TRIM(appendix)//'.txt')
DO i=1,n
WRITE(50,*) lumdistance_clustercore(i),lumtot_clustercore(i),mass_clustercore(i)
END DO
 CLOSE(50)
 
OPEN(50,file='distdep_catplot_multi_'//TRIM(appendix)//'.txt')
DO i=1,n
IF (members_clustercore(i)>1) THEN
WRITE(50,*) lumdistance_clustercore(i),mdyn_clustercore(i),&
LOG10(sigma_clustercore(i)),LOG10(radius_clustercore(i)),members_clustercore(i)
END IF
END DO
 CLOSE(50)
 

 
 
 
 
 
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((lumtot_clustercore(i)-7.D0)*20.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_lum_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,200
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/20.D0)+7.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

  
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((lumobs_clustercore(i)-7.D0)*20.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_lumobs_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,200
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/20.D0)+7.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

   
   
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((mstar_cluster(i)-6.D0)*20.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_mstar_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,200
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/20.D0)+7.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 


 
 
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((mass_clustercore(i)-10.D0)*25.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_mass_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/25.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 
  
 
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
IF (members_clustercore(i)>1) THEN
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((mdyn_clustercore(i)-7.D0)*20.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_mdyn_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/20.D0)+7.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
   
 
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
IF (members_clustercore(i)>1) THEN
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((LOG10(sigma_clustercore(i))-1.4D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_sigma_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/50.D0)+1.4D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
   
 
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
IF (members_clustercore(i)>1) THEN
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((LOG10(radius_clustercore(i)))*40.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_rad_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/40.D0)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
    
 
DO i=0,1000
DO ii=0,1000
mapmap(i,ii)=0
END DO
END DO

DO i=1,n
help_x=NINT((lumdistance_clustercore(i))/2.D0)
help_y=NINT((LOG10(DBLE(members_clustercore(i)))+0.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<200)) THEN
IF (mapmap(help_x,help_y)<1.5D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+(2.D0/mapmap(help_x,help_y))
END IF
END IF
END IF
END DO

OPEN(61,file='distdep_nvis_'//TRIM(appendix)//'.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)*2.D0)
help_y2=(DBLE(ii)/50.D0)-0.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 
 
 
 
 

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


 
 
 
 
 