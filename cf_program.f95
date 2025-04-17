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
double precision, allocatable :: lumdistance_clustercore(:),mass_clustercore(:),colour_clustercore(:)
double precision, allocatable :: lumotherband(:),lum_g(:),sigma_gap(:),lumobs_clustercore(:),mdyn_clustercore(:)
integer(kind=8), allocatable :: members_clustercore(:),clusterindex(:),clustergalnmember(:)
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
WRITE(*,*) '    programme CLUSTER FINDER started'
WRITE(*,*) '============================================================'



! define constants
PI=ACOS(-1.D0)

OPEN(50,file='input_cf.txt')

READ(50,*) Omega_m,Omega_l,H0
READ(50,*) b_factor,R_factor,lum_power
READ(50,*) appendix
 CLOSE(50)
 appendix=TRIM(adjustl(appendix))
 n_multi_cut=n_multi_cut-1

 
 filename=TRIM(appendix)//'_real/'//TRIM(appendix)//'_obs_galaxies_radecz.txt'
 
 
!  OPEN(50,file='fit_lumass_'//TRIM(appendix)//'.txt')
! READ(50,*) mfit_a,mfit_b,mfit_c
! READ(50,*)
! READ(50,*)
!  CLOSE(50)
 
  OPEN(50,file='globalparameters_'//TRIM(appendix)//'.txt')
READ(50,*) basic_link
READ(50,*)
READ(50,*) av_rad_pec_vel
READ(50,*) avlogmass
 CLOSE(50)

! 
  OPEN(61,file='mass/mfit_'//TRIM(appendix)//'_s.txt')
READ(61,*) rms_s
DO i=1,7
READ(61,*) mass_coeff_s(i),dummy_var
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/mfit_'//TRIM(appendix)//'_ml.txt')
READ(61,*) rms_ml
DO i=1,9
READ(61,*) mass_coeff_ml(i),dummy_var
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/mfit_'//TRIM(appendix)//'_mh.txt')
READ(61,*) rms_mh
DO i=1,10
READ(61,*) mass_coeff_mh(i),dummy_var
END DO
 CLOSE(61)
 
 OPEN(61,file='mass/lumdev_'//TRIM(appendix)//'.txt')
 READ(61,*) sigma_lumtot_s,sigma_lumtot_ml,sigma_lumtot_mh 
 CLOSE(61)
  
  
  

 
IF (TRIM(appendix)=='SDSS') THEN
area_cover=9376.D0/41253.D0
mag_sol_r=4.71D0
mag_sol_g=5.31D0     
red_measure=30.D0
mag_limit=17.77D0
magmin=-15.D0
magmax=-30.D0
ELSE
area_cover=0.91D0
mag_sol_r=3.28D0 !Ks
mag_sol_g=3.64D0 !J
red_measure=32.D0 !31.96D0
mag_limit=11.75D0
magmin=-18.D0
magmax=-30.D0
END IF




q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)

 cV_mill=(380.D6/h)**3
m_expected=cV_mill*rho_crit*Omega_m
part_exp=(2160.D0**3)*((330.D0/500.D0)**3)



mock_cover=1.D0/8.D0




appendix=TRIM(adjustl(appendix))

deltamag=(-magmax+magmin)/1500.D0




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



allocate(ra(1:n))
allocate(dec(1:n))
allocate(z(1:n))
allocate(mag_g(1:n))
allocate(mag_r(1:n))
allocate(id1(1:n))
allocate(id2(1:n))
allocate(numberofgal(1:n))

! read file
OPEN(50,file=filename)
DO i=1,n
READ(50,*) ra(i),dec(i),z(i),mag_g(i),mag_r(i)
END DO
CLOSE(50)



   OPEN(50,file=TRIM(appendix)//'_real/'//TRIM(appendix)//'_id_file.txt')
 DO i=1,n  
 numberofgal(i)=i
 IF (TRIM(appendix)=='SDSS') THEN
   READ(50,*) id1(i) 
ELSE
   READ(50,*) dummy_var,dummy_var2
   id1(i)=NINT(dummy_var)
   id2(i)=NINT(dummy_var2)
END IF
   END DO
 CLOSE(50)




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


DO i=1,n
 clusterindex(i)=0
ra_rad(i)=ra(i)*PI/180.D0
dec_rad(i)=dec(i)*PI/180.D0
END DO

WRITE(*,*) 'initialisation done'

! calculate volume weights to corret for Malmquist bias
V_sum=0.D0
DO i=1,n

active(i)=.TRUE.
volume_weight(i)=0.D0

IF (mag_r(i)<magmax) THEN
active(i)=.FALSE.
END IF

IF (mag_r(i)>magmin) THEN
active(i)=.FALSE.
END IF

IF (mag_g(i)<(magmax-1.5D0)) THEN
active(i)=.FALSE.
END IF

IF (mag_g(i)>(magmin+1.5D0)) THEN
active(i)=.FALSE.
END IF


IF ((mag_g(i)-mag_r(i))>2.D0) THEN
active(i)=.FALSE.
END IF

IF ((mag_g(i)-mag_r(i))<-1.D0) THEN
active(i)=.FALSE.
END IF

END DO




!consider size difference of halos depending on the luminosity of the central galaxy

DO i=1,n
IF (active(i)) THEN
lum(i)=10.D0**(-0.4D0*(mag_r(i)-mag_sol_r)) 
lum_g(i)=10.D0**(-0.4D0*(mag_g(i)-mag_sol_g)) 
! mass_halo_fit=mfit_a*(LOG10(lum(i)))**2+mfit_b*(LOG10(lum(i)))+mfit_c
vred(i)=((((1.D0+z(i))**2)-1.D0)/(((1.D0+z(i))**2)+1.D0))*light
avmass=10.D0**avlogmass
! mass_halo_fit=10.D0**mass_halo_fit
lum_gal(i)=lum(i)
! mass_stretch(i)=mass_halo_fit/avmass
END IF
END DO
! WRITE(*,*) 'rescaling linking length according to the masses'



IF (TRIM(appendix)=='SDSS') THEN
DO i=1,n
IF (active(i)) THEN
DO ii=1,i-1
IF (ii.NE.i) THEN

IF (id1(i)==id1(ii)) THEN
active(i)=.FALSE.
END IF

END IF
END DO
END IF
END DO
END IF


DO i=1,n

IF (active(i)) THEN
D_limit=-0.2D0*mag_r(i)+((mag_limit+5.D0)/5.D0)
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


volume_weight(i)=1.D0/V_limit

V_sum=V_sum+volume_weight(i)

lum_gal(i)=lum_gal(i)*volume_weight(i)
lum_sum=lum_sum+lum_gal(i)


END IF

END DO 




DO i=1,n
IF (active(i)) THEN
lum_gal(i)=lum_gal(i)/lum_sum
END IF
END DO

! 
DO i=1,n
IF (active(i)) THEN
divisor=(SQRT(1.D0+2.D0*q0*z(i))+1.D0+q0*z(i))
lum_dist(i)=light/H0*z(i)*(1.D0+((z(i)*(1.D0-q0))/divisor))
lum_dist(i)=lum_dist(i)*1.D6
lum_dist(i)=LOG10(lum_dist(i))
END IF
END DO





! rescale linking length according to the completeness of the luminosity function 


DO i=1,n 
IF (active(i)) THEN
vol_int(i)=1.D-15 !about zero, but to avoid division by zero later we make it a little larger
mabs_vollim=mag_limit+5.D0-5.D0*lum_dist(i)
DO ii=1,n
IF (active(ii)) THEN
IF (mag_r(ii)<mabs_vollim) THEN
vol_int(i)=vol_int(i)+volume_weight(ii)
END IF
END IF
END DO


END IF
END DO




WRITE(*,*) 'incompletness effects corrected'



! adopt linking lenght for redshift space distortions and trigonemtry due to distance
DO i=1,n
IF (active(i)) THEN
angular_dist(i)=(10.D0**lum_dist(i))*((1.D0+z(i))**(-2.D0))/1.D6
basic_link_a(i)=ATAN(basic_link/angular_dist(i))!*180.D0/PI
IF (basic_link_a(i)<0.D0) THEN
basic_link_a(i)=-basic_link_a(i)
END IF
basic_link_R(i)=(basic_link*H0+av_rad_pec_vel*2.D0)/light
END IF
END DO
WRITE(*,*) 'redshift distortions and projection effects considered'




! modify linking lenght according to the modification calculated before
DO i=1,n
IF (active(i)) THEN
b_eff(i)=b_factor*basic_link_a(i)*((vol_int(i)/V_sum)**(-lum_power/3.D0))
R_eff(i)=R_factor*basic_link_R(i)*((vol_int(i)/V_sum)**(-lum_power/3.D0))

END IF
END DO

WRITE(*,*) 'linking length adapted'




OPEN(50,file='linking_conditions.txt')
DO i=1,n
IF (active(i)) THEN
WRITE(50,*) i,b_eff(i),R_eff(i),basic_link_a(i),basic_link_R(i),z(i),&
((vol_int(i)/V_sum)**(-lum_power/3.D0))
END IF
END DO
 CLOSE(50)

 n_used=0
DO i=1,n 
IF (active(i)) THEN
oldact(i)=.TRUE.
n_used=n_used+1
ELSE
oldact(i)=.FALSE.
END IF
END DO
WRITE(*,*) 'using',n_used,'of',n,'galaxies'
 
! run FOF algoritm

 counter=0

WRITE(*,*) 'start FOF algoritm'
DO i=1,n
WRITE(*,*) 'doing ',i,' of ',n

IF (active(i)) THEN

active(i)=.FALSE.
! WRITE(*,*) i

DO ii=1,n
fof(ii)=.FALSE.
fof_recursive(ii)=.FALSE.
END DO

fof(i)=.TRUE.
fof_recursive(i)=.TRUE.
fof_n=1
old_fof_n=0

DO WHILE (fof_n>old_fof_n)

old_fof_n=fof_n

DO iii=1,n

IF (fof_recursive(iii)) THEN

fof_recursive(iii)=.FALSE.

DO ii=1,n 
IF (active(ii)) THEN

angular_sep=COS(dec_rad(iii))*COS(dec_rad(ii))*COS((ra_rad(ii)-ra_rad(iii)))
angular_sep=angular_sep+SIN(dec_rad(iii))*SIN(dec_rad(ii))
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

delta_z=z(iii)-z(ii)

IF (delta_z<0.D0) THEN
delta_z=-delta_z
END IF
 
! av_z=(z(iii)+z(ii))/2.D0

 

IF ((delta_z<R_eff(iii)).OR.(delta_z<R_eff(ii))) THEN

IF ((angular_sep<b_eff(iii)).OR.(angular_sep<b_eff(ii))) THEN

IF (fof(ii).EQV..FALSE.) THEN
fof_recursive(ii)=.TRUE.
END IF

fof(ii)=.TRUE.
active(ii)=.FALSE.

END IF
END IF


END IF
END DO

END IF
END DO

fof_n=0
DO ii=1,n 
IF (fof(ii)) THEN
fof_n=fof_n+1
END IF
END DO

END DO







!CLUSTER berechnen und hinzufügen
 counter=counter+1


ra_clustercore(counter)=0.D0
dec_clustercore(counter)=0.D0
z_clustercore(counter)=0.D0
lumtot_clustercore(counter)=0.D0
lumotherband(counter)=0.D0
sigma_clustercore(counter)=0.D0
angradius_clustercore(counter)=0.D0
lumdistance_clustercore(counter)=0.D0
mass_clustercore(counter)=0.D0
members_clustercore(counter)=0
mdyn_clustercore(counter)=0.D0
 colour_clustercore(i)=0.D0

! using iter method from Robotham to find cluster centre


help_index=0

DO ii=1,n
IF (fof(ii)) THEN
members_clustercore(counter)=members_clustercore(counter)+1
 clusterindex(ii)=counter
lumtot_clustercore(counter)=lumtot_clustercore(counter)+(lum(ii))
lumotherband(counter)=lumotherband(counter)+(lum_g(ii))
help_index=ii
END IF
END DO


! WRITE(*,*) members_clustercore(counter),'members in cluster ',counter

! DO ii=1,n
! IF (fof(ii)) THEN
! WRITE(*,*) ii
! END IF
! END DO


IF (members_clustercore(counter)>1) THEN


DO ii=1,n
iter_help(ii)=fof(ii)
END DO




help_count=10

DO WHILE (help_count>1) 
ra_clustercore(counter)=0.D0
dec_clustercore(counter)=0.D0
help_lum=0.D0
DO ii=1,n
IF (iter_help(ii)) THEN
ra_clustercore(counter)=ra_clustercore(counter)+(ra_rad(ii)*lum(ii))
dec_clustercore(counter)=dec_clustercore(counter)+(dec_rad(ii)*lum(ii))
help_lum=help_lum+lum(ii)
END IF
END DO
ra_clustercore(counter)=ra_clustercore(counter)/help_lum
dec_clustercore(counter)=dec_clustercore(counter)/help_lum



DO ii=1,n
IF (iter_help(ii)) THEN
angular_sep=COS(dec_clustercore(counter))*COS(dec_rad(ii))*COS((ra_rad(ii)-ra_clustercore(counter)))
angular_sep=angular_sep+SIN(dec_clustercore(counter))*SIN(dec_rad(ii))
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
DO ii=1,n
IF (iter_help(ii)) THEN
IF (angular_sep<dummy_core(ii)) THEN
angular_sep=dummy_core(ii)
help_index=ii
END IF
END IF
END DO
iter_help(help_index)=.FALSE.

help_count=0
DO ii=1,n
IF (iter_help(ii)) THEN
help_count=help_count+1
END IF
END DO


END DO

DO ii=1,n
IF (iter_help(ii)) THEN
ra_clustercore(counter)=ra_rad(ii)
dec_clustercore(counter)=dec_rad(ii)
 central_index=ii
END IF
END DO



! use median redshift for centre
help_count=0
DO ii=1,n
IF (fof(ii)) THEN
help_count=help_count+1
dummy_core(help_count)=z(ii)
END IF
END DO


DO ii=1,help_count
DO iii=ii,help_count

IF (dummy_core(ii)<dummy_core(iii)) THEN
help_lum=dummy_core(ii)
dummy_core(ii)=dummy_core(iii)
dummy_core(iii)=help_lum
END IF

END DO
END DO


IF (MOD(help_count,2)==0) THEN
help_index=help_count/2
z_clustercore(counter)=(dummy_core(help_index)+dummy_core(help_index+1))/2.D0
ELSE
help_index=help_count/2
z_clustercore(counter)=dummy_core(help_index+1)
END IF



help_count=0
DO ii=1,n
IF (fof(ii)) THEN
help_count=help_count+1

angular_sep=COS(dec_clustercore(counter))*COS(dec_rad(ii))*COS((ra_rad(ii)-ra_clustercore(counter)))
angular_sep=angular_sep+SIN(dec_clustercore(counter))*SIN(dec_rad(ii))
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

dummy_core(help_count)=angular_sep
END IF
END DO


DO ii=1,help_count
DO iii=ii,help_count

IF (dummy_core(ii)<dummy_core(iii)) THEN
help_lum=dummy_core(ii)
dummy_core(ii)=dummy_core(iii)
dummy_core(iii)=help_lum
END IF

END DO
END DO


IF (MOD(help_count,2)==0) THEN
help_index=help_count/2
angradius_clustercore(counter)=(dummy_core(help_index)+dummy_core(help_index+1))/2.D0
ELSE
help_index=help_count/2
angradius_clustercore(counter)=dummy_core(help_index+1)
END IF





n_help=help_count

help_count=0
DO ii=1,n
IF (fof(ii)) THEN
help_count=help_count+1
dummy_core(help_count)=vred(ii)
END IF
END DO


DO ii=1,n_help
DO iii=ii,n_help
IF (dummy_core(ii)<dummy_core(iii)) THEN
dummy_var=dummy_core(ii)
dummy_core(ii)=dummy_core(iii)
dummy_core(iii)=dummy_var
END IF
END DO
END DO

sigma_gap(counter)=0.D0

DO ii=1,(n_help-1)
weight=DBLE(ii*(n_help-ii))
sigma_gap(counter)=sigma_gap(counter)+weight*ABS(dummy_core(ii+1)-dummy_core(ii))
END DO

sigma_gap(counter)=sigma_gap(counter)*SQRT(PI)/DBLE(n_help*(n_help-1))

IF (sigma_gap(counter)<red_measure) THEN
sigma_gap(counter)=red_measure
END IF

sigma_clustercore(counter)=SQRT(DBLE(n_help)/DBLE(n_help-1)*(sigma_gap(counter)**2)-(red_measure**2))


IF (sigma_clustercore(counter)<red_measure) THEN
sigma_clustercore(counter)=red_measure
END IF





ra_clustercore(counter)=ra(central_index)
dec_clustercore(counter)=dec(central_index)



ELSE 


ra_clustercore(counter)=ra(help_index)
dec_clustercore(counter)=dec(help_index)
z_clustercore(counter)=z(help_index)
END IF


END IF
END DO



WRITE(*,*) counter,'clusters found'

DO i=1,counter
!luminosity distance in Mpc!
divisor=(SQRT(1.D0+2.D0*q0*z_clustercore(i))+1.D0+q0*z_clustercore(i))
lumdistance_clustercore(i)=light/H0*z_clustercore(i)*(1.D0+((z_clustercore(i)*(1.D0-q0))/divisor))

END DO



DO i=1,counter
lum_int(i)=0.D0
mabs_vollim=mag_limit+5.D0-5.D0*LOG10(lumdistance_clustercore(i)*1.D6)
DO ii=1,n
IF (oldact(ii)) THEN
IF (mag_r(ii)<mabs_vollim) THEN
lum_int(i)=lum_int(i)+lum_gal(ii)
END IF
END IF
END DO
! WRITE (*,*) lum_int(i)
lumobs_clustercore(i)=lumtot_clustercore(i)
lumtot_clustercore(i)=lumtot_clustercore(i)/lum_int(i)
END DO



DO i=1,counter
angdist_help=lumdistance_clustercore(i)*((1.D0+z_clustercore(i))**(-2.D0))*1000.D0  !ang. dist in kpc
mdyn_clustercore(i)=1.D0
loglum=LOG10(lumtot_clustercore(i))
loglumdist=LOG10(lumdistance_clustercore(i))
 
 IF (members_clustercore(i)>1) THEN
 radius_clustercore(i)=TAN(angradius_clustercore(i))*angdist_help
 mdyn_clustercore(i)=10.D0*(3.D0*(sigma_clustercore(i)**2))*radius_clustercore(i)/G
 logmdyn=LOG10(mdyn_clustercore(i))
 logmember=LOG10(DBLE(members_clustercore(i)))
 lograd=LOG10(radius_clustercore(i))
 logsigma=LOG10(sigma_clustercore(i))
 END IF
 
mass_clustercore(i)=0.D0
 ! actually log10(mass)
IF (members_clustercore(i)==1) THEN
! WRITE(*,*) loglum,loglumdist
mass_clustercore(i)=mass_coeff_s(1)*loglum+mass_coeff_s(2)*(loglum**2)+mass_coeff_s(3)*(loglum**3)
! WRITE(*,*) mass_clustercore(i)
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_s(4)*loglumdist+mass_coeff_s(5)*(loglumdist**2)
! WRITE(*,*) mass_clustercore(i)
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_s(6)*(loglumdist**3)+mass_coeff_s(7)
! WRITE(*,*) mass_coeff_s(1:7)

mass_err(i)=rms_s
lumtot_err(i)=sigma_lumtot_s

ELSE 
IF (members_clustercore(i)<5) THEN
mass_clustercore(i)=mass_coeff_ml(1)*loglum+mass_coeff_ml(2)*(loglum**2)+mass_coeff_ml(3)*(loglum**3)
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_ml(4)*loglumdist+mass_coeff_ml(5)*(loglumdist**2)
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_ml(6)*(loglumdist**3)+mass_coeff_ml(7)*logsigma
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_ml(8)*lograd+mass_coeff_ml(9)

mass_err(i)=rms_ml
lumtot_err(i)=sigma_lumtot_ml
ELSE
mass_clustercore(i)=mass_coeff_mh(1)*loglum+mass_coeff_mh(2)*(loglum**2)+mass_coeff_mh(3)*(loglum**3)
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_mh(4)*loglumdist+mass_coeff_mh(5)*(loglumdist**2)
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_mh(6)*(loglumdist**3)+mass_coeff_mh(7)*logsigma
mass_clustercore(i)=mass_clustercore(i)+mass_coeff_mh(8)*lograd+mass_coeff_mh(9)*logmember+mass_coeff_mh(10)

mass_err(i)=rms_mh
lumtot_err(i)=sigma_lumtot_mh

END IF
END IF



END DO





OPEN(50,file='catalogues/cluster_list_all_'//TRIM(appendix)//'.txt')
DO i=1,counter
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
LOG10(lumtot_clustercore(i)),lumtot_err(i),LOG10(lumobs_clustercore(i)),&
mass_clustercore(i),mass_err(i),LOG10(mdyn_clustercore(i)),&
sigma_clustercore(i),radius_clustercore(i),(angradius_clustercore(i)/PI*180.D0),lumdistance_clustercore(i),&
members_clustercore(i)
END DO
 CLOSE(50)

OPEN(50,file='catalogues/cluster_list_single_'//TRIM(appendix)//'.txt')
DO i=1,counter
IF (members_clustercore(i)==1) THEN
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
LOG10(lumtot_clustercore(i)),lumtot_err(i),LOG10(lumobs_clustercore(i)),&
mass_clustercore(i),mass_err(i),lumdistance_clustercore(i)
END IF
END DO
 CLOSE(50)

OPEN(50,file='catalogues/cluster_list_multi_'//TRIM(appendix)//'.txt')
DO i=1,counter
IF (members_clustercore(i)>1) THEN
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
LOG10(lumtot_clustercore(i)),lumtot_err(i),LOG10(lumobs_clustercore(i)),&
mass_clustercore(i),mass_err(i),LOG10(mdyn_clustercore(i)),&
sigma_clustercore(i),radius_clustercore(i),(angradius_clustercore(i)/PI*180.D0),lumdistance_clustercore(i),&
members_clustercore(i)
END IF
END DO
 CLOSE(50)
 
 OPEN(50,file='catalogues/cluster_list_massive_'//TRIM(appendix)//'.txt')
DO i=1,counter
IF ((mass_clustercore(i)>15.D0).OR.(members_clustercore(i)>100)) THEN
WRITE(50,*) i,ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
LOG10(lumtot_clustercore(i)),lumtot_err(i),LOG10(lumobs_clustercore(i)),&
mass_clustercore(i),mass_err(i),LOG10(mdyn_clustercore(i)),&
sigma_clustercore(i),radius_clustercore(i),(angradius_clustercore(i)/PI*180.D0),lumdistance_clustercore(i),&
members_clustercore(i)
END IF
END DO
 CLOSE(50)
 
 
 
OPEN(50,file='catalogues/cluster_members_'//TRIM(appendix)//'.txt')
DO i=1,n
IF (oldact(i)) THEN
 IF (TRIM(appendix)=='SDSS') THEN
   WRITE(50,*) numberofgal(i),clusterindex(i),id1(i) 
ELSE
   WRITE(50,*)  numberofgal(i),clusterindex(i),id1(i),id2(i)
END IF
END IF
END DO
 CLOSE(50)
 
 
 
 
OPEN(50,file='catalogues/generaldata_'//TRIM(appendix)//'.txt')
WRITE(50,*) 'number of clusters:',counter
WRITE(50,*) 'with',n_used,'galaxies'
 CLOSE(50)
 

 
 
DO i=1,102
binclusterfof(i)=0
END DO

DO i=1,100
DO ii=1,counter
IF (members_clustercore(ii)==i) THEN
binclusterfof(i)=binclusterfof(i)+1
END IF
END DO
END DO


DO ii=1,counter
IF ((members_clustercore(ii)>100).AND.(members_clustercore(ii)<1000)) THEN
binclusterfof(101)=binclusterfof(101)+1
END IF
END DO

DO ii=1,counter
IF (members_clustercore(ii)>1000) THEN
binclusterfof(102)=binclusterfof(102)+1
END IF
END DO


OPEN(50,file='catalogues/binclusterfof_'//TRIM(appendix)//'.txt')
DO i=1,102
WRITE(50,*) i,binclusterfof(i)
END DO
 CLOSE(50)
 

 

 
 

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


 
 
 
 
 
