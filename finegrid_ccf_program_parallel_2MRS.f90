PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
include 'mpif.h'

! double precision :: simplex_alpha,simplex_beta,simplex_gamma,simplex_delta
! double precision, dimension(1:3,1:4) :: x_simplex
double precision, dimension(1:3) :: ax_c,x_oldlow,x_h,x_l,x_c,x_alpha,x_gamma,x_help,x_beta,err_X
double precision, dimension(0:10,0:10) :: S_function
integer e,ee,eee,eeee
 
integer :: itercount,highest,second,lowest
logical :: continueloop,notskip
double precision :: S_h,S_s,S_l,com_change,S_alpha,S_gamma,S_beta,S_old_low,s_help
double precision ::  b_factor_opt,R_factor_opt,lum_power_opt,S_tot_opt
double precision :: b_used,R_used,s_used,lum_used

integer, dimension(1:8) :: u_array
integer, dimension(1:1) :: u_local
double precision, dimension(1:8) :: S_array
double precision, dimension(1:1) :: S_local

character(20) :: appendix
double precision :: calculategroupcostfunction
double precision :: uuu,S_tot_av,dummy_var,errorsum
!parallization variables		
integer ierr, rankrank, nbnodes, namelen
 character (len=MPI_MAX_PROCESSOR_NAME) :: name

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme CLUSTER FINDER started'
WRITE(*,*) '============================================================'



call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rankrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nbnodes, ierr)
call MPI_GET_PROCESSOR_NAME(name, namelen, ierr)



! OPEN(50,file='cal_input.txt')
! READ(50,*) dummy_var
! READ(50,*) appendix
! READ(50,*) dummy_var,dummy_var
! 
!  CLOSE(50)
 appendix='2MRS'
 
appendix=TRIM(adjustl(appendix))


!  simplex_alpha=1.0
!  simplex_beta=0.35
!  simplex_gamma=2.0
!  simplex_delta=0.5

!  b_factor_opt=0.4
!  R_factor_opt=0.6
!  s_power_opt=0.1
!  lum_power_opt=1.0


 

DO e=1,8
u_array(e)=e
END DO


! s_used=0.023D0
! lum_used=1.D0
lum_used=0.58289910917027388

DO e=0,10
DO ee=0,10

b_used=(DBLE(e)/100.D0)+0.55D0
R_used=(DBLE(ee)/100.D0)+0.65D0

x_help(1)=b_used
x_help(2)=R_used
! x_help(3)=s_used
x_help(3)=lum_used
 
WRITE(*,*) ((e-1)*10+(ee-1))+1

call MPI_Scatter(u_array,1,MPI_INTEGER,u_local,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

uuu=DBLE(u_local(1))
S_local(1)=calculategroupcostfunction(x_help,uuu)
  
call MPI_AllGather(S_local,1,MPI_DOUBLE_PRECISION,S_array,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)



S_tot_av=0.D0
DO eee=1,8
IF (S_array(eee)<0.D0) THEN
S_array(eee)=0.D0
END IF
END DO

DO eee=1,8
DO eeee=1,8
IF (S_array(eee)>S_array(eeee)) THEN
s_help=S_array(eee)
S_array(eee)=S_array(eeee)
S_array(eeee)=s_help
END IF
END DO
END DO


! S_tot_av=(S_array(4)+S_array(5))/2.D0
  S_tot_av=0.D0
 DO eee=1,8
 S_tot_av=S_tot_av+S_array(eee)
 END DO
 S_tot_av=S_tot_av/8.D0
 
 S_function(e,ee)=S_tot_av
 
 END DO
 END DO
 
 OPEN(50,file='grid_fine_2MRS.txt')
DO e=0,10
DO ee=0,10
b_used=(DBLE(e)/100.D0)+0.55D0
R_used=(DBLE(ee)/100.D0)+0.65D0

WRITE(50,*) b_used,R_used,S_function(e,ee)

END DO
END DO
 CLOSE(50)


  
!finalize parallization
call MPI_FINALIZE(ierr)

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'


END PROGRAM


 
 
 
 
 
 
 
 
 
 
 


! function for Gaussian error
double precision FUNCTION calculategroupcostfunction(x_pos,u_trans)
IMPLICIT NONE


double precision :: PI,H0,q0,light,G,tophat,rho_crit,area_cover,dummy_var,P_help
integer(kind=8) :: io_err,n,i,ii,iii,hint,counter,old_fof_n,fof_n,n_dreieck,help_count,help_index
character(200) :: filename,intial_values,listname,fulllistname
integer(kind=8) :: n_l,n_fl,n_truegroups_mock,n_truegroups_all,n2plus_truemock,n2plus_trueall,n_used
logical :: flipvar
integer(kind=8) :: n_g2_fof,n_g2_mock,n_g2_bij,n_multi_cut,n_help,n_multi_cut_higher,central_index

integer :: uuu

logical, allocatable :: group_bij(:),samecluster(:)
double precision, allocatable ::  P_mock(:),P_fof(:)

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor,R_factor,basic_link
double precision :: d_red,d_ang_dist,angular_sep,delta_z,hreal,l_next,hl_z,hl_a,l_old,mabs_vollim
double precision ::mass_halo_fit,avlogmass,avmass,mag_sol_r,help_lum,s_power,lum_sum

double precision, allocatable ::  ra(:),dec(:),z(:),mag_g(:),mag_r(:)
! integer, allocatable :: neighbours(:)
logical, allocatable :: active(:),fof(:),fof_recursive(:),iter_help(:),mockfofid_act(:),oldact(:)
integer(kind=8), allocatable :: mmockid(:),mockid(:),mockfofid(:)

double precision, allocatable :: ra_clustercore(:),dec_clustercore(:),z_clustercore(:),lumtot_clustercore(:)
double precision, allocatable :: dummy_core(:),vred(:)
integer(kind=8), allocatable :: members_clustercore(:),clusterindex(:),clustergalnmember(:)
integer(kind=8), allocatable :: trueclusterlist(:),trueclustermembers(:),trueclusterfofid(:),truegalnmember(:)
integer(kind=8), allocatable :: clusterlist_ntrue(:),truelist_ncluster(:)

double precision, allocatable ::  angular_dist(:),ra_rad(:),dec_rad(:),lum(:),volume_weight(:),lum_gal(:)
double precision, allocatable ::  lum_dist(:),vol_int(:),b_eff(:),R_eff(:)!,mass_stretch(:)
double precision, allocatable ::  basic_link_a(:),basic_link_R(:)
double precision :: mag_limit,D_limit,t1,t2,t3,z_limit,V_limit,V_sum,binhelp,mag_vislim,magmin,magmax
double precision :: deltamag,strechdistant,av_rad_pec_vel,mfit_a,mfit_b,mfit_c
double precision :: sdss_cover,twomass_cover,mock_cover,lum_power

double precision :: E_fof,E_mock,E_tot,Q_fof,Q_mock,Q_tot,S_tot,P_helpm,u_trans

double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y
integer, dimension(1:102) :: binclustermock,binclusterfof

character(20) :: apx,nm,ecs


integer :: u,exec_counter

double precision, dimension(1:3) :: x_pos



!parallization procedures
! uuu=1
b_factor=x_pos(1)
R_factor=x_pos(2)
! s_power=x_pos(3)
lum_power=x_pos(3)

uuu=NINT(u_trans)


! define constants
PI=ACOS(-1.D0)



OPEN(50,file='cal_input_2MRS.txt')
READ(50,*) mag_limit
READ(50,*) apx
READ(50,*) magmin,magmax
READ(50,*) n_multi_cut
 CLOSE(50)
 
 Omega_m=0.25D0
 Omega_l=0.75D0
 H0=73.D0
 
!  n_multi_cut=5
 n_multi_cut=n_multi_cut-1
 
 
 
 OPEN(50,file='fit_lumass_2mass.txt')
READ(50,*) mfit_a,mfit_b,mfit_c
READ(50,*)
READ(50,*)
 CLOSE(50)
 
  OPEN(50,file='globalparameters_2MRS.txt')
READ(50,*) basic_link
READ(50,*)
READ(50,*) av_rad_pec_vel
READ(50,*) avlogmass
 CLOSE(50)

 
 
q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)
 cV_mill=(380.D6/h)**3
m_expected=cV_mill*rho_crit*Omega_m
part_exp=(2160.D0**3)*((330.D0/500.D0)**3)
mag_sol_r=4.76D0
sdss_cover=9274.D0/41253.D0
twomass_cover=0.91D0
mock_cover=1.D0/8.D0
area_cover=mock_cover


apx=TRIM(adjustl(apx))
deltamag=(-magmax+magmin)/1500.D0

OPEN(50,file='exec_counter_grid.txt')
READ(50,*) exec_counter
 CLOSE(50)


WRITE(ecs,*) exec_counter
ecs=adjustl(ecs)


WRITE(nm,*) uuu
nm=adjustl(nm)
 
 

 
!   WRITE(*,*) 'now in CPU',uuu
 
 
 

filename=TRIM(apx)//'_mock'//TRIM(nm)//'/'//TRIM(apx)//'_galaxies_mock'//TRIM(nm)//'_final.txt'
! get length of file
OPEN(50,file=TRIM(filename))
io_err=0
n=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n=n+1
END DO
 CLOSE(50)
n=n-1

! WRITE(*,*) n,'galaxies will be used in',uuu


listname=TRIM(apx)//'_mock'//TRIM(nm)//'/'//TRIM(apx)//'_galaxies_mock'//TRIM(nm)//'_final_simid.txt'
OPEN(50,file=TRIM(listname))
io_err=0
n_l=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_l=n_l+1
END DO
 CLOSE(50)
n_l=n_l-1




allocate(mmockid(1:n_l))
allocate(mockid(1:n_l))
allocate(mockfofid(1:n_l))
allocate(mockfofid_act(1:n_l))
allocate(trueclusterlist(1:n_l))
allocate(trueclustermembers(1:n_l))
allocate(trueclusterfofid(1:n_l))
allocate(group_bij(1:n_l))
allocate(samecluster(1:n_l))
allocate(P_mock(1:n_l))
allocate(P_fof(1:n_l))
allocate(clustergalnmember(1:n))
allocate(truegalnmember(1:n_l))
allocate(clusterlist_ntrue(1:n_l))
allocate(truelist_ncluster(1:n_l))

allocate(ra(1:n))
allocate(dec(1:n))
allocate(z(1:n))
allocate(mag_g(1:n))
allocate(mag_r(1:n))

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
allocate(vred(1:n))
allocate(members_clustercore(1:n))
allocate(fof(1:n))
allocate(fof_recursive(1:n))
allocate(iter_help(1:n))
allocate(dummy_core(1:n))
allocate(clusterindex(1:n))
allocate(lum_gal(1:n))



! read file
OPEN(50,file=TRIM(filename))
DO i=1,n
READ(50,*) ra(i),dec(i),z(i),mag_g(i),mag_r(i)
END DO
CLOSE(50)




OPEN(50,file=TRIM(listname))
DO i=1,n_l
READ(50,*) mmockid(i),mockid(i),mockfofid(i)
END DO
CLOSE(50)



DO i=1,n_l
mockfofid_act(i)=.TRUE.
END DO



n_truegroups_mock=0
n_truegroups_all=0
n2plus_truemock=0
n2plus_trueall=0
DO i=1,n_l
trueclusterlist(i)=0
trueclustermembers(i)=0
trueclusterfofid(i)=0
END DO



DO i=1,n_l
IF (mockfofid_act(i)) THEN
n_truegroups_mock=n_truegroups_mock+1
trueclusterlist(i)=n_truegroups_mock
trueclustermembers(n_truegroups_mock)=trueclustermembers(n_truegroups_mock)+1
trueclusterfofid(n_truegroups_mock)=mockfofid(n_truegroups_mock)
mockfofid_act(i)=.FALSE.

flipvar=.FALSE.

IF (i<n_l) THEN
DO ii=i+1,n_l
IF (mockfofid(i)==mockfofid(ii)) THEN
trueclusterlist(ii)=n_truegroups_mock
trueclustermembers(n_truegroups_mock)=trueclustermembers(n_truegroups_mock)+1
mockfofid_act(ii)=.FALSE.
flipvar=.TRUE.
END IF 
END DO
END IF

IF (flipvar) THEN
n2plus_truemock=n2plus_truemock+1
END IF

END IF
END DO



! WRITE(*,*) n2plus_truemock,'of',n_truegroups_mock,'mock groups have more than 2 members'





DO i=1,n
 clusterindex(i)=0
ra_rad(i)=ra(i)*PI/180.D0
dec_rad(i)=dec(i)*PI/180.D0
END DO

! WRITE(*,*) 'initialisation done'

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
! lum_g(i)=10.D0**(-0.4D0*(mag_g(i)-mag_sol_g)) 
mass_halo_fit=mfit_a*(LOG10(lum(i)))**2+mfit_b*(LOG10(lum(i)))+mfit_c
vred(i)=((((1.D0+z(i))**2)-1.D0)/(((1.D0+z(i))**2)+1.D0))*light
avmass=10.D0**avlogmass
mass_halo_fit=10.D0**mass_halo_fit
lum_gal(i)=lum(i)
! mass_stretch(i)=mass_halo_fit/avmass
END IF
END DO
WRITE(*,*) 'rescaling linking length according to the masses'






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


! (mass_stretch(i)**s_power), (mass_stretch(i)**s_power)

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
 
 IF (uuu==1) THEN
 OPEN(50,file='exec_counter_grid.txt')
WRITE(50,*) (exec_counter+1)
 CLOSE(50)
 END IF
 
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
! lumotherband(counter)=0.D0
! sigma_clustercore(counter)=0.D0
! angradius_clustercore(counter)=0.D0
! lumdistance_clustercore(counter)=0.D0
! mass_clustercore(counter)=0.D0
members_clustercore(counter)=0
! mdyn_clustercore(counter)=0.D0
!  colour_clustercore(i)=0.D0

! using iter method from Robotham to find cluster centre


help_index=0

DO ii=1,n
IF (fof(ii)) THEN
members_clustercore(counter)=members_clustercore(counter)+1
 clusterindex(ii)=counter
lumtot_clustercore(counter)=lumtot_clustercore(counter)+(lum(ii))
! lumotherband(counter)=lumotherband(counter)+(lum_g(ii))
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
! angradius_clustercore(counter)=(dummy_core(help_index)+dummy_core(help_index+1))/2.D0
ELSE
help_index=help_count/2
! angradius_clustercore(counter)=dummy_core(help_index+1)
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

! sigma_gap(counter)=0.D0
! 
! DO ii=1,(n_help-1)
! weight=DBLE(ii*(n_help-ii))
! sigma_gap(counter)=sigma_gap(counter)+weight*ABS(dummy_core(ii+1)-dummy_core(ii))
! END DO
! 
! sigma_gap(counter)=sigma_gap(counter)*SQRT(PI)/DBLE(n_help*(n_help-1))
! 
! IF (sigma_gap(counter)<red_measure) THEN
! sigma_gap(counter)=red_measure
! END IF
! 
! sigma_clustercore(counter)=SQRT(DBLE(n_help)/DBLE(n_help-1)*(sigma_gap(counter)**2)-(red_measure**2))
! 
! 
! IF (sigma_clustercore(counter)<red_measure) THEN
! sigma_clustercore(counter)=red_measure
! END IF
! 
! 
! 


ra_clustercore(counter)=ra(central_index)
dec_clustercore(counter)=dec(central_index)



ELSE 


ra_clustercore(counter)=ra(help_index)
dec_clustercore(counter)=dec(help_index)
z_clustercore(counter)=z(help_index)
END IF


END IF
END DO




! 
! 
! OPEN(50,file='results_'//TRIM(apx)//'/cluster_list_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
! DO i=1,counter
! WRITE(50,*) ra_clustercore(i),dec_clustercore(i),z_clustercore(i),&
! lumtot_clustercore(i),sigma_clustercore(i),members_clustercore(i)
! END DO
!  CLOSE(50)




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


OPEN(50,file='gridquality_'//TRIM(apx)//'/fine_binclusterfof_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
DO i=1,102
WRITE(50,*) i,binclusterfof(i)
END DO
 CLOSE(50)
 

 



DO i=1,102
binclustermock(i)=0
END DO

DO i=1,100
DO ii=1,n_truegroups_mock
IF (trueclustermembers(ii)==i) THEN
binclustermock(i)=binclustermock(i)+1
END IF
END DO
END DO


DO ii=1,n_truegroups_mock
IF ((trueclustermembers(ii)>100).AND.(trueclustermembers(ii)<1000)) THEN
binclustermock(101)=binclustermock(101)+1
END IF
END DO

DO ii=1,n_truegroups_mock
IF (trueclustermembers(ii)>1000) THEN
binclustermock(102)=binclustermock(102)+1
END IF
END DO

OPEN(50,file='gridquality_'//TRIM(apx)//'/fine_binclustermock_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
DO i=1,102
WRITE(50,*) i,binclustermock(i)
END DO
 CLOSE(50)


 n_g2_fof=0
 n_g2_mock=0
 n_g2_bij=0
 
DO i=(n_multi_cut+1),102
 n_g2_fof=n_g2_fof+binclusterfof(i)
 n_g2_mock=n_g2_mock+binclustermock(i)
END DO
 
 
 
 DO i=1,n_l
 group_bij(i)=.FALSE.
 END DO 

 
!   OPEN(50,file='helpcount.txt')
  
  
 DO i=1,n 
 IF (oldact(i)) THEN
 IF (members_clustercore(clusterindex(i))>n_multi_cut) THEN
 
 help_count=0
 
 
 
 
 DO ii=1,n
 IF (oldact(ii)) THEN
 IF (trueclustermembers(trueclusterlist(ii))>n_multi_cut) THEN 
 
 
 IF (clusterindex(i)==clusterindex(ii)) THEN
 IF (trueclusterlist(i)==trueclusterlist(ii)) THEN 

 help_count=help_count+1
 
 END IF
 END IF
 
 

 END IF
 END IF
 END DO
 
!  WRITE(50,*) help_count
 IF ((2*help_count)>trueclustermembers(trueclusterlist(i))) THEN
 IF ((2*help_count)>members_clustercore(clusterindex(i))) THEN
 group_bij(clusterindex(i))=.TRUE. 
 END IF
 END IF
 

 END IF
 END IF
 END DO
 
!   CLOSE(50)
 
 DO i=1,n_l
 IF (group_bij(i)) THEN
 n_g2_bij=n_g2_bij+1
 END IF
 END DO
 
 OPEN(50,file='gridquality_'//TRIM(apx)//'/fine_n_g_multi_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
WRITE(50,*) n_g2_bij,n_g2_fof,n_g2_mock
WRITE(50,*) counter,n_truegroups_mock
 CLOSE(50)
 
 E_fof=DBLE(n_g2_bij)/DBLE(n_g2_fof)
 E_mock=DBLE(n_g2_bij)/DBLE(n_g2_mock)
 E_tot=E_fof*E_mock
 
 
 
   DO ii=1,n
 IF (oldact(ii)) THEN
 clustergalnmember(ii)=members_clustercore(clusterindex(ii))
 truegalnmember(ii)=trueclustermembers(trueclusterlist(ii))
 END IF
END DO
 
 
  DO i=1,n 
 clusterlist_ntrue(i)=-1
 truelist_ncluster(i)=-1
 
  END DO
  
 DO i=1,n 
 IF (oldact(i)) THEN
 IF (clusterlist_ntrue(i)==-1) THEN
 
 n_help=0
 DO ii=i,n
 IF (oldact(ii)) THEN
 IF (clusterlist_ntrue(ii)==-1) THEN
 IF (clusterindex(i)==clusterindex(ii)) THEN
 IF (trueclusterlist(i)==trueclusterlist(ii)) THEN
 n_help=n_help+1
 clusterlist_ntrue(ii)=-2
 END IF
 END IF
 END IF
 END IF
 END DO
 
 DO ii=i,n
 IF (clusterlist_ntrue(ii)==-2) THEN
 clusterlist_ntrue(ii)=n_help
 END IF
 END DO
 
 END IF
 END IF
 END DO

 

 DO i=1,counter
  P_fof(i)=0.D0
  IF (members_clustercore(i)>n_multi_cut) THEN
   
  DO ii=1,n
  IF (oldact(ii)) THEN
  IF (clusterindex(ii)==i) THEN
  
  P_help=(DBLE(clusterlist_ntrue(ii))**2)/&
  (DBLE(members_clustercore(clusterindex(ii)))*DBLE(trueclustermembers(trueclusterlist(ii))))
   IF (P_help>P_fof(i)) THEN
   P_fof(i)=P_help
   END IF   
   END IF
  END IF
  END DO 

  END IF
  END DO
 
  DO i=1,n_truegroups_mock
  P_mock(i)=0.D0
  IF (trueclustermembers(i)>n_multi_cut) THEN
   
  DO ii=1,n
  IF (oldact(ii)) THEN
  IF (trueclusterlist(ii)==i) THEN
  
  P_help=(DBLE(clusterlist_ntrue(ii))**2)/&
  (DBLE(members_clustercore(clusterindex(ii)))*DBLE(trueclustermembers(trueclusterlist(ii))))
   IF (P_help>P_mock(i)) THEN
   P_mock(i)=P_help
   END IF   
   END IF
  END IF
  END DO 

  END IF
  END DO
 
 

    Q_fof=0.D0
    divisor=0.D0 
    DO i=1,counter
    IF (members_clustercore(i)>n_multi_cut) THEN    
    Q_fof=Q_fof+P_fof(i)*DBLE(members_clustercore(i))
    divisor=divisor+DBLE(members_clustercore(i))
    END IF
    END DO
    Q_fof=Q_fof/divisor

   
   
   
    Q_mock=0.D0
    divisor=0.D0
    DO i=1,n_truegroups_mock
    IF (trueclustermembers(i)>n_multi_cut) THEN
    Q_mock=Q_mock+P_mock(i)*DBLE(trueclustermembers(i))
    divisor=divisor+DBLE(trueclustermembers(i))
    END IF
    END DO  
    Q_mock=Q_mock/divisor
   
   
!    
! OPEN(50,file='results_'//TRIM(apx)//'/P_fof_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
! DO i=1,counter
! WRITE(50,*) P_fof(i)
! END DO
! CLOSE(50)
! 
! OPEN(50,file='results_'//TRIM(apx)//'/P_mock_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
! DO i=1,n_truegroups_mock
! WRITE(50,*) P_mock(i)
! END DO
! CLOSE(50)
! 

   Q_tot=Q_fof*Q_mock
   
   
   S_tot=E_tot*Q_tot
 
 

OPEN(50,file='gridquality_'//TRIM(apx)//'/fine_group_cost_function_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
WRITE(50,*) b_factor,R_factor,lum_power
WRITE(50,*) E_fof,E_mock,E_tot
WRITE(50,*) Q_fof,Q_mock,Q_tot
WRITE(50,*) S_tot
 CLOSE(50)
 
 
 
!  WRITE(*,*) 'quality estimated'
 

 
deallocate(mmockid)
deallocate(mockid)
deallocate(mockfofid)
deallocate(mockfofid_act)
deallocate(trueclusterlist)
deallocate(trueclustermembers)
deallocate(trueclusterfofid)
deallocate(group_bij)
deallocate(samecluster)
deallocate(P_mock)
deallocate(P_fof)
deallocate(clustergalnmember)
deallocate(truegalnmember)
deallocate(clusterlist_ntrue)
deallocate(truelist_ncluster)

deallocate(ra)
deallocate(dec)
deallocate(z)
deallocate(mag_g)
deallocate(mag_r)

deallocate(active)
deallocate(oldact)
deallocate(volume_weight)
deallocate(vol_int)
deallocate(lum_dist)
deallocate(b_eff)
deallocate(R_eff)
! deallocate(mass_stretch)
deallocate(angular_dist)
deallocate(basic_link_a)
deallocate(basic_link_R)
deallocate(ra_rad)
deallocate(dec_rad)
deallocate(lum)
deallocate(ra_clustercore)
deallocate(dec_clustercore)
deallocate(z_clustercore)
deallocate(lumtot_clustercore)
deallocate(vred)
deallocate(members_clustercore)
deallocate(fof)
deallocate(fof_recursive)
deallocate(iter_help)
deallocate(dummy_core)
deallocate(clusterindex)
deallocate(lum_gal)
 
 
 calculategroupcostfunction=S_tot






RETURN
END FUNCTION calculategroupcostfunction

 
 
 