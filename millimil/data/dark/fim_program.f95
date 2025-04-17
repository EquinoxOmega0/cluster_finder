PROGRAM millimilhalo
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,n_h,i,ii,n_p,uuu,anpfi,npartalt,n_f,n_g,n_sub,idnotfound,maxloop,n_rand,n_2mrs,n_snap
integer :: n_inside,n_fi,n_critdens,n_200,n_100,n_rhm10,doshells,help_count,galactive,ee,n_extra
character(200) :: filename,intial_values
integer(kind=8) :: dummy_int

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor,mparticles,part_tot
double precision :: dist,ts_modificator,d,mnew,part_in_empty,mag_sol_r,darkmass,avlogmass
double precision :: nextdist,avnextdist,helpdist,av_velrad,avdist,g_part_tot,cg_part_tot,mag_sol_Ks

double precision, allocatable ::  hx(:),hy(:),hz(:),h_np(:),h_r(:),px(:),py(:),pz(:),hm(:),hr2(:),h_mcrit(:)
double precision, allocatable ::r_fi(:),hmfi(:),hmfi_first(:),p_weights(:),hmass_sum(:)
logical, allocatable :: inside(:),in_fi(:),active(:)
integer(kind=8), allocatable :: h_fofid(:),g_id(:),g_fofid(:),sub_fofid(:),sub_hid(:)
double precision, allocatable :: gx(:),gy(:),gz(:),gvelX(:),gvelY(:),gvelZ(:),gnp(:)
double precision, allocatable :: ggDust(:),grDust(:),giDust(:),gKsmag(:),gJmag(:)

double precision, allocatable :: fof_group_totlum(:),g_lum(:),cg_lum(:),cg_fof_mass(:),cg_fof_rad(:)
double precision, allocatable :: gal_rad_fi(:),magapp(:)
double precision, allocatable :: gal_mass(:),gmfi(:)
logical, allocatable :: gal_in_fi(:)
integer, allocatable :: fof_group_members(:)
logical, allocatable :: gal_fof_sel(:),fof_unused(:),near_central_10(:),activeg(:)
double precision, allocatable ::  extra_x(:),extra_y(:),extra_z(:),extra_r(:)
double precision, allocatable :: extra_r_fi(:)
double precision, allocatable :: M_fi_rad(:)

integer, dimension(0:500,0:500) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y,hcount,no_iter,i1,i2,i3

double precision, dimension(0:110) :: bin
double precision :: d_J,dummy_coeffient,d_Ks,e_J,e_Ks,f_J,f_Ks,sidelength

double precision, allocatable :: rand_x(:),rand_y(:),rand_z(:)
logical, allocatable :: rand_in_fi(:)
logical, allocatable :: rand_in_fi_iter(:)
integer :: nrand_fi,nrand_fi_iter

integer, dimension(0:500) :: bin_numbers
double precision :: binhelp,redlimit,D_lum_limit
double precision, dimension(1:3,1:3,1:3) :: cubeshift_x,cubeshift_y,cubeshift_z

double precision, allocatable ::  deviation(:),otherside(:)
logical, allocatable :: activef(:)
double precision :: detA,detB1,detB2,detB3
double precision :: a_fit,b_fit,c_fit,rms
double precision :: delta_a_fit,delta_b_fit,delta_c_fit
double precision :: a_err,b_err,c_err
logical :: continueloop
double precision, dimension(1:3,1:3) :: Amatrix,Bmatrix1,Bmatrix2,Bmatrix3
double precision, dimension(1:3) :: Yvector
integer, dimension (1:6) :: snapnum
double precision, dimension (1:6) :: redatsnap


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme FINITE INFINITY CORRECTIONS started'
WRITE(*,*) '============================================================'

 CALL SYSTEM_CLOCK(hcount)
 hcount=hcount-INT(DBLE(hcount)/100000.D0)*100000
 CALL srand(hcount) 

!  OPEN(52,file='2mass_sdss_transformation3.txt')
! READ(52,*) 
! READ(52,*) d_J,dummy_coeffient,d_Ks
! READ(52,*) e_J,dummy_coeffient,e_Ks
! READ(52,*) f_J,dummy_coeffient,f_Ks
!  CLOSE(52)
 
   OPEN(52,file='n.txt')
READ(52,*) n_snap
 CLOSE(52)
 
  OPEN(52,file='snapshot_redshift_list.txt')
DO i=1,6
READ(52,*) snapnum(i),redatsnap(i)
END DO
 CLOSE(52)

 redlimit=0.D0
 DO i=1,6
 IF (n_snap==snapnum(i)) THEN
 redlimit=redatsnap(i)
 END IF
 END DO
 

! define constants
PI=ACOS(-1.D0)

maxloop=10
doshells=1
no_iter=0
Omega_m=0.25D0
Omega_l=0.75D0

sidelength=62.5D0
q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
H0=73.D0
h=H0/100.D0
G=4.302D-9 !in Mpc/Msol * (km/s)**2
rho_crit=3.D0*((H0)**2)/(8.D0*PI*G)
! WRITE(*,*) rho_crit
ts_modificator=(50.1D0/61.7D0)**2

 cV_mill=(sidelength/h)**3
m_expected=cV_mill*rho_crit*Omega_m

mag_sol_r=4.76D0
mag_sol_Ks=3.27D0
! WRITE(*,*) m_expected


part_exp=(270.D0**3)
mparticles=8.61D8/h

m_expected=mparticles*part_exp



! WRITE(*,*) m_expected

! get length of file
OPEN(50,file='fof.txt')
io_err=0
n_h=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_h=n_h+1
END DO
 CLOSE(50)
n_h=n_h-2

WRITE(*,*) n_h,'fof groups will be used'



! get length of file
OPEN(50,file='p.txt')
io_err=0
n_p=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_p=n_p+1
END DO
 CLOSE(50)
n_p=n_p-2

WRITE(*,*) n_p,'particles will be used'


OPEN(50,file='gal.txt')
io_err=0
n_g=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_g=n_g+1
END DO
 CLOSE(50)
n_g=n_g-2

WRITE(*,*) n_g,'galaxies will be used'



n_rand=1000000
allocate(rand_x(1:n_rand))
allocate(rand_y(1:n_rand))
allocate(rand_z(1:n_rand))
allocate(rand_in_fi(1:n_rand))
allocate(rand_in_fi_iter(1:n_rand))


allocate(h_fofid(1:n_h))
allocate(hx(1:n_h))
allocate(hy(1:n_h))
allocate(hz(1:n_h))
allocate(h_np(1:n_h))
allocate(h_r(1:n_h))
allocate(hm(1:n_h))
allocate(hr2(1:n_h))
allocate(h_mcrit(1:n_h))

allocate(r_fi(1:n_h))
allocate(hmass_sum(1:n_h))

allocate(fof_group_totlum(1:n_h))
allocate(fof_group_members(1:n_h))

allocate(gal_fof_sel(1:n_g))
allocate(g_lum(1:n_g))
allocate(magapp(1:n_g))


allocate(in_fi(1:n_p))


allocate(px(1:n_p))
allocate(py(1:n_p))
allocate(pz(1:n_p))
allocate(inside(1:n_p))
allocate(p_weights(1:n_p))

allocate(active(1:n_h))
allocate(hmfi(1:n_h))
allocate(hmfi_first(1:n_h))


allocate(g_id(1:n_g))
allocate(g_fofid(1:n_g))
allocate(gx(1:n_g))
allocate(gy(1:n_g))
allocate(gz(1:n_g))
allocate(gvelX(1:n_g))
allocate(gvelY(1:n_g))
allocate(gvelZ(1:n_g))
allocate(gnp(1:n_g))

allocate(ggDust(1:n_g))
allocate(grDust(1:n_g))
allocate(giDust(1:n_g))
allocate(gKsmag(1:n_g))
allocate(gJmag(1:n_g))

allocate(M_fi_rad(1:n_h))



allocate(fof_unused(1:n_h))

allocate(gal_rad_fi(1:n_g))

allocate(gal_mass(1:n_g))

allocate(gal_in_fi(1:n_p))

allocate(activeg(1:n_g))
allocate(gmfi(1:n_g))

WRITE(*,*) 'arrays allocated'

DO i=1,n_rand
rand_x(i)=RAND()*sidelength/h
rand_y(i)=RAND()*sidelength/h
rand_z(i)=RAND()*sidelength/h

rand_in_fi(i)=.FALSE.
rand_in_fi_iter(i)=.FALSE.

END DO






! read file
OPEN(50,file='fof.txt')
READ(50,*)
DO i=1,n_h
READ(50,*) h_fofid(i),hx(i),hy(i),hz(i),h_np(i),h_mcrit(i),h_r(i)
END DO
CLOSE(50)


OPEN(50,file='gal.txt')
READ(50,*)
DO i=1,n_g
READ(50,*) g_id(i),g_fofid(i),gx(i),gy(i),gz(i),gvelX(i),gvelY(i),gvelZ(i),gnp(i),&
ggDust(i),grDust(i),giDust(i)
END DO
CLOSE(50)






! read file
OPEN(50,file='p.txt')
READ(50,*)
DO i=1,n_p
READ(50,*) px(i),py(i),pz(i)
p_weights(i)=1.D0
END DO
CLOSE(50)

WRITE(*,*) '------------------------------------------------------------'
WRITE(*,*) 'all data read in'
g_part_tot=0.D0
DO i=1,n_g

gx(i)=gx(i)/h
gy(i)=gy(i)/h
gz(i)=gz(i)/h

gal_mass(i)=gnp(i)*mparticles
g_part_tot=g_part_tot+h_np(i)

gal_rad_fi(i)=(3.D0/(4.D0*PI)*gal_mass(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)
gal_rad_fi(i)=gal_rad_fi(i)**2


END DO




DO i=1,n_p
gal_in_fi(i)=.FALSE.

 END DO
 
 
 
 
part_tot=0.D0
DO i=1,n_h
hx(i)=hx(i)/h
hy(i)=hy(i)/h
hz(i)=hz(i)/h
part_tot=part_tot+h_np(i)
hm(i)=h_np(i)*mparticles
h_r(i)=h_r(i)/h




r_fi(i)=(3.D0/(4.D0*PI)*hm(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)
! r_fi(i)=r_fi(i)**2




END DO


DO i=1,n_p
px(i)=px(i)/h
py(i)=py(i)/h
pz(i)=pz(i)/h
inside(i)=.FALSE.
in_fi(i)=.FALSE.

END DO

WRITE(*,*) part_tot/part_exp*100.D0,'% of all particles are in halos >=20 particles'



WRITE(*,*) g_part_tot/part_exp*100.D0,'% of all particles are in galaxies'



divisor=(SQRT(1.D0+2.D0*q0*redlimit)+1.D0+q0*redlimit)
!luminosity distance
D_lum_limit=light/H0*redlimit*(1.D0+((redlimit*(1.D0-q0))/divisor))


DO i=1,n_g
gal_fof_sel(i)=.TRUE.
g_lum(i)=10.D0**(-0.4D0*(grDust(i)-mag_sol_r))
magapp(i)=grDust(i)+5.D0*LOG10(D_lum_limit*1.D6)-5.D0
END DO

DO i=1,n_g

IF (magapp(i)>17.77D0) THEN
gal_fof_sel(i)=.FALSE.
END IF

END DO
 
 
DO i=1,n_h
fof_group_totlum(i)=0.D0
fof_group_members(i)=0
DO ii=1,n_g
IF (gal_fof_sel(ii)) THEN

IF (h_fofid(i)==g_fofid(ii)) THEN

gal_fof_sel(ii)=.FALSE.
fof_group_members(i)=fof_group_members(i)+1
fof_group_totlum(i)=fof_group_totlum(i)+g_lum(ii)

END IF

END IF
END DO




END DO

help_count=n_h
darkmass=0.D0

DO i=1,n_h
IF (fof_group_members(i)==0) THEN
help_count=help_count-1
darkmass=darkmass+hm(i)
END IF
END DO


WRITE(*,*) 'mass in dark halos:', darkmass/m_expected,'%'
WRITE(*,*) help_count,'luminous halos'
DO i=1,n_h
fof_unused(i)=.TRUE.
END DO



 DO i1=1,3
 DO i2=1,3
 DO i3=1,3
 cubeshift_x(i1,i2,i3)=DBLE(i1-2)*sidelength/h
 cubeshift_y(i1,i2,i3)=DBLE(i2-2)*sidelength/h
 cubeshift_z(i1,i2,i3)=DBLE(i3-2)*sidelength/h
 END DO
 END DO
 END DO


DO i=1,n_h
DO ii=1,n_p

dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)


IF (in_fi(ii).EQV..FALSE.) THEN
IF (dist<(r_fi(i)**2)) THEN
in_fi(ii)=.TRUE.
END IF
END IF


END DO


DO ii=1,n_rand

 DO i1=1,3
 DO i2=1,3
 DO i3=1,3
  
dist=((hx(i)-rand_x(ii)-cubeshift_x(i1,i2,i3))**2)+&
((hy(i)-rand_y(ii)-cubeshift_y(i1,i2,i3))**2)+((hz(i)-rand_z(ii)-cubeshift_z(i1,i2,i3))**2)


IF (rand_in_fi(ii).EQV..FALSE.) THEN
IF (dist<(r_fi(i)**2)) THEN
rand_in_fi(ii)=.TRUE.
END IF
END IF

 END DO
 END DO
 END DO


 
 

END DO


END DO




n_fi=0

DO i=1,n_p

IF (in_fi(i)) THEN
n_fi=n_fi+1
END IF

END DO



nrand_fi=0


DO i=1,n_rand


IF (rand_in_fi(i)) THEN
nrand_fi=nrand_fi+1
END IF

END DO





OPEN(50,file='results_fi.txt')
WRITE(50,*) part_tot/part_exp*100.D0,'% of all particles are in halos >=20 particles'
WRITE(50,*) '-------------------------------------------'
WRITE(50,*) 'finite infinity: percentage',n_fi,(100.D0*DBLE(n_fi)/n_p)
WRITE(50,*) 'volume percentage:',nrand_fi,(100.D0*DBLE(nrand_fi)/n_rand)

WRITE(*,*) part_tot/part_exp*100.D0,'% of all particles are in halos >=20 particles'
WRITE(*,*) '-------------------------------------------'
WRITE(*,*) 'finite infinity: percentage',n_fi,(100.D0*DBLE(n_fi)/n_p)
WRITE(*,*) 'volume percentage:',nrand_fi,(100.D0*DBLE(nrand_fi)/n_rand)


DO i=1,n_h

active(i)=.TRUE. 

! r_fi(i)=SQRT(r_fi(i))
hmfi(i)=hm(i)
hmfi_first(i)=hm(i)
hmass_sum(i)=hm(i)
END DO


anpfi=0
npartalt=part_tot




uuu=0
DO WHILE (uuu<10)
uuu=uuu+1
WRITE(*,*) '-------------------------------------------'
WRITE(50,*) '-------------------------------------------'
WRITE(*,*) uuu,'iteration for fi'
WRITE(50,*) uuu,'iteration for fi'

IF (uuu>2) THEN
IF ((DBLE(anpfi)/DBLE(n_fi))<1.001D0) THEN
IF ((DBLE(anpfi)/DBLE(n_fi))>0.999D0) THEN
uuu=11
END IF
END IF
END IF


! filter fi regions which are inside of other fi regions
DO i=1,n_h
IF (active(i)) THEN
!WRITE(*,*) i

ii=0
DO WHILE (ii<n_h)
ii=ii+1
IF (active(ii)) THEN


IF (i/=ii) THEN
d=(hx(i)-hx(ii))**2+(hy(i)-hy(ii))**2+(hz(i)-hz(ii))**2
d=SQRT(d)
IF ((d+r_fi(ii))<r_fi(i)) THEN
active(ii)=.FALSE.

hmass_sum(i)=hmass_sum(i)+hmass_sum(ii)
mnew=hmfi(i)+hmfi(ii)
hmfi_first(i)=hmfi_first(i)+hmfi_first(ii)
hx(i)=(hx(i)*hmfi(i)+hx(ii)*hmfi(ii))/mnew
hy(i)=(hy(i)*hmfi(i)+hy(ii)*hmfi(ii))/mnew
hz(i)=(hz(i)*hmfi(i)+hz(ii)*hmfi(ii))/mnew
hmfi(i)=mnew



r_fi(i)=(3.D0/(4.D0*PI)*hmfi(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)

ii=0
END IF
END IF
END IF

END DO

END IF
END DO


ii=0
DO i=1,n_h
IF (active(i)) THEN
ii=ii+1
! r_fi(i)=r_fi(i)**2
END IF
END DO


WRITE(*,*) ii,'halos will be used after the filtering'
WRITE(50,*) ii,'halos will be used after the filtering'







DO i=1,n_p
in_fi(i)=.FALSE.
p_weights(i)=0.D0
END DO


DO i=1,n_h
IF (active(i)) THEN
!WRITE(*,*) i
DO ii=1,n_p
! IF (in_fi(ii).EQV..FALSE.) THEN
dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)

IF (dist<(r_fi(i)**2)) THEN
in_fi(ii)=.TRUE.
p_weights(ii)=p_weights(ii)+1.D0
END IF
! END IF

END DO
END IF
END DO


anpfi=n_fi
n_fi=0

DO i=1,n_p

IF (in_fi(i)) THEN
n_fi=n_fi+1
p_weights(i)=1.D0/p_weights(i)
END IF


END DO




WRITE(*,*) 'inside finite infinity:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
WRITE(50,*) 'inside finite infinity:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))


DO i=1,n_h
IF (active(i)) THEN

hmfi(i)=0.D0
DO ii=1,n_p
IF (in_fi(ii)) THEN
dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)
IF (dist<(r_fi(i)**2)) THEN
hmfi(i)=hmfi(i)+p_weights(ii)*mparticles
END IF
END IF
END DO
r_fi(i)=(3.D0/(4.D0*PI)*hmfi(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)
END IF
END DO

! npartalt=n_fi

IF (uuu==1) THEN
OPEN(51,file='first_masses_fi.txt')
DO i=1,n_h
IF (active(i)) THEN
WRITE(51,*) hmass_sum(i),hmfi(i)
hmfi_first(i)=hmfi(i)
END IF
END DO
CLOSE(51)
END IF





END DO



WRITE(*,*) '-------------------------------------------'
WRITE(50,*) '-------------------------------------------'




DO i=1,n_h
IF (active(i)) THEN

DO ii=1,n_rand

 DO i1=1,3
 DO i2=1,3
 DO i3=1,3
 
 


dist=((hx(i)-rand_x(ii)-cubeshift_x(i1,i2,i3))**2)+&
((hy(i)-rand_y(ii)-cubeshift_y(i1,i2,i3))**2)+((hz(i)-rand_z(ii)-cubeshift_z(i1,i2,i3))**2)



IF (rand_in_fi_iter(ii).EQV..FALSE.) THEN
IF (dist<(r_fi(i)**2)) THEN
rand_in_fi_iter(ii)=.TRUE.
END IF
END IF

END DO
END DO
END DO

END DO
END IF
END DO
nrand_fi_iter=0


DO i=1,n_rand

IF (rand_in_fi_iter(i)) THEN
nrand_fi_iter=nrand_fi_iter+1
END IF


END DO


OPEN(51,file='final_masses_fi.txt')
DO i=1,n_h
IF (active(i)) THEN
WRITE(51,*) hmass_sum(i),hmfi(i)
END IF
END DO
CLOSE(51)

OPEN(51,file='final_from_first_masses_fi.txt')
DO i=1,n_h
IF (active(i)) THEN
WRITE(51,*) hmfi_first(i),hmfi(i)
END IF
END DO
CLOSE(51)




CLOSE(50)





OPEN(51,file='volume_percentage_fi.txt')
WRITE(51,*) 'in fi:',(DBLE(nrand_fi)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in fi iter:',(DBLE(nrand_fi_iter)/DBLE(n_rand)*100.D0),'%'
CLOSE(51)







WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


