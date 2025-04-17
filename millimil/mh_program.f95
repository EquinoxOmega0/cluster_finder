PROGRAM millimilhalo
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,n_h,i,ii,n_p,uuu,anpfi,npartalt,n_f,n_g,n_cg,n_sub,idnotfound,maxloop,n_rand,n_2mrs
integer :: n_inside,n_fi,n_critdens,n_200,n_100,n_rhm10,doshells,help_count,galactive,ee,n_extra
character(200) :: filename,intial_values
integer(kind=8) :: dummy_int

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor,mparticles,part_tot
double precision :: dist,ts_modificator,d,mnew,part_in_empty,mag_sol_r,darkmass,avlogmass
double precision :: nextdist,avnextdist,helpdist,av_velrad,avdist,g_part_tot,cg_part_tot,mag_sol_Ks


double precision, allocatable ::  hx(:),hy(:),hz(:),h_np(:),h_r(:),px(:),py(:),pz(:),hm(:),hr2(:),h_mcrit(:)
double precision, allocatable :: r_crit(:),r_fi(:),r_200(:),r_100(:),rhm10(:),hmfi(:)
logical, allocatable :: inside(:),in_fi(:),in_critdensity(:),in_200(:),in_100(:),in_rhm10(:),active(:)
integer(kind=8), allocatable :: h_fofid(:),g_id(:),g_fofid(:),cg_id(:),cg_fofid(:),sub_fofid(:),sub_hid(:)
double precision, allocatable :: fofnp(:),gx(:),gy(:),gz(:),gvelX(:),gvelY(:),gvelZ(:),gnp(:)
double precision, allocatable :: ggDust(:),grDust(:),giDust(:),gKsmag(:),gJmag(:)
double precision, allocatable :: cgx(:),cgy(:),cgz(:),cgvelX(:),cgvelY(:),cgvelZ(:),cgnp(:)
double precision, allocatable :: cggDust(:),cgrDust(:),cgiDust(:),cgKsmag(:),cgJmag(:)
double precision, allocatable :: fof_group_totlum(:),g_lum(:),cg_lum(:),cg_fof_mass(:),cg_fof_rad(:)
double precision, allocatable :: gal_rad_crit(:),gal_rad_fi(:),gal_rad_200(:),gal_rad_100(:),galcentral_crit(:)
double precision, allocatable :: gal_mass(:),gmfi(:)
logical, allocatable :: gal_in_fi(:),gal_in_critdensity(:),gal_in_200(:),gal_in_100(:),galcentral_in_crit(:)
integer, allocatable :: fof_group_members(:)
logical, allocatable :: gal_fof_sel(:),fof_unused(:),near_central_10(:),activeg(:)
double precision, allocatable ::  extra_x(:),extra_y(:),extra_z(:),extra_r(:)
double precision, allocatable ::  extra_r_crit(:),extra_r_fi(:),extra_r_200(:),extra_r_100(:),extra_rhm10(:)
double precision, allocatable :: M_crit_rad(:),M_fi_rad(:),M_100_rad(:),M_200_rad(:),Mr_rad(:),Mr10_rad(:)

integer, dimension(0:500,0:500) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y,hcount,no_iter,i1,i2,i3

double precision, dimension(0:110) :: bin
double precision :: d_J,dummy_coeffient,d_Ks,e_J,e_Ks,f_J,f_Ks,sidelength

double precision, allocatable :: rand_x(:),rand_y(:),rand_z(:)
logical, allocatable :: rand_in_fi(:),rand_in_critdensity(:),rand_in_200(:),rand_in_100(:)
logical, allocatable :: rand_in_fi_iter(:),rand_in_critdensity_iter(:),rand_in_rhm10(:),rand_inside(:)
integer :: nrand_inside,nrand_fi,nrand_critdens,nrand_200,nrand_100,nrand_rhm10,nrand_critdens_iter,nrand_fi_iter

integer, dimension(0:500) :: bin_numbers
double precision :: binhelp
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

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme MILLIMIL HALOS started'
WRITE(*,*) '============================================================'

 CALL SYSTEM_CLOCK(hcount)
 hcount=hcount-INT(hcount/100000)*100000
 CALL srand(hcount) 

 OPEN(52,file='2mass_sdss_transformation3.txt')
READ(52,*) 
READ(52,*) d_J,dummy_coeffient,d_Ks
READ(52,*) e_J,dummy_coeffient,e_Ks
READ(52,*) f_J,dummy_coeffient,f_Ks
 CLOSE(52)
 

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
ts_modificator=(48.2D0/61.7D0)**2

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
OPEN(50,file='fof_millimil.txt')
io_err=0
n_h=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_h=n_h+1
END DO
 CLOSE(50)
n_h=n_h-2

WRITE(*,*) n_h,'fof groups will be used'

OPEN(50,file='fof_empty.txt')
io_err=0
n_f=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_f=n_f+1
END DO
 CLOSE(50)
n_f=n_f-2

WRITE(*,*) n_f,'fof groups are empty'

! get length of file
OPEN(50,file='particles_millimil.txt')
io_err=0
n_p=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_p=n_p+1
END DO
 CLOSE(50)
n_p=n_p-2

WRITE(*,*) n_p,'particles will be used'


OPEN(50,file='allgal_millimil.txt')
io_err=0
n_g=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_g=n_g+1
END DO
 CLOSE(50)
n_g=n_g-2

WRITE(*,*) n_g,'galaxies will be used'

OPEN(50,file='centralgal_millimil.txt')
io_err=0
n_cg=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_cg=n_cg+1
END DO
 CLOSE(50)
n_cg=n_cg-2

WRITE(*,*) n_cg,'central galaxies will be used'


n_rand=1000000
allocate(rand_x(1:n_rand))
allocate(rand_y(1:n_rand))
allocate(rand_z(1:n_rand))
allocate(rand_in_fi(1:n_rand))
allocate(rand_in_critdensity(1:n_rand))
allocate(rand_in_200(1:n_rand))
allocate(rand_in_100(1:n_rand))
allocate(rand_in_fi_iter(1:n_rand))
allocate(rand_in_critdensity_iter(1:n_rand))
allocate(rand_in_rhm10(1:n_rand))
allocate(rand_inside(1:n_rand))


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
allocate(r_crit(1:n_h))
allocate(r_100(1:n_h))
allocate(r_200(1:n_h))
allocate(rhm10(1:n_h))

allocate(fof_group_totlum(1:n_h))
allocate(fof_group_members(1:n_h))

allocate(gal_fof_sel(1:n_g))
allocate(g_lum(1:n_g))

allocate(in_fi(1:n_p))
allocate(in_critdensity(1:n_p))
allocate(in_100(1:n_p))
allocate(in_200(1:n_p))
allocate(in_rhm10(1:n_p))


allocate(px(1:n_p))
allocate(py(1:n_p))
allocate(pz(1:n_p))
allocate(inside(1:n_p))


allocate(active(1:n_h))
allocate(hmfi(1:n_h))

allocate(fofnp(1:n_f))

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
allocate(M_crit_rad(1:n_h))
allocate(M_fi_rad(1:n_h))
allocate(M_100_rad(1:n_h))
allocate(M_200_rad(1:n_h))
allocate(Mr_rad(1:n_h))
allocate(Mr10_rad(1:n_h))

allocate(cg_id(1:n_cg))
allocate(cg_fofid(1:n_cg))
allocate(cgx(1:n_cg))
allocate(cgy(1:n_cg))
allocate(cgz(1:n_cg))
allocate(cgvelX(1:n_cg))
allocate(cgvelY(1:n_cg))
allocate(cgvelZ(1:n_cg))
allocate(cgnp(1:n_cg))

allocate(cggDust(1:n_cg))
allocate(cgrDust(1:n_cg))
allocate(cgiDust(1:n_cg))
allocate(cgKsmag(1:n_cg))
allocate(cgJmag(1:n_cg))

allocate(cg_lum(1:n_cg))
allocate(cg_fof_mass(1:n_cg))
allocate(cg_fof_rad(1:n_cg))
allocate(fof_unused(1:n_h))
allocate(near_central_10(1:n_p))

allocate(gal_rad_crit(1:n_g))
allocate(gal_rad_fi(1:n_g))
allocate(gal_rad_200(1:n_g))
allocate(gal_rad_100(1:n_g))
allocate(galcentral_crit(1:n_cg))
allocate(gal_mass(1:n_g))

allocate(gal_in_fi(1:n_p))
allocate(gal_in_critdensity(1:n_p))
allocate(gal_in_200(1:n_p))
allocate(gal_in_100(1:n_p))
allocate(galcentral_in_crit(1:n_p))

allocate(activeg(1:n_g))
allocate(gmfi(1:n_g))



allocate(deviation(1:n_cg))
allocate(otherside(1:n_cg))
allocate(activef(1:n_cg))

DO i=1,n_rand
rand_x(i)=RAND()*sidelength/h
rand_y(i)=RAND()*sidelength/h
rand_z(i)=RAND()*sidelength/h

rand_in_fi(i)=.FALSE.
rand_in_critdensity(i)=.FALSE.
rand_in_200(i)=.FALSE.
rand_in_100(i)=.FALSE.
rand_in_fi_iter(i)=.FALSE.
rand_in_critdensity_iter(i)=.FALSE.
rand_in_rhm10(i)=.FALSE.
rand_inside(i)=.FALSE.
END DO






! read file
OPEN(50,file='fof_millimil.txt')
READ(50,*)
DO i=1,n_h
READ(50,*) h_fofid(i),hx(i),hy(i),hz(i),h_np(i),h_mcrit(i),h_r(i)
END DO
CLOSE(50)


OPEN(50,file='allgal_millimil.txt')
READ(50,*)
DO i=1,n_g
READ(50,*) g_id(i),g_fofid(i),gx(i),gy(i),gz(i),gvelX(i),gvelY(i),gvelZ(i),gnp(i),&
ggDust(i),grDust(i),giDust(i)
END DO
CLOSE(50)

OPEN(50,file='centralgal_millimil.txt')
READ(50,*)
DO i=1,n_cg
READ(50,*) cg_id(i),cg_fofid(i),cgx(i),cgy(i),cgz(i),cgvelX(i),cgvelY(i),cgvelZ(i),cgnp(i),&
 cggDust(i),cgrDust(i),cgiDust(i)
END DO
CLOSE(50)



OPEN(50,file='fof_empty.txt')
READ(50,*)
DO i=1,n_f
READ(50,*) dummy_int,fofnp(i)
END DO
CLOSE(50)

part_in_empty=0.D0
DO i=1,n_f
part_in_empty=part_in_empty+fofnp(i)
END DO

! read file
OPEN(50,file='particles_millimil.txt')
READ(50,*)
DO i=1,n_p
READ(50,*) px(i),py(i),pz(i)
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

gal_rad_crit(i)=(3.D0/(4.D0*PI)*gal_mass(i)/(rho_crit))**(1.D0/3.D0)
gal_rad_crit(i)=gal_rad_crit(i)**2

gal_rad_100(i)=(3.D0/(4.D0*PI)*gal_mass(i)/(rho_crit*100.D0))**(1.D0/3.D0)
gal_rad_100(i)=gal_rad_100(i)**2

gal_rad_200(i)=(3.D0/(4.D0*PI)*gal_mass(i)/(rho_crit*200.D0))**(1.D0/3.D0)
gal_rad_200(i)=gal_rad_200(i)**2

END DO


 cg_part_tot=0.D0
DO i=1,n_cg
 cgx(i)=cgx(i)/h
 cgy(i)=cgy(i)/h
 cgz(i)=cgz(i)/h
 cg_part_tot=cg_part_tot+cgnp(i)
 
 galcentral_crit(i)=(3.D0/(4.D0*PI)*cgnp(i)*mparticles/(rho_crit))**(1.D0/3.D0)
 galcentral_crit(i)=galcentral_crit(i)**2
 
END DO





DO i=1,n_p
gal_in_fi(i)=.FALSE.
gal_in_critdensity(i)=.FALSE.
gal_in_200(i)=.FALSE.
gal_in_100(i)=.FALSE.
 galcentral_in_crit(i)=.FALSE.
 END DO
 
 
 
 
part_tot=0.D0
DO i=1,n_h
hx(i)=hx(i)/h
hy(i)=hy(i)/h
hz(i)=hz(i)/h
part_tot=part_tot+h_np(i)
hm(i)=h_np(i)*mparticles
h_r(i)=h_r(i)/h
hr2(i)=h_r(i)**2

rhm10(i)=h_r(i)*10.D0
rhm10(i)=rhm10(i)**2

r_fi(i)=(3.D0/(4.D0*PI)*hm(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)
r_fi(i)=r_fi(i)**2

r_crit(i)=(3.D0/(4.D0*PI)*hm(i)/(rho_crit))**(1.D0/3.D0)
r_crit(i)=r_crit(i)**2

r_100(i)=(3.D0/(4.D0*PI)*hm(i)/(rho_crit*100.D0))**(1.D0/3.D0)
r_100(i)=r_100(i)**2

r_200(i)=(3.D0/(4.D0*PI)*hm(i)/(rho_crit*200.D0))**(1.D0/3.D0)
r_200(i)=r_200(i)**2


END DO


DO i=1,n_p
px(i)=px(i)/h
py(i)=py(i)/h
pz(i)=pz(i)/h
inside(i)=.FALSE.
in_fi(i)=.FALSE.
in_critdensity(i)=.FALSE.
in_200(i)=.FALSE.
in_100(i)=.FALSE.
in_rhm10(i)=.FALSE.
END DO

WRITE(*,*) part_tot/part_exp*100.D0,'% of all particles are in halos >=20 particles'

WRITE(*,*) part_in_empty/part_exp*100.D0,'% of all particles are in <<empty>> halos'

WRITE(*,*) g_part_tot/part_exp*100.D0,'% of all particles are in galaxies'

WRITE(*,*) cg_part_tot/part_exp*100.D0,'% of all particles are in central galaxies'

avnextdist=0.D0
avdist=0.D0
DO i=1,n_g
nextdist=100.D0
DO ii=1,n_g
IF (ii.NE.i) THEN
helpdist=((gx(i)-gx(ii))**2)+((gy(i)-gy(ii))**2)+((gz(i)-gz(ii))**2)
avdist=avdist+SQRT(helpdist)
IF (helpdist<(nextdist**2)) THEN
nextdist=SQRT(helpdist)
END IF
END IF
END DO
avnextdist=avnextdist+nextdist
END DO
avnextdist=avnextdist/DBLE(n_g)
avdist=avdist/(DBLE(n_g)*DBLE(n_g-1))

av_velrad=0.D0
DO i=1,n_g
av_velrad=av_velrad+(gvelX(i)**2)+(gvelY(i)**2)+(gvelZ(i)**2)
END DO
av_velrad=SQRT(av_velrad/(3.D0*DBLE(n_g)))




DO i=0,400
bin_numbers(i)=0
END DO

DO i=1,n_g

DO ii=0,400
binhelp=DBLE(ii-200)*10.D0
IF ((gvelX(i)>(binhelp-5.D0)).AND.(gvelX(i)<(binhelp+5.D0))) THEN
bin_numbers(ii)=bin_numbers(ii)+1
END IF
IF ((gvelY(i)>(binhelp-5.D0)).AND.(gvelY(i)<(binhelp+5.D0))) THEN
bin_numbers(ii)=bin_numbers(ii)+1
END IF
IF ((gvelZ(i)>(binhelp-5.D0)).AND.(gvelZ(i)<(binhelp+5.D0))) THEN
bin_numbers(ii)=bin_numbers(ii)+1
END IF
END DO

END DO

OPEN(70,file='vel_bin_numbers.txt')
DO i=0,400
binhelp=DBLE(i-200)*10.D0
WRITE(70,*) binhelp,bin_numbers(i)
END DO
 CLOSE (70)





DO i=1,n_g
gal_fof_sel(i)=.TRUE.
g_lum(i)=10.D0**(-0.4D0*(grDust(i)-mag_sol_r))
END DO

OPEN(70,file='log_mass_log_lum.txt')
DO i=1,n_g
WRITE(70,*) LOG10(g_lum(i)),LOG10(gal_mass(i)),(ggDust(i)-grDust(i))
END DO
 CLOSE (70)

 
 
 
 
 
 
 

 
 
 
 
 
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

help_count=0
darkmass=0.D0

OPEN(51,file='fof_gal_relations.txt')
DO i=1,n_h
IF (fof_group_members(i)>0) THEN
WRITE(51,*) LOG10(hm(i)),LOG10(fof_group_totlum(i)),fof_group_members(i)
ELSE
help_count=help_count+1
darkmass=darkmass+hm(i)
END IF
END DO
CLOSE(51)


WRITE(*,*) help_count,'of ',n_h,'groups are dark'
WRITE(*,*) 'totalling',(darkmass*100.D0/m_expected),'% of the entire mass'

DO i=1,n_h
fof_unused(i)=.TRUE.
END DO



DO i=1,n_cg
 cg_lum(i)=10.D0**(-0.4D0*(cgrDust(i)-mag_sol_r)) 
 cg_fof_rad(i)=0.D0
 cg_fof_mass(i)=0.D0

DO ii=1,n_h
! IF (fof_unused(ii)) THEN
IF (h_fofid(ii)==cg_fofid(i)) THEN

! fof_unused(ii)=.FALSE.
 cg_fof_mass(i)=hm(ii)
 cg_fof_rad(i)=h_r(ii)
END IF
! END IF
END DO

END DO


help_count=0
OPEN(51,file='central_size_lum.txt')
DO i=1,n_cg
IF (cg_fof_mass(i)>0.D0) THEN
WRITE(51,*)  LOG10(cg_fof_rad(i)),LOG10(cg_lum(i)),LOG10(cg_fof_mass(i)),LOG10(cgnp(i)*mparticles)
ELSE
help_count=help_count+1
END IF
END DO
CLOSE(51)


WRITE(*,*) 'for',help_count,'of ',n_cg,'central galaxies, we do not find FOF groups'

DO i=1,n_cg
IF (cg_fof_mass(i)>0.D0) THEN
activef(i)=.TRUE.
else
activef(i)=.false.
END IF
IF (LOG10(cgnp(i)*mparticles)>12.5D0) THEN
activef(i)=.false.
END IF
IF (LOG10(cgnp(i)*mparticles)<10.5D0) THEN
activef(i)=.false.
END IF
END DO


galactive=0
DO i=1,n_cg
IF (activef(i)) THEN
galactive=galactive+1
END IF
END DO
WRITE(*,*) galactive

 a_fit=0.25D0
 b_fit=-3.75D0
 c_fit=25.D0


 continueloop=.TRUE.
ee=0

delta_a_fit=0.D0
delta_b_fit=0.D0
delta_c_fit=0.D0


! do loop for opitmization
DO WHILE (continueloop)
ee=ee+1



!preparing matrizes
DO i=1,3 
Yvector(i)=0.D0
DO ii=1,3
Amatrix(i,ii)=0.D0
Bmatrix1(i,ii)=0.D0
Bmatrix2(i,ii)=0.D0
Bmatrix3(i,ii)=0.D0
END DO
END DO

!calculating basic matrix und vector
DO i=1,n_cg
IF (activef(i)) THEN

Amatrix(1,1)=Amatrix(1,1)+(LOG10(cg_lum(i))**4)
Amatrix(2,2)=Amatrix(2,2)+(LOG10(cg_lum(i))**2)
Amatrix(3,3)=Amatrix(3,3)+1.D0
Amatrix(1,2)=Amatrix(1,2)+(LOG10(cg_lum(i))**3)
Amatrix(2,1)=Amatrix(2,1)+(LOG10(cg_lum(i))**3)
Amatrix(1,3)=Amatrix(1,3)+(LOG10(cg_lum(i))**2)
Amatrix(3,1)=Amatrix(3,1)+(LOG10(cg_lum(i))**2)
Amatrix(2,3)=Amatrix(2,3)+LOG10(cg_lum(i))
Amatrix(3,2)=Amatrix(3,2)+LOG10(cg_lum(i))
Yvector(1)=Yvector(1)+LOG10(cgnp(i)*mparticles)*(LOG10(cg_lum(i))**2)
Yvector(2)=Yvector(2)+LOG10(cgnp(i)*mparticles)*LOG10(cg_lum(i))
Yvector(3)=Yvector(3)+LOG10(cgnp(i)*mparticles)

END IF
END DO



detA=Amatrix(1,1)*Amatrix(2,2)*Amatrix(3,3)
detA=detA+Amatrix(2,1)*Amatrix(3,2)*Amatrix(1,3)
detA=detA+Amatrix(1,2)*Amatrix(2,3)*Amatrix(3,1)
detA=detA-Amatrix(1,3)*Amatrix(2,2)*Amatrix(3,1)
detA=detA-Amatrix(2,1)*Amatrix(1,2)*Amatrix(3,3)
detA=detA-Amatrix(2,3)*Amatrix(3,2)*Amatrix(1,1)


!calculate B matrizes
DO i=1,3
DO ii=1,3
Bmatrix1(i,ii)=Amatrix(i,ii)
Bmatrix2(i,ii)=Amatrix(i,ii)
Bmatrix3(i,ii)=Amatrix(i,ii)
END DO
END DO

DO i=1,3
Bmatrix1(1,i)=Yvector(i)
Bmatrix2(2,i)=Yvector(i)
Bmatrix3(3,i)=Yvector(i)
END DO



! get determinat of Bmatrizes



detB1=Bmatrix1(1,1)*Bmatrix1(2,2)*Bmatrix1(3,3)
detB1=detB1+Bmatrix1(2,1)*Bmatrix1(3,2)*Bmatrix1(1,3)
detB1=detB1+Bmatrix1(1,2)*Bmatrix1(2,3)*Bmatrix1(3,1)
detB1=detB1-Bmatrix1(1,3)*Bmatrix1(2,2)*Bmatrix1(3,1)
detB1=detB1-Bmatrix1(2,1)*Bmatrix1(1,2)*Bmatrix1(3,3)
detB1=detB1-Bmatrix1(2,3)*Bmatrix1(3,2)*Bmatrix1(1,1)
detB2=Bmatrix2(1,1)*Bmatrix2(2,2)*Bmatrix2(3,3)
detB2=detB2+Bmatrix2(2,1)*Bmatrix2(3,2)*Bmatrix2(1,3)
detB2=detB2+Bmatrix2(1,2)*Bmatrix2(2,3)*Bmatrix2(3,1)
detB2=detB2-Bmatrix2(1,3)*Bmatrix2(2,2)*Bmatrix2(3,1)
detB2=detB2-Bmatrix2(2,1)*Bmatrix2(1,2)*Bmatrix2(3,3)
detB2=detB2-Bmatrix2(2,3)*Bmatrix2(3,2)*Bmatrix2(1,1)
detB3=Bmatrix3(1,1)*Bmatrix3(2,2)*Bmatrix3(3,3)
detB3=detB3+Bmatrix3(2,1)*Bmatrix3(3,2)*Bmatrix3(1,3)
detB3=detB3+Bmatrix3(1,2)*Bmatrix3(2,3)*Bmatrix3(3,1)
detB3=detB3-Bmatrix3(1,3)*Bmatrix3(2,2)*Bmatrix3(3,1)
detB3=detB3-Bmatrix3(2,1)*Bmatrix3(1,2)*Bmatrix3(3,3)
detB3=detB3-Bmatrix3(2,3)*Bmatrix3(3,2)*Bmatrix3(1,1)


 a_fit=detB1/detA
 b_fit=detB2/detA
 c_fit=detB3/detA


!calculate right side of fundamental plane
DO i=1,n_cg
IF (activef(i)) THEN


otherside(i)=a_fit*((LOG10(cg_lum(i)))**2)+b_fit*LOG10(cg_lum(i))+c_fit


END IF
END DO

! propagation of error
!root mean square
rms=0.D0

DO i=1,n_cg
IF (activef(i)) THEN



deviation(i)=otherside(i)-LOG10(cgnp(i)*mparticles)

rms=rms+((deviation(i))**2)

END IF
END DO



rms=SQRT(rms/DBLE(galactive))


! error of fitting parameters


 a_err=(Amatrix(2,2)*Amatrix(3,3)-Amatrix(2,3)*Amatrix(3,2))/detA
 a_err=SQRT(a_err)
 b_err=(Amatrix(1,1)*Amatrix(3,3)-Amatrix(1,3)*Amatrix(3,1))/detA
 b_err=SQRT(b_err)
 c_err=(Amatrix(1,1)*Amatrix(2,2)-Amatrix(1,2)*Amatrix(2,1))/detA
 c_err=SQRT(c_err)


!see if another loop is still nessecary


delta_a_fit=delta_a_fit-a_fit
delta_b_fit=delta_b_fit-b_fit
delta_c_fit=delta_c_fit-c_fit
delta_a_fit=delta_a_fit*100.D0
delta_b_fit=delta_b_fit*100.D0
delta_c_fit=delta_c_fit*100.D0
delta_a_fit=delta_a_fit**2
delta_b_fit=delta_b_fit**2
delta_c_fit=delta_c_fit**2


 continueloop=.FALSE.



IF (delta_a_fit>(a_err**2)) THEN
 continueloop=.TRUE.
END IF
IF (delta_b_fit>(b_err**2)) THEN
 continueloop=.TRUE.
END IF
IF (delta_c_fit>(c_err**2)) THEN
 continueloop=.TRUE.
END IF



delta_a_fit=a_fit
delta_b_fit=b_fit
delta_c_fit=c_fit



IF (ee>maxloop) THEN
 continueloop=.FALSE.
END IF




IF (continueloop) THEN

!3-sigma-clipping
DO i=1,n_cg
IF (activef(i)) THEN

! WRITE(*,*) deviation(i),rms

IF ((deviation(i)**2.D0)>((3.D0*rms)**2.D0)) THEN
activef(i)=.FALSE.
END IF


END IF
END DO


! count galaxies still in use
galactive=0
DO i=1,n_cg
IF (activef(i)) THEN
galactive=galactive+1
END IF
END DO

WRITE(*,*) galactive

END IF

! end long loop
END DO



avlogmass=0.D0
DO i=1,n_g
avlogmass=avlogmass+LOG10(gal_mass(i))
END DO
avlogmass=avlogmass/DBLE(n_g)




OPEN(51,file='fit_lumass.txt')

WRITE(51,*) a_fit,b_fit,c_fit
WRITE(51,*) a_err,b_err,c_err
WRITE(51,*) rms,galactive

CLOSE(51)





OPEN(51,file='plot_fit_lumass.txt')

WRITE(51,*) 'f(x)=a_fit*x*x+b_fit*x+c_fit'
WRITE(51,*) 'a_fit=',a_fit
WRITE(51,*) 'b_fit=',b_fit
WRITE(51,*) 'c_fit=',c_fit

CLOSE(51)



DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

!maps for contoure plots
DO i=1,n_cg
IF (cg_fof_mass(i)>0.D0) THEN

help_x=NINT((LOG10(cg_lum(i))-7.7D0)*100.D0)
help_y=NINT((LOG10(cgnp(i)*mparticles)-10.D0)*80.D0)

IF ((help_x>0).AND.(help_x<400)) THEN
IF ((help_y>0).AND.(help_y<360)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF

END IF
END DO

OPEN(61,file='central_lum_mass.txt')
DO i=0,400
DO ii=0,360
help_x2=(DBLE(i)/100.D0)+7.7D0
help_y2=(DBLE(ii)/80.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)



DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

!maps for contoure plots
DO i=1,n_g


help_x=NINT((LOG10(g_lum(i))-7.7D0)*100.D0)
help_y=NINT((LOG10(gnp(i)*mparticles)-10.D0)*80.D0)

IF ((help_x>0).AND.(help_x<400)) THEN
IF ((help_y>0).AND.(help_y<360)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO

OPEN(61,file='gal_lum_mass.txt')
DO i=0,400
DO ii=0,360
help_x2=(DBLE(i)/100.D0)+7.7D0
help_y2=(DBLE(ii)/80.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)




DO i=1,n_g
gJmag(i)=ggDust(i)-(d_J*(ggDust(i)-grDust(i))+e_J*(grDust(i)-giDust(i))+f_J)
gKsmag(i)=ggDust(i)-(d_Ks*(ggDust(i)-grDust(i))+e_Ks*(grDust(i)-giDust(i))+f_Ks)
END DO

DO i=1,n_cg
 cgJmag(i)=cggDust(i)-(d_J*(cggDust(i)-cgrDust(i))+e_J*(cgrDust(i)-cgiDust(i))+f_J)
 cgKsmag(i)=cggDust(i)-(d_Ks*(cggDust(i)-cgrDust(i))+e_Ks*(cgrDust(i)-cgiDust(i))+f_Ks)
END DO






DO i=1,n_h
fof_unused(i)=.TRUE.
END DO

DO i=1,n_cg
g_lum(i)=10.D0**(-0.4D0*(gKsmag(i)-mag_sol_Ks))
END DO

DO i=1,n_cg
 cg_lum(i)=10.D0**(-0.4D0*(cgKsmag(i)-mag_sol_Ks)) 
 cg_fof_rad(i)=0.D0
 cg_fof_mass(i)=0.D0

DO ii=1,n_h
! IF (fof_unused(ii)) THEN
IF (h_fofid(ii)==cg_fofid(i)) THEN

! fof_unused(ii)=.FALSE.
 cg_fof_mass(i)=hm(ii)
 cg_fof_rad(i)=h_r(ii)
END IF
! END IF
END DO

END DO


help_count=0
OPEN(51,file='central_size_lum_2mass.txt')
DO i=1,n_cg
IF (cg_fof_mass(i)>0.D0) THEN
WRITE(51,*)  LOG10(cg_fof_rad(i)),LOG10(cg_lum(i)),LOG10(cg_fof_mass(i)),LOG10(cgnp(i)*mparticles)
ELSE
help_count=help_count+1
END IF
END DO
CLOSE(51)


WRITE(*,*) 'for',help_count,'of ',n_cg,'2mass central galaxies, we do not find FOF groups'

DO i=1,n_cg
IF (cg_fof_mass(i)>0.D0) THEN
activef(i)=.TRUE.
else
activef(i)=.false.
END IF
IF (LOG10(cgnp(i)*mparticles)>12.5D0) THEN
activef(i)=.false.
END IF
IF (LOG10(cgnp(i)*mparticles)<10.5D0) THEN
activef(i)=.false.
END IF
END DO


galactive=0
DO i=1,n_cg
IF (activef(i)) THEN
galactive=galactive+1
END IF
END DO
WRITE(*,*) galactive

 a_fit=0.25D0
 b_fit=-3.75D0
 c_fit=25.D0


 continueloop=.TRUE.
ee=0

delta_a_fit=0.D0
delta_b_fit=0.D0
delta_c_fit=0.D0


! do loop for opitmization
DO WHILE (continueloop)
ee=ee+1



!preparing matrizes
DO i=1,3 
Yvector(i)=0.D0
DO ii=1,3
Amatrix(i,ii)=0.D0
Bmatrix1(i,ii)=0.D0
Bmatrix2(i,ii)=0.D0
Bmatrix3(i,ii)=0.D0
END DO
END DO

!calculating basic matrix und vector
DO i=1,n_cg
IF (activef(i)) THEN

Amatrix(1,1)=Amatrix(1,1)+(LOG10(cg_lum(i))**4)
Amatrix(2,2)=Amatrix(2,2)+(LOG10(cg_lum(i))**2)
Amatrix(3,3)=Amatrix(3,3)+1.D0
Amatrix(1,2)=Amatrix(1,2)+(LOG10(cg_lum(i))**3)
Amatrix(2,1)=Amatrix(2,1)+(LOG10(cg_lum(i))**3)
Amatrix(1,3)=Amatrix(1,3)+(LOG10(cg_lum(i))**2)
Amatrix(3,1)=Amatrix(3,1)+(LOG10(cg_lum(i))**2)
Amatrix(2,3)=Amatrix(2,3)+LOG10(cg_lum(i))
Amatrix(3,2)=Amatrix(3,2)+LOG10(cg_lum(i))
Yvector(1)=Yvector(1)+LOG10(cgnp(i)*mparticles)*(LOG10(cg_lum(i))**2)
Yvector(2)=Yvector(2)+LOG10(cgnp(i)*mparticles)*LOG10(cg_lum(i))
Yvector(3)=Yvector(3)+LOG10(cgnp(i)*mparticles)

END IF
END DO



detA=Amatrix(1,1)*Amatrix(2,2)*Amatrix(3,3)
detA=detA+Amatrix(2,1)*Amatrix(3,2)*Amatrix(1,3)
detA=detA+Amatrix(1,2)*Amatrix(2,3)*Amatrix(3,1)
detA=detA-Amatrix(1,3)*Amatrix(2,2)*Amatrix(3,1)
detA=detA-Amatrix(2,1)*Amatrix(1,2)*Amatrix(3,3)
detA=detA-Amatrix(2,3)*Amatrix(3,2)*Amatrix(1,1)


!calculate B matrizes
DO i=1,3
DO ii=1,3
Bmatrix1(i,ii)=Amatrix(i,ii)
Bmatrix2(i,ii)=Amatrix(i,ii)
Bmatrix3(i,ii)=Amatrix(i,ii)
END DO
END DO

DO i=1,3
Bmatrix1(1,i)=Yvector(i)
Bmatrix2(2,i)=Yvector(i)
Bmatrix3(3,i)=Yvector(i)
END DO



! get determinat of Bmatrizes



detB1=Bmatrix1(1,1)*Bmatrix1(2,2)*Bmatrix1(3,3)
detB1=detB1+Bmatrix1(2,1)*Bmatrix1(3,2)*Bmatrix1(1,3)
detB1=detB1+Bmatrix1(1,2)*Bmatrix1(2,3)*Bmatrix1(3,1)
detB1=detB1-Bmatrix1(1,3)*Bmatrix1(2,2)*Bmatrix1(3,1)
detB1=detB1-Bmatrix1(2,1)*Bmatrix1(1,2)*Bmatrix1(3,3)
detB1=detB1-Bmatrix1(2,3)*Bmatrix1(3,2)*Bmatrix1(1,1)
detB2=Bmatrix2(1,1)*Bmatrix2(2,2)*Bmatrix2(3,3)
detB2=detB2+Bmatrix2(2,1)*Bmatrix2(3,2)*Bmatrix2(1,3)
detB2=detB2+Bmatrix2(1,2)*Bmatrix2(2,3)*Bmatrix2(3,1)
detB2=detB2-Bmatrix2(1,3)*Bmatrix2(2,2)*Bmatrix2(3,1)
detB2=detB2-Bmatrix2(2,1)*Bmatrix2(1,2)*Bmatrix2(3,3)
detB2=detB2-Bmatrix2(2,3)*Bmatrix2(3,2)*Bmatrix2(1,1)
detB3=Bmatrix3(1,1)*Bmatrix3(2,2)*Bmatrix3(3,3)
detB3=detB3+Bmatrix3(2,1)*Bmatrix3(3,2)*Bmatrix3(1,3)
detB3=detB3+Bmatrix3(1,2)*Bmatrix3(2,3)*Bmatrix3(3,1)
detB3=detB3-Bmatrix3(1,3)*Bmatrix3(2,2)*Bmatrix3(3,1)
detB3=detB3-Bmatrix3(2,1)*Bmatrix3(1,2)*Bmatrix3(3,3)
detB3=detB3-Bmatrix3(2,3)*Bmatrix3(3,2)*Bmatrix3(1,1)


 a_fit=detB1/detA
 b_fit=detB2/detA
 c_fit=detB3/detA


!calculate right side of fundamental plane
DO i=1,n_cg
IF (activef(i)) THEN


otherside(i)=a_fit*((LOG10(cg_lum(i)))**2)+b_fit*LOG10(cg_lum(i))+c_fit


END IF
END DO

! propagation of error
!root mean square
rms=0.D0

DO i=1,n_cg
IF (activef(i)) THEN



deviation(i)=otherside(i)-LOG10(cgnp(i)*mparticles)

rms=rms+((deviation(i))**2)

END IF
END DO



rms=SQRT(rms/DBLE(galactive))


! error of fitting parameters


 a_err=(Amatrix(2,2)*Amatrix(3,3)-Amatrix(2,3)*Amatrix(3,2))/detA
 a_err=SQRT(a_err)
 b_err=(Amatrix(1,1)*Amatrix(3,3)-Amatrix(1,3)*Amatrix(3,1))/detA
 b_err=SQRT(b_err)
 c_err=(Amatrix(1,1)*Amatrix(2,2)-Amatrix(1,2)*Amatrix(2,1))/detA
 c_err=SQRT(c_err)


!see if another loop is still nessecary


delta_a_fit=delta_a_fit-a_fit
delta_b_fit=delta_b_fit-b_fit
delta_c_fit=delta_c_fit-c_fit
delta_a_fit=delta_a_fit*100.D0
delta_b_fit=delta_b_fit*100.D0
delta_c_fit=delta_c_fit*100.D0
delta_a_fit=delta_a_fit**2
delta_b_fit=delta_b_fit**2
delta_c_fit=delta_c_fit**2


 continueloop=.FALSE.



IF (delta_a_fit>(a_err**2)) THEN
 continueloop=.TRUE.
END IF
IF (delta_b_fit>(b_err**2)) THEN
 continueloop=.TRUE.
END IF
IF (delta_c_fit>(c_err**2)) THEN
 continueloop=.TRUE.
END IF



delta_a_fit=a_fit
delta_b_fit=b_fit
delta_c_fit=c_fit



IF (ee>maxloop) THEN
 continueloop=.FALSE.
END IF




IF (continueloop) THEN

!3-sigma-clipping
DO i=1,n_cg
IF (activef(i)) THEN

! WRITE(*,*) deviation(i),rms

IF ((deviation(i)**2.D0)>((3.D0*rms)**2.D0)) THEN
activef(i)=.FALSE.
END IF


END IF
END DO


! count galaxies still in use
galactive=0
DO i=1,n_cg
IF (activef(i)) THEN
galactive=galactive+1
END IF
END DO

WRITE(*,*) galactive

END IF

! end long loop
END DO



avlogmass=0.D0
DO i=1,n_g
avlogmass=avlogmass+LOG10(gal_mass(i))
END DO
avlogmass=avlogmass/DBLE(n_g)




OPEN(51,file='fit_lumass_2mass.txt')

WRITE(51,*) a_fit,b_fit,c_fit
WRITE(51,*) a_err,b_err,c_err
WRITE(51,*) rms,galactive

CLOSE(51)





OPEN(51,file='plot_fit_lumass_2mass.txt')

WRITE(51,*) 'f(x)=a_fit*x*x+b_fit*x+c_fit'
WRITE(51,*) 'a_fit=',a_fit
WRITE(51,*) 'b_fit=',b_fit
WRITE(51,*) 'c_fit=',c_fit

CLOSE(51)






DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

!maps for contoure plots
DO i=1,n_cg
IF (cg_fof_mass(i)>0.D0) THEN

help_x=NINT((LOG10(cg_lum(i))-7.7D0)*100.D0)
help_y=NINT((LOG10(cgnp(i)*mparticles)-10.D0)*80.D0)

IF ((help_x>0).AND.(help_x<400)) THEN
IF ((help_y>0).AND.(help_y<360)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF

END IF
END DO

OPEN(61,file='central_lum_mass_2mass.txt')
DO i=0,400
DO ii=0,360
help_x2=(DBLE(i)/100.D0)+7.7D0
help_y2=(DBLE(ii)/80.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)



DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

!maps for contoure plots
DO i=1,n_g


help_x=NINT((LOG10(g_lum(i))-7.7D0)*100.D0)
help_y=NINT((LOG10(gnp(i)*mparticles)-10.D0)*80.D0)

IF ((help_x>0).AND.(help_x<400)) THEN
IF ((help_y>0).AND.(help_y<360)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO

OPEN(61,file='gal_lum_mass_2mass.txt')
DO i=0,400
DO ii=0,360
help_x2=(DBLE(i)/100.D0)+7.7D0
help_y2=(DBLE(ii)/80.D0)+10.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

















OPEN(51,file='globalparameters.txt')

WRITE(51,*) avnextdist
WRITE(51,*) avdist
WRITE(51,*) av_velrad
WRITE(51,*) avlogmass

CLOSE(51)
WRITE(*,*) 'mystuffdone'

IF (doshells==4) THEN
OPEN(51,file='globalparameters_SDSS.txt')

WRITE(51,*) avnextdist
WRITE(51,*) avdist
WRITE(51,*) av_velrad
WRITE(51,*) avlogmass

CLOSE(51)



n_2mrs=0
DO i=1,n_g
IF (gKsmag(i)<-18.D0) THEN
n_2mrs=n_2mrs+1
END IF
END DO

avnextdist=0.D0
avdist=0.D0
DO i=1,n_g
IF (gKsmag(i)<-18.D0) THEN
nextdist=100.D0
DO ii=1,n_g
IF (gKsmag(i)<-18.D0) THEN
IF (ii.NE.i) THEN
helpdist=((gx(i)-gx(ii))**2)+((gy(i)-gy(ii))**2)+((gz(i)-gz(ii))**2)
avdist=avdist+SQRT(helpdist)
IF (helpdist<(nextdist**2)) THEN
nextdist=SQRT(helpdist)
END IF
END IF
END IF
END DO
avnextdist=avnextdist+nextdist
END IF
END DO
avnextdist=avnextdist/DBLE(n_2mrs)
avdist=avdist/(DBLE(n_2mrs)*DBLE(n_2mrs-1))

av_velrad=0.D0
DO i=1,n_g
IF (gKsmag(i)<-18.D0) THEN
av_velrad=av_velrad+(gvelX(i)**2)+(gvelY(i)**2)+(gvelZ(i)**2)
END IF
END DO
av_velrad=SQRT(av_velrad/(3.D0*DBLE(n_2mrs)))





OPEN(51,file='globalparameters_2MRS.txt')

WRITE(51,*) avnextdist
WRITE(51,*) avdist
WRITE(51,*) av_velrad
WRITE(51,*) avlogmass

CLOSE(51)





END IF



 IF (doshells<3) THEN
 
 
DO i=1,n_p
near_central_10(i)=.TRUE.
END DO


DO i=1,n_cg 
DO ii=1,n_p
IF (near_central_10(ii)) THEN
dist=((cgx(i)-px(ii))**2)+((cgy(i)-py(ii))**2)+((cgz(i)-pz(ii))**2)

IF (dist<rhm10(i)) THEN
near_central_10(ii)=.FALSE.
END IF

END IF
END DO
END DO


help_count=n_p
DO i=1,n_p
IF (near_central_10(i)) THEN
help_count=help_count-1
END IF
END DO

WRITE(*,*) DBLE(help_count)/part_exp*100.D0,'% of all particles are near the central galaxies'




END IF

IF (doshells==1) THEN


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
!WRITE(*,*) i


DO ii=1,n_p

dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)

IF (inside(ii).EQV..FALSE.) THEN
IF (dist<hr2(i)) THEN
inside(ii)=.TRUE.
END IF
END IF

IF (in_rhm10(ii).EQV..FALSE.) THEN
IF (dist<rhm10(i)) THEN
in_rhm10(ii)=.TRUE.
END IF
END IF

IF (in_fi(ii).EQV..FALSE.) THEN
IF (dist<r_fi(i)) THEN
in_fi(ii)=.TRUE.
END IF
END IF

IF (in_critdensity(ii).EQV..FALSE.) THEN
IF (dist<r_crit(i)) THEN
in_critdensity(ii)=.TRUE.
END IF
END IF

IF (in_100(ii).EQV..FALSE.) THEN
IF (dist<r_100(i)) THEN
in_100(ii)=.TRUE.
END IF
END IF

IF (in_200(ii).EQV..FALSE.) THEN
IF (dist<r_200(i)) THEN
in_200(ii)=.TRUE.
END IF
END IF


END DO

DO ii=1,n_rand



 DO i1=1,3
 DO i2=1,3
 DO i3=1,3
 
 


dist=((hx(i)-rand_x(ii)-cubeshift_x(i1,i2,i3))**2)+&
((hy(i)-rand_y(ii)-cubeshift_y(i1,i2,i3))**2)+((hz(i)-rand_z(ii)-cubeshift_z(i1,i2,i3))**2)




IF (rand_inside(ii).EQV..FALSE.) THEN
IF (dist<hr2(i)) THEN
rand_inside(ii)=.TRUE.
END IF
END IF

IF (rand_in_rhm10(ii).EQV..FALSE.) THEN
IF (dist<rhm10(i)) THEN
rand_in_rhm10(ii)=.TRUE.
END IF
END IF

IF (rand_in_fi(ii).EQV..FALSE.) THEN
IF (dist<r_fi(i)) THEN
rand_in_fi(ii)=.TRUE.
END IF
END IF

IF (rand_in_critdensity(ii).EQV..FALSE.) THEN
IF (dist<r_crit(i)) THEN
rand_in_critdensity(ii)=.TRUE.
END IF
END IF

IF (rand_in_100(ii).EQV..FALSE.) THEN
IF (dist<r_100(i)) THEN
rand_in_100(ii)=.TRUE.
END IF
END IF

IF (rand_in_200(ii).EQV..FALSE.) THEN
IF (dist<r_200(i)) THEN
rand_in_200(ii)=.TRUE.
END IF
END IF


 END DO
 END DO
 END DO




END DO

END DO




n_inside=0
n_fi=0
n_critdens=0
n_200=0
n_100=0
n_rhm10=0


DO i=1,n_p


IF (inside(i)) THEN
n_inside=n_inside+1
END IF


IF (in_rhm10(i)) THEN
n_rhm10=n_rhm10+1
END IF


IF (in_fi(i)) THEN
n_fi=n_fi+1
END IF


IF (in_critdensity(i)) THEN
n_critdens=n_critdens+1
END IF


IF (in_100(i)) THEN
n_100=n_100+1
END IF


IF (in_200(i)) THEN
n_200=n_200+1
END IF




END DO



nrand_inside=0
nrand_fi=0
nrand_critdens=0
nrand_200=0
nrand_100=0
nrand_rhm10=0
nrand_critdens_iter=0
nrand_fi_iter=0



DO i=1,n_rand


IF (rand_inside(i)) THEN
nrand_inside=nrand_inside+1
END IF


IF (rand_in_rhm10(i)) THEN
nrand_rhm10=nrand_rhm10+1
END IF


IF (rand_in_fi(i)) THEN
nrand_fi=nrand_fi+1
END IF


IF (rand_in_critdensity(i)) THEN
nrand_critdens=nrand_critdens+1
END IF


IF (rand_in_100(i)) THEN
nrand_100=nrand_100+1
END IF


IF (rand_in_200(i)) THEN
nrand_200=nrand_200+1
END IF



END DO


DO i=1,n_h
M_crit_rad(i)=0.D0
M_fi_rad(i)=0.D0
M_100_rad(i)=0.D0
M_200_rad(i)=0.D0
Mr_rad(i)=0.D0
Mr10_rad(i)=0.D0


DO ii=1,n_p

dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)


IF (dist<r_crit(i)) THEN
M_crit_rad(i)=M_crit_rad(i)+mparticles
END IF


IF (dist<r_fi(i)) THEN
M_fi_rad(i)=M_fi_rad(i)+mparticles
END IF


IF (dist<hr2(i)) THEN
Mr_rad(i)=Mr_rad(i)+mparticles
END IF


IF (dist<rhm10(i)) THEN
Mr10_rad(i)=Mr10_rad(i)+mparticles
END IF


IF (dist<r_100(i)) THEN
M_100_rad(i)=M_100_rad(i)+mparticles
END IF


IF (dist<r_200(i)) THEN
M_200_rad(i)=M_200_rad(i)+mparticles
END IF






END DO

END DO


OPEN(51,file='Mass_radii.txt')
DO i=1,n_h
WRITE(51,*) LOG10(fof_group_totlum(i)),LOG10(hm(i)),LOG10(Mr_rad(i)),LOG10(Mr10_rad(i)),&
LOG10(M_crit_rad(i)),LOG10(M_fi_rad(i)),LOG10(M_200_rad(i)),LOG10(M_100_rad(i)),LOG10(DBLE(fof_group_members(i)))
END DO
CLOSE(51)

OPEN(51,file='Mass_radii_m1.txt')
DO i=1,n_h
IF (fof_group_members(i)==1) THEN
WRITE(51,*) LOG10(fof_group_totlum(i)),LOG10(hm(i)),LOG10(Mr_rad(i)),LOG10(Mr10_rad(i)),&
LOG10(M_crit_rad(i)),LOG10(M_fi_rad(i)),LOG10(M_200_rad(i)),LOG10(M_100_rad(i))
END IF
END DO
CLOSE(51)

OPEN(51,file='Mass_radii_m2-5.txt')
DO i=1,n_h
IF ((fof_group_members(i)>1).AND.(fof_group_members(i)<6)) THEN
WRITE(51,*) LOG10(fof_group_totlum(i)),LOG10(hm(i)),LOG10(Mr_rad(i)),LOG10(Mr10_rad(i)),&
LOG10(M_crit_rad(i)),LOG10(M_fi_rad(i)),LOG10(M_200_rad(i)),LOG10(M_100_rad(i))
END IF
END DO
CLOSE(51)

OPEN(51,file='Mass_radii_m6-10.txt')
DO i=1,n_h
IF ((fof_group_members(i)>5).AND.(fof_group_members(i)<11)) THEN
WRITE(51,*) LOG10(fof_group_totlum(i)),LOG10(hm(i)),LOG10(Mr_rad(i)),LOG10(Mr10_rad(i)),&
LOG10(M_crit_rad(i)),LOG10(M_fi_rad(i)),LOG10(M_200_rad(i)),LOG10(M_100_rad(i))
END IF
END DO
CLOSE(51)

OPEN(51,file='Mass_radii_m11-99.txt')
DO i=1,n_h
IF ((fof_group_members(i)>10).AND.(fof_group_members(i)<100)) THEN
WRITE(51,*) LOG10(fof_group_totlum(i)),LOG10(hm(i)),LOG10(Mr_rad(i)),LOG10(Mr10_rad(i)),&
LOG10(M_crit_rad(i)),LOG10(M_fi_rad(i)),LOG10(M_200_rad(i)),LOG10(M_100_rad(i))
END IF
END DO
CLOSE(51)

OPEN(51,file='Mass_radii_m100+.txt')
DO i=1,n_h
IF (fof_group_members(i)>100) THEN
WRITE(51,*) LOG10(fof_group_totlum(i)),LOG10(hm(i)),LOG10(Mr_rad(i)),LOG10(Mr10_rad(i)),&
LOG10(M_crit_rad(i)),LOG10(M_fi_rad(i)),LOG10(M_200_rad(i)),LOG10(M_100_rad(i))
END IF
END DO
CLOSE(51)
! WRITE(*,*)


! read file
OPEN(50,file='results.txt')
WRITE(50,*) part_tot/part_exp*100.D0,'% of all particles are in halos >=20 particles'
WRITE(50,*) '-------------------------------------------'
WRITE(50,*) 'radius particles percentage'
WRITE(50,*) 'half light:',n_inside,(100.D0*DBLE(n_inside)/n_p)
WRITE(50,*) '10 times half light:',n_rhm10,(100.D0*DBLE(n_rhm10)/n_p)
WRITE(50,*) 'finite infinity:',n_fi,(100.D0*DBLE(n_fi)/n_p)
WRITE(50,*) 'critical density:',n_critdens,(100.D0*DBLE(n_critdens)/n_p)
WRITE(50,*) '100 times critical density:',n_100,(100.D0*DBLE(n_100)/n_p)
WRITE(50,*) '200 times critical density:',n_200,(100.D0*DBLE(n_200)/n_p)




WRITE(*,*) 'half light:',n_inside,(100.D0*DBLE(n_inside)/n_p)
WRITE(*,*) '10 times half light:',n_rhm10,(100.D0*DBLE(n_rhm10)/n_p)
WRITE(*,*) 'finite infinity:',n_fi,(100.D0*DBLE(n_fi)/n_p)
WRITE(*,*) 'critical density:',n_critdens,(100.D0*DBLE(n_critdens)/n_p)
WRITE(*,*) '100 times critical density:',n_100,(100.D0*DBLE(n_100)/n_p)
WRITE(*,*) '200 times critical density:',n_200,(100.D0*DBLE(n_200)/n_p)


OPEN(51,file='firstrun.txt')
DO i=1,n_h
WRITE(51,*) hx(i),hy(i),hz(i),hm(i),SQRT(r_fi(i)),SQRT(r_crit(i)),SQRT(r_100(i)),SQRT(r_200(i)),SQRT(rhm10(i)),h_r(i)
END DO
CLOSE(51)

DO i=1,n_h
active(i)=.TRUE. 
r_fi(i)=SQRT(r_fi(i))
hmfi(i)=hm(i)
END DO



IF (no_iter==0) THEN



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


mnew=hmfi(i)+hmfi(ii)
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
r_fi(i)=r_fi(i)**2
END IF
END DO


WRITE(*,*) ii,'halos will be used after the filtering'
WRITE(50,*) ii,'halos will be used after the filtering'
















DO i=1,n_p
in_fi(i)=.FALSE.
END DO


DO i=1,n_h
IF (active(i)) THEN
!WRITE(*,*) i
DO ii=1,n_p
IF (in_fi(ii).EQV..FALSE.) THEN
dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)

IF (dist<r_fi(i)) THEN
in_fi(ii)=.TRUE.
END IF
END IF

END DO
END IF
END DO



n_fi=0



DO i=1,n_p

IF (in_fi(i)) THEN
n_fi=n_fi+1
END IF


END DO




WRITE(*,*) 'inside finite infinity:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
WRITE(50,*) 'inside finite infinity:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
anpfi=n_fi

DO i=1,n_h
IF (active(i)) THEN
r_fi(i)=SQRT(r_fi(i))
hmfi(i)=hmfi(i)*(DBLE(n_fi)/DBLE(n_p))/(npartalt/part_exp)
r_fi(i)=(3.D0/(4.D0*PI)*hmfi(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)
END IF
END DO

npartalt=n_fi



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


! dist=((hx(i)-rand_x(ii))**2)+((hy(i)-rand_y(ii))**2)+((hz(i)-rand_z(ii))**2)
! 



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

DO i=1,n_rand

IF (rand_in_fi_iter(i)) THEN
nrand_fi_iter=nrand_fi_iter+1
END IF


END DO


OPEN(51,file='final_fi.txt')
DO i=1,n_h
IF (active(i)) THEN
WRITE(51,*) hx(i),hy(i),hz(i),hmfi(i),r_fi(i)
END IF
END DO
CLOSE(51)













DO i=1,n_h
active(i)=.TRUE. 
r_crit(i)=SQRT(r_crit(i))
hmfi(i)=hm(i)
END DO







anpfi=0
npartalt=part_tot




uuu=0
DO WHILE (uuu<10)
uuu=uuu+1
WRITE(*,*) '-------------------------------------------'
WRITE(50,*) '-------------------------------------------'
WRITE(*,*) uuu,'iteration for critical'
WRITE(50,*) uuu,'iteration for critical'
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
IF ((d+r_crit(ii))<r_crit(i)) THEN
active(ii)=.FALSE.


mnew=hmfi(i)+hmfi(ii)
hx(i)=(hx(i)*hmfi(i)+hx(ii)*hmfi(ii))/mnew
hy(i)=(hy(i)*hmfi(i)+hy(ii)*hmfi(ii))/mnew
hz(i)=(hz(i)*hmfi(i)+hz(ii)*hmfi(ii))/mnew
hmfi(i)=mnew



r_crit(i)=(3.D0/(4.D0*PI)*hmfi(i)/(rho_crit))**(1.D0/3.D0)

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
r_crit(i)=r_crit(i)**2
END IF
END DO


WRITE(*,*) ii,'halos will be used after the filtering'
WRITE(50,*) ii,'halos will be used after the filtering'
















DO i=1,n_p
in_fi(i)=.FALSE.
END DO


DO i=1,n_h
IF (active(i)) THEN
!WRITE(*,*) i
DO ii=1,n_p
IF (in_fi(ii).EQV..FALSE.) THEN
dist=((hx(i)-px(ii))**2)+((hy(i)-py(ii))**2)+((hz(i)-pz(ii))**2)

IF (dist<r_crit(i)) THEN
in_fi(ii)=.TRUE.
END IF
END IF

END DO
END IF
END DO



n_fi=0



DO i=1,n_p

IF (in_fi(i)) THEN
n_fi=n_fi+1
END IF


END DO








WRITE(*,*) 'inside critical:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
WRITE(50,*) 'inside critical:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
anpfi=n_fi

DO i=1,n_h
IF (active(i)) THEN
r_crit(i)=SQRT(r_crit(i))
hmfi(i)=hmfi(i)*(DBLE(n_fi)/DBLE(n_p))/(npartalt/part_exp)
r_crit(i)=(3.D0/(4.D0*PI)*hmfi(i)/(rho_crit))**(1.D0/3.D0)

END IF
END DO

npartalt=n_fi



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

! dist=((hx(i)-rand_x(ii))**2)+((hy(i)-rand_y(ii))**2)+((hz(i)-rand_z(ii))**2)




IF (rand_in_critdensity_iter(ii).EQV..FALSE.) THEN
IF (dist<(r_crit(i)**2)) THEN
rand_in_critdensity_iter(ii)=.TRUE.
END IF
END IF

END DO
END DO
END DO

END DO
END IF
END DO

DO i=1,n_rand

IF (rand_in_critdensity_iter(i)) THEN
nrand_critdens_iter=nrand_critdens_iter+1
END IF


END DO



OPEN(51,file='final_critical.txt')
DO i=1,n_h
IF (active(i)) THEN
WRITE(51,*) hx(i),hy(i),hz(i),hmfi(i),r_crit(i)
END IF
END DO
CLOSE(51)


OPEN(51,file='volume_percentage.txt')
WRITE(51,*) 'in group:',(DBLE(nrand_inside)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in fi:',(DBLE(nrand_fi)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in crit_dens:',(DBLE(nrand_critdens)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in crit_200:',(DBLE(nrand_200)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in crit 100:',(DBLE(nrand_100)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in 10r group:',(DBLE(nrand_rhm10)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in crit_dens iter:',(DBLE(nrand_critdens_iter)/DBLE(n_rand)*100.D0),'%'
WRITE(51,*) 'in fi iter:',(DBLE(nrand_fi_iter)/DBLE(n_rand)*100.D0),'%'
CLOSE(51)




CLOSE(50)


END IF

END IF



IF (doshells==2) THEN




DO i=1,n_g
!WRITE(*,*) i
DO ii=1,n_p

dist=((gx(i)-px(ii))**2)+((gy(i)-py(ii))**2)+((gz(i)-pz(ii))**2)



IF (gal_in_fi(ii).EQV..FALSE.) THEN
IF (dist<gal_rad_fi(i)) THEN
gal_in_fi(ii)=.TRUE.
END IF
END IF

IF (gal_in_critdensity(ii).EQV..FALSE.) THEN
IF (dist<gal_rad_crit(i)) THEN
gal_in_critdensity(ii)=.TRUE.
END IF
END IF

IF (gal_in_100(ii).EQV..FALSE.) THEN
IF (dist<gal_rad_100(i)) THEN
gal_in_100(ii)=.TRUE.
END IF
END IF

IF (gal_in_200(ii).EQV..FALSE.) THEN
IF (dist<gal_rad_200(i)) THEN
gal_in_200(ii)=.TRUE.
END IF
END IF




END DO
END DO




n_inside=0
n_fi=0
n_critdens=0
n_200=0
n_100=0
n_rhm10=0


DO i=1,n_p



IF (gal_in_fi(i)) THEN
n_fi=n_fi+1
END IF


IF (gal_in_critdensity(i)) THEN
n_critdens=n_critdens+1
END IF


IF (gal_in_100(i)) THEN
n_100=n_100+1
END IF


IF (gal_in_200(i)) THEN
n_200=n_200+1
END IF




END DO



! WRITE(*,*)


! read file
OPEN(50,file='results_gal.txt')
! WRITE(50,*) part_tot/part_exp*100.D0,'% of all particles are near galaxies'
WRITE(50,*) '-------------------------------------------'
WRITE(50,*) 'radius particles percentage'
WRITE(50,*) 'finite infinity:',n_fi,(100.D0*DBLE(n_fi)/n_p)
WRITE(50,*) 'critical density:',n_critdens,(100.D0*DBLE(n_critdens)/n_p)
WRITE(50,*) '100 times critical density:',n_100,(100.D0*DBLE(n_100)/n_p)
WRITE(50,*) '200 times critical density:',n_200,(100.D0*DBLE(n_200)/n_p)


WRITE(*,*) 'finite infinity:',n_fi,(100.D0*DBLE(n_fi)/n_p)
WRITE(*,*) 'critical density:',n_critdens,(100.D0*DBLE(n_critdens)/n_p)
WRITE(*,*) '100 times critical density:',n_100,(100.D0*DBLE(n_100)/n_p)
WRITE(*,*) '200 times critical density:',n_200,(100.D0*DBLE(n_200)/n_p)


OPEN(51,file='firstrun_gal.txt')
DO i=1,n_g
WRITE(51,*) gx(i),gy(i),gz(i),gal_mass(i),SQRT(gal_rad_fi(i)),&
SQRT(gal_rad_crit(i)),SQRT(gal_rad_100(i)),SQRT(gal_rad_200(i))
END DO
CLOSE(51)


DO i=1,n_g
activeg(i)=.TRUE. 
gal_rad_fi(i)=SQRT(gal_rad_fi(i))
gmfi(i)=gal_mass(i)
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
DO i=1,n_g
IF (activeg(i)) THEN
!WRITE(*,*) i

ii=0
DO WHILE (ii<n_g)
ii=ii+1
IF (activeg(ii)) THEN


IF (i/=ii) THEN
d=(gx(i)-gx(ii))**2+(gy(i)-gy(ii))**2+(gz(i)-gz(ii))**2
d=SQRT(d)
IF ((d+gal_rad_fi(ii))<gal_rad_fi(i)) THEN
activeg(ii)=.FALSE.


mnew=gmfi(i)+gmfi(ii)
gx(i)=(gx(i)*gmfi(i)+gx(ii)*gmfi(ii))/mnew
gy(i)=(gy(i)*gmfi(i)+gy(ii)*gmfi(ii))/mnew
gz(i)=(gz(i)*gmfi(i)+gz(ii)*gmfi(ii))/mnew
gmfi(i)=mnew



gal_rad_fi(i)=(3.D0/(4.D0*PI)*gmfi(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)

ii=0
END IF
END IF
END IF

END DO

END IF
END DO


ii=0
DO i=1,n_g
IF (activeg(i)) THEN
ii=ii+1
gal_rad_fi(i)=gal_rad_fi(i)**2
END IF
END DO


WRITE(*,*) ii,'galaxies will be used after the filtering'
WRITE(50,*) ii,'galaxies will be used after the filtering'
















DO i=1,n_p
gal_in_fi(i)=.FALSE.
END DO


DO i=1,n_h
IF (activeg(i)) THEN
!WRITE(*,*) i
DO ii=1,n_p
IF (gal_in_fi(ii).EQV..FALSE.) THEN
dist=((gx(i)-px(ii))**2)+((gy(i)-py(ii))**2)+((gz(i)-pz(ii))**2)

IF (dist<gal_rad_fi(i)) THEN
gal_in_fi(ii)=.TRUE.
END IF
END IF

END DO
END IF
END DO



n_fi=0



DO i=1,n_p

IF (gal_in_fi(i)) THEN
n_fi=n_fi+1
END IF


END DO


WRITE(*,*) 'inside finite infinity:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
WRITE(50,*) 'inside finite infinity:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
anpfi=n_fi

DO i=1,n_g
IF (activeg(i)) THEN
gal_rad_fi(i)=SQRT(gal_rad_fi(i))
gmfi(i)=gmfi(i)*(DBLE(n_fi)/DBLE(n_p))/(npartalt/part_exp)
gal_rad_fi(i)=(3.D0/(4.D0*PI)*gmfi(i)/(rho_crit*ts_modificator))**(1.D0/3.D0)
END IF
END DO

npartalt=n_fi



END DO



WRITE(*,*) '-------------------------------------------'
WRITE(50,*) '-------------------------------------------'





OPEN(51,file='final_fi_gal.txt')
DO i=1,n_g
IF (activeg(i)) THEN
WRITE(51,*) gx(i),gy(i),gz(i),gmfi(i),gal_rad_fi(i)
END IF
END DO
CLOSE(51)













DO i=1,n_g
activeg(i)=.TRUE. 
gal_rad_crit(i)=SQRT(gal_rad_crit(i))
gmfi(i)=gal_mass(i)
END DO







anpfi=0
npartalt=part_tot




uuu=0
DO WHILE (uuu<10)
uuu=uuu+1
WRITE(*,*) '-------------------------------------------'
WRITE(50,*) '-------------------------------------------'
WRITE(*,*) uuu,'iteration for critical'
WRITE(50,*) uuu,'iteration for critical'
IF (uuu>2) THEN
IF ((DBLE(anpfi)/DBLE(n_fi))<1.001D0) THEN
IF ((DBLE(anpfi)/DBLE(n_fi))>0.999D0) THEN
uuu=11
END IF
END IF
END IF

! filter fi regions which are inside of other fi regions
DO i=1,n_g
IF (activeg(i)) THEN
!WRITE(*,*) i

ii=0
DO WHILE (ii<n_g)
ii=ii+1
IF (activeg(ii)) THEN


IF (i/=ii) THEN
d=(gx(i)-gx(ii))**2+(gy(i)-gy(ii))**2+(gz(i)-gz(ii))**2
d=SQRT(d)
IF ((d+gal_rad_crit(ii))<gal_rad_crit(i)) THEN
activeg(ii)=.FALSE.


mnew=gmfi(i)+gmfi(ii)
gx(i)=(gx(i)*gmfi(i)+gx(ii)*gmfi(ii))/mnew
gy(i)=(gy(i)*gmfi(i)+gy(ii)*gmfi(ii))/mnew
gz(i)=(gz(i)*gmfi(i)+gz(ii)*gmfi(ii))/mnew
gmfi(i)=mnew



gal_rad_crit(i)=(3.D0/(4.D0*PI)*gmfi(i)/(rho_crit))**(1.D0/3.D0)

ii=0
END IF
END IF
END IF

END DO

END IF
END DO


ii=0
DO i=1,n_h
IF (activeg(i)) THEN
ii=ii+1
gal_rad_crit(i)=gal_rad_crit(i)**2
END IF
END DO


WRITE(*,*) ii,'galaxies will be used after the filtering'
WRITE(50,*) ii,'galaxies will be used after the filtering'
















DO i=1,n_p
gal_in_fi(i)=.FALSE.
END DO


DO i=1,n_g
IF (activeg(i)) THEN
!WRITE(*,*) i
DO ii=1,n_p
IF (gal_in_fi(ii).EQV..FALSE.) THEN
dist=((gx(i)-px(ii))**2)+((gy(i)-py(ii))**2)+((gz(i)-pz(ii))**2)

IF (dist<gal_rad_crit(i)) THEN
gal_in_fi(ii)=.TRUE.
END IF
END IF

END DO
END IF
END DO



n_fi=0



DO i=1,n_p

IF (gal_in_fi(i)) THEN
n_fi=n_fi+1
END IF


END DO


WRITE(*,*) 'inside critical:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
WRITE(50,*) 'inside critical:',n_fi,(100.D0*DBLE(n_fi)/DBLE(n_p))
anpfi=n_fi

DO i=1,n_g
IF (activeg(i)) THEN
gal_rad_crit(i)=SQRT(gal_rad_crit(i))
gmfi(i)=gmfi(i)*(DBLE(n_fi)/DBLE(n_p))/(npartalt/part_exp)
gal_rad_crit(i)=(3.D0/(4.D0*PI)*gmfi(i)/(rho_crit))**(1.D0/3.D0)

END IF
END DO

npartalt=n_fi



END DO



WRITE(*,*) '-------------------------------------------'
WRITE(50,*) '-------------------------------------------'





OPEN(51,file='final_critical_gal.txt')
DO i=1,n_g
IF (activeg(i)) THEN
WRITE(51,*) gx(i),gy(i),gz(i),gmfi(i),gal_rad_crit(i)
END IF
END DO
CLOSE(51)









CLOSE(50)




END IF


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


