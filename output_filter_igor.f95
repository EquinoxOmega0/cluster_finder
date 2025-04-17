PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE

! include 'mpif.h'

real :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,n,i,ii,good_mag,n_2mass,best_2mrs
real, allocatable :: ra(:),dec(:),redshift(:),MagU(:),MagB(:),MagV(:),MagRx(:),MagK(:),Magg(:),Magr(:) 
integer, allocatable :: numbergal_2mrs(:),internalid_2mrs(:)
real, allocatable :: twomass_ra(:),twomass_dec(:),c_redshift_2mrs(:)
logical, allocatable :: active(:) 
real :: dummy_var,ang_min,z_min,angular_sep,delta_z,dec1,dec2,ra1,ra2,asold
character(200), allocatable :: mrs_id(:),id_name(:)
character(40) :: str1,str2,pref1,pref2,signum



WRITE(*,*) '============================================================'
WRITE(*,*) '    programme SPECIAL OUTPUT started'
WRITE(*,*) '============================================================'




! define constants
PI=ACOS(-1.D0)



! get length of file
OPEN(50,file='simbad_all.txt')
io_err=0
n=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n=n+1
END DO
 CLOSE(50)
n=n-2 ! remove header


allocate(ra(1:n))
allocate(dec(1:n))
allocate(redshift(1:n))
allocate(MagU(1:n))
allocate(MagB(1:n))
allocate(MagV(1:n))
allocate(MagRx(1:n))
allocate(MagK(1:n))
allocate(Magg(1:n))
allocate(Magr(1:n))
allocate(active(1:n))
allocate(internalid_2mrs(1:n))
allocate(id_name(1:n))
 

! read file
OPEN(50,file='simbad_all.txt')
READ(50,*)
DO i=1,n
! WRITE(*,*) i
READ(50,*) dummy_var,ra(i),dec(i),redshift(i),MagU(i),MagB(i),MagV(i),MagRx(i),MagK(i),Magg(i),Magr(i)

active(i)=.TRUE.

END DO

 CLOSE(50)


DO i=1,n
good_mag=0

IF (MagU(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagB(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagV(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagRx(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagK(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (Magg(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (Magr(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (good_mag<3) THEN
active(i)=.FALSE.
END IF

! IF (redshift(i)<0.D0) THEN
! active(i)=.FALSE.
! END IF

END DO





! get length of file
OPEN(50,file='igor/2MRS_igor.txt')
io_err=0
n_2mass=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mass=n_2mass+1
END DO
 CLOSE(50)
n_2mass=n_2mass-1 


allocate(numbergal_2mrs(1:n_2mass))
allocate(mrs_id(1:n_2mass))
allocate(twomass_ra(1:n_2mass))
allocate(twomass_dec(1:n_2mass))
allocate(c_redshift_2mrs(1:n_2mass))

OPEN(50,file='igor/2MRS_igor.txt')

! WRITE(50,*) "internal_2MRS_ID 2MRS_ID RA DEC Z_group K_t H_t J_t E_BV"
DO i=1,n_2mass 
READ(50,*) numbergal_2mrs(i),mrs_id(i),twomass_ra(i),twomass_dec(i),&
 c_redshift_2mrs(i),dummy_var,dummy_var,dummy_var,dummy_var

END DO

 CLOSE(50)

 
 
 
ang_min=1.D0/60.D0
z_min=1.D-2


DO i=1,n
IF (active(i)) THEN

best_2mrs=0
asold=1000.D0
DO ii=1,n_2mass

delta_z=c_redshift_2mrs(ii)-redshift(i)
IF (delta_z<0.D0) THEN
delta_z=-delta_z
END IF

IF (delta_z<z_min) THEN

ra1=ra(i)*PI/180.D0
ra2=twomass_ra(ii)*PI/180.D0
dec1=dec(i)*PI/180.D0
dec2=twomass_dec(ii)*PI/180.D0

angular_sep=COS(dec1)*COS(dec2)*COS((ra2-ra1))+SIN(dec1)*SIN(dec2)

IF (angular_sep>1.D0) THEN
angular_sep=1.D0
END IF
IF (angular_sep<-1.D0) THEN
angular_sep=-1.D0
END IF
angular_sep=ACOS(angular_sep)*180.D0/PI

IF (angular_sep<0.D0) THEN
angular_sep=-angular_sep
END IF

IF (angular_sep<ang_min) THEN

IF (angular_sep<asold) THEN
best_2mrs=ii
asold=angular_sep
END IF

END IF


END IF
END DO

IF (best_2mrs==0) THEN
active(i)=.FALSE.
ELSE
! WRITE(*,*) asold
internalid_2mrs(i)=numbergal_2mrs(best_2mrs)
id_name(i)=mrs_id(best_2mrs)
redshift(i)=c_redshift_2mrs(best_2mrs)
END IF


END IF
END DO






OPEN(50,file='igor/2MRS_simbad.txt')

! WRITE(50,*) "internal_2MRS_ID 2MRS_ID RA DEC Z_group K_t H_t J_t E_BV"
DO i=1,n 
IF (active(i)) THEN
WRITE(50,*) internalid_2mrs(i),',',TRIM(id_name(i)),',',ra(i),',',dec(i),',',redshift(i),',',MagU(i),&
',',MagB(i),',',MagV(i),',',MagRx(i),',',MagK(i),',',Magg(i),',',Magr(i)
END IF

END DO

 CLOSE(50)



DO i=1,n
good_mag=0

IF (MagU(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagB(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagV(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagRx(i)>0.D0) THEN
good_mag=good_mag+1
END IF

IF (MagK(i)>0.D0) THEN
good_mag=good_mag+1
END IF


IF (good_mag<3) THEN
active(i)=.FALSE.
END IF



END DO

 OPEN(50,file='igor/2MRS_simbad_noSDSS.txt')

! WRITE(50,*) "internal_2MRS_ID 2MRS_ID RA DEC Z_group K_t H_t J_t E_BV"
DO i=1,n 
IF (active(i)) THEN
WRITE(50,*) internalid_2mrs(i),',',TRIM(id_name(i)),',',ra(i),',',dec(i),',',redshift(i),',',MagU(i),&
',',MagB(i),',',MagV(i),',',MagRx(i),',',MagK(i)!,',',Magg(i),',',Magr(i)
END IF

END DO

 CLOSE(50)




WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'





END PROGRAM


