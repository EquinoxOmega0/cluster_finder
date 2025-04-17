PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,i,ii,iii,iiii,active_count,hcount,n_gal_max,nall_gal_max,n_gal_max2,nall_gal_max2
integer, dimension(1:4) :: n_mock,nall_mock,n_mock2,nall_mock2
double precision :: dummy_id,divisor,helpflip,av_sigma,sigma_sigma,binhelp,zmax,zmin,disthelp1,disthelp2,cmv
double precision :: shiftlength

double precision, allocatable :: gx(:,:),gy(:,:),gz(:,:)
integer(kind=8), allocatable ::  gid(:,:),mockid(:,:)
logical, allocatable :: active(:,:)
double precision, allocatable :: gx2(:,:),gy2(:,:),gz2(:,:)
integer(kind=8), allocatable ::  gid2(:,:),mockid2(:,:)
logical, allocatable :: active2(:,:)
double precision, allocatable :: minigx(:,:),minigy(:,:),minigz(:,:)
logical, allocatable :: miniactive(:,:)
double precision, allocatable :: minigx2(:,:),minigy2(:,:),minigz2(:,:)
logical, allocatable :: miniactive2(:,:)
integer(kind=8) :: dummy_int

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp



integer, dimension(0:1500,0:1500) :: mapmap,map_a,map_b,map_ol
double precision :: help_x2,help_y2
integer :: help_x,help_y


OPEN(33,file='logfile.txt')


WRITE(*,*) '============================================================'
WRITE(*,*) '    programme POINT DISTRIBUTION CREATOR started'
WRITE(*,*) '============================================================'
WRITE(33,*) '============================================================'
WRITE(33,*) '    programme POINT DISTRIBUTION CREATOR started'
WRITE(33,*) '============================================================'

! random seed
 CALL SYSTEM_CLOCK(hcount)
 hcount=hcount-INT(hcount/100000)*100000
 CALL srand(hcount) 


! define constants
PI=ACOS(-1.D0)

Omega_m=0.25D0
Omega_l=0.75D0

q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
H0=73.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)

shiftlength=380.D0



! get length of file
OPEN(50,file='SDSS_mock1/SDSS_galaxies_mock1_final_simid.txt')
io_err=0
n_mock(1)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock(1)=n_mock(1)+1
END DO
 CLOSE(50)
n_mock(1)=n_mock(1)-1 ! remove header

OPEN(50,file='SDSS_mock1/SDSS_allgalaxies_mock1_id_file.txt')
io_err=0
nall_mock(1)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock(1)=nall_mock(1)+1
END DO
 CLOSE(50)
nall_mock(1)=nall_mock(1)-1 ! remove header



OPEN(50,file='SDSS_mock2/SDSS_galaxies_mock2_final_simid.txt')
io_err=0
n_mock(2)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock(2)=n_mock(2)+1
END DO
 CLOSE(50)
n_mock(2)=n_mock(2)-1 ! remove header

OPEN(50,file='SDSS_mock2/SDSS_allgalaxies_mock2_id_file.txt')
io_err=0
nall_mock(2)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock(2)=nall_mock(2)+1
END DO
 CLOSE(50)
nall_mock(2)=nall_mock(2)-1 ! remove header

OPEN(50,file='SDSS_mock5/SDSS_galaxies_mock5_final_simid.txt')
io_err=0
n_mock(3)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock(3)=n_mock(3)+1
END DO
 CLOSE(50)
n_mock(3)=n_mock(3)-1 ! remove header

OPEN(50,file='SDSS_mock5/SDSS_allgalaxies_mock5_id_file.txt')
io_err=0
nall_mock(3)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock(3)=nall_mock(3)+1
END DO
 CLOSE(50)
nall_mock(3)=nall_mock(3)-1 ! remove header


OPEN(50,file='SDSS_mock8/SDSS_galaxies_mock8_final_simid.txt')
io_err=0
n_mock(4)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock(4)=n_mock(4)+1
END DO
 CLOSE(50)
n_mock(4)=n_mock(4)-1 ! remove header

OPEN(50,file='SDSS_mock8/SDSS_allgalaxies_mock8_id_file.txt')
io_err=0
nall_mock(4)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock(4)=nall_mock(4)+1
END DO
 CLOSE(50)
nall_mock(4)=nall_mock(4)-1 ! remove header


n_gal_max=0
nall_gal_max=0


DO i=1,4
IF (n_mock(i)>n_gal_max) THEN
n_gal_max=n_mock(i)
END IF
END DO


DO i=1,4
IF (nall_mock(i)>nall_gal_max) THEN
nall_gal_max=nall_mock(i)
END IF
END DO





! get length of file
OPEN(50,file='2MRS_mock1/2MRS_galaxies_mock1_final_simid.txt')
io_err=0
n_mock2(1)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock2(1)=n_mock2(1)+1
END DO
 CLOSE(50)
n_mock2(1)=n_mock2(1)-1 ! remove header

OPEN(50,file='2MRS_mock1/2MRS_allgalaxies_mock1_id_file.txt')
io_err=0
nall_mock2(1)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock2(1)=nall_mock2(1)+1
END DO
 CLOSE(50)
nall_mock2(1)=nall_mock2(1)-1 ! remove header


OPEN(50,file='2MRS_mock2/2MRS_galaxies_mock2_final_simid.txt')
io_err=0
n_mock2(2)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock2(2)=n_mock2(2)+1
END DO
 CLOSE(50)
n_mock2(2)=n_mock2(2)-1 ! remove header

OPEN(50,file='2MRS_mock2/2MRS_allgalaxies_mock2_id_file.txt')
io_err=0
nall_mock2(2)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock2(2)=nall_mock2(2)+1
END DO
 CLOSE(50)
nall_mock2(2)=nall_mock2(2)-1 ! remove header

OPEN(50,file='2MRS_mock5/2MRS_galaxies_mock5_final_simid.txt')
io_err=0
n_mock2(3)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock2(3)=n_mock2(3)+1
END DO
 CLOSE(50)
n_mock2(3)=n_mock2(3)-1 ! remove header

OPEN(50,file='2MRS_mock5/2MRS_allgalaxies_mock5_id_file.txt')
io_err=0
nall_mock2(3)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock2(3)=nall_mock2(3)+1
END DO
 CLOSE(50)
nall_mock2(3)=nall_mock2(3)-1 ! remove header


OPEN(50,file='2MRS_mock8/2MRS_galaxies_mock8_final_simid.txt')
io_err=0
n_mock2(4)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mock2(4)=n_mock2(4)+1
END DO
 CLOSE(50)
n_mock2(4)=n_mock2(4)-1 ! remove header

OPEN(50,file='2MRS_mock8/2MRS_allgalaxies_mock8_id_file.txt')
io_err=0
nall_mock2(4)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
nall_mock2(4)=nall_mock2(4)+1
END DO
 CLOSE(50)
nall_mock2(4)=nall_mock2(4)-1 ! remove header


n_gal_max2=0
nall_gal_max2=0


DO i=1,4
IF (n_mock2(i)>n_gal_max2) THEN
n_gal_max2=n_mock2(i)
END IF
END DO


DO i=1,4
IF (nall_mock2(i)>nall_gal_max2) THEN
nall_gal_max2=nall_mock2(i)
END IF
END DO



allocate(gx(1:nall_gal_max,1:4))
allocate(gy(1:nall_gal_max,1:4))
allocate(gz(1:nall_gal_max,1:4))
allocate(gid(1:nall_gal_max,1:4))
! allocate(gido(1:nall_gal_max,1:4))

allocate(mockid(1:n_gal_max,1:4))
! allocate(mockido(1:n_gal_max,1:4))

allocate(active(1:nall_gal_max,1:9))


! read file
OPEN(50,file='SDSS_mock1/SDSS_allgalaxies_mock1_pos3D.txt')
DO i=1,nall_mock(1)

READ(50,*)  gx(i,1),gy(i,1),gz(i,1),dummy_id,dummy_id

END DO

 CLOSE(50)
 
 WRITE(*,*) 'read',nall_mock(1),'galaxies from the mock 1'

OPEN(50,file='SDSS_mock2/SDSS_allgalaxies_mock2_pos3D.txt')
DO i=1,nall_mock(2)

READ(50,*)  gx(i,2),gy(i,2),gz(i,2),dummy_id,dummy_id

END DO

 CLOSE(50)
 
  WRITE(*,*) 'read',nall_mock(1),'galaxies from the mock 2'

OPEN(50,file='SDSS_mock5/SDSS_allgalaxies_mock5_pos3D.txt')
DO i=1,nall_mock(3)

READ(50,*)  gx(i,3),gy(i,3),gz(i,3),dummy_id,dummy_id

END DO

 CLOSE(50)
 WRITE(*,*) 'read',nall_mock(1),'galaxies from the mock 5'

OPEN(50,file='SDSS_mock8/SDSS_allgalaxies_mock8_pos3D.txt')
DO i=1,nall_mock(4)

READ(50,*)  gx(i,4),gy(i,4),gz(i,4),dummy_id,dummy_id

END DO

 CLOSE(50)

 WRITE(*,*) 'read',nall_mock(1),'galaxies from the mock 8'
 
 
 

OPEN(50,file='SDSS_mock1/SDSS_allgalaxies_mock1_id_file.txt')
DO i=1,nall_mock(1)
READ(50,*)  dummy_int,gid(i,1),dummy_int
END DO
 CLOSE(50)
 
OPEN(50,file='SDSS_mock1/SDSS_galaxies_mock1_final_simid.txt')
DO i=1,n_mock(1)
READ(50,*)  dummy_int,mockid(i,1),dummy_int
END DO
 CLOSE(50)

 

 
OPEN(50,file='SDSS_mock2/SDSS_allgalaxies_mock2_id_file.txt')
DO i=1,nall_mock(2)
READ(50,*)  dummy_int,gid(i,2),dummy_int
END DO
 CLOSE(50)
OPEN(50,file='SDSS_mock2/SDSS_galaxies_mock2_final_simid.txt')
DO i=1,n_mock(2)
READ(50,*)  dummy_int,mockid(i,2),dummy_int
END DO
 CLOSE(50)

 


 

OPEN(50,file='SDSS_mock5/SDSS_allgalaxies_mock5_id_file.txt')
DO i=1,nall_mock(3)
READ(50,*)  dummy_int,gid(i,3),dummy_int
END DO
 CLOSE(50) 
OPEN(50,file='SDSS_mock5/SDSS_galaxies_mock5_final_simid.txt')
DO i=1,n_mock(3)
READ(50,*)  dummy_int,mockid(i,3),dummy_int
END DO
 CLOSE(50)

 

 
 
  
OPEN(50,file='SDSS_mock8/SDSS_allgalaxies_mock8_id_file.txt')
DO i=1,nall_mock(4)
READ(50,*)  dummy_int,gid(i,4),dummy_int
END DO
 CLOSE(50)
OPEN(50,file='SDSS_mock8/SDSS_galaxies_mock8_final_simid.txt')
DO i=1,n_mock(4)
READ(50,*)  dummy_int,mockid(i,4),dummy_int
END DO
 CLOSE(50)

 

 
 WRITE(*,*) 'all IDs read in'
 


allocate(gx2(1:nall_gal_max2,1:4))
allocate(gy2(1:nall_gal_max2,1:4))
allocate(gz2(1:nall_gal_max2,1:4))
allocate(gid2(1:nall_gal_max2,1:4))

allocate(mockid2(1:n_gal_max2,1:4))

allocate(active2(1:nall_gal_max2,1:9))


! read file
OPEN(50,file='2MRS_mock1/2MRS_allgalaxies_mock1_pos3D.txt')
DO i=1,nall_mock2(1)

READ(50,*)  gx2(i,1),gy2(i,1),gz2(i,1),dummy_id,dummy_id

END DO

 CLOSE(50)
 
 WRITE(*,*) 'read',nall_mock2(1),'galaxies from the mock 1'

OPEN(50,file='2MRS_mock2/2MRS_allgalaxies_mock2_pos3D.txt')
DO i=1,nall_mock2(2)

READ(50,*)  gx2(i,2),gy2(i,2),gz2(i,2),dummy_id,dummy_id

END DO

 CLOSE(50)
 
  WRITE(*,*) 'read',nall_mock2(1),'galaxies from the mock 2'

OPEN(50,file='2MRS_mock5/2MRS_allgalaxies_mock5_pos3D.txt')
DO i=1,nall_mock2(3)

READ(50,*)  gx2(i,3),gy2(i,3),gz2(i,3),dummy_id,dummy_id

END DO

 CLOSE(50)
 WRITE(*,*) 'read',nall_mock2(1),'galaxies from the mock 5'

OPEN(50,file='2MRS_mock8/2MRS_allgalaxies_mock8_pos3D.txt')
DO i=1,nall_mock2(4)

READ(50,*)  gx2(i,4),gy2(i,4),gz2(i,4),dummy_id,dummy_id

END DO

 CLOSE(50)

 WRITE(*,*) 'read',nall_mock2(1),'galaxies from the mock 8'
 
 
 

OPEN(50,file='2MRS_mock1/2MRS_galaxies_mock1_final_simid.txt')
DO i=1,n_mock2(1)
READ(50,*)  dummy_int,mockid2(i,1),dummy_int

END DO
 CLOSE(50)

 
OPEN(50,file='2MRS_mock1/2MRS_allgalaxies_mock1_id_file.txt')
DO i=1,nall_mock2(1)
READ(50,*)  dummy_int,gid2(i,1),dummy_int
END DO
 CLOSE(50)

 
 
 OPEN(50,file='2MRS_mock2/2MRS_galaxies_mock2_final_simid.txt')
DO i=1,n_mock2(2)
READ(50,*)  dummy_int,mockid2(i,2),dummy_int
END DO
 CLOSE(50)

 
OPEN(50,file='2MRS_mock2/2MRS_allgalaxies_mock2_id_file.txt')
DO i=1,nall_mock2(2)
READ(50,*)  dummy_int,gid2(i,2),dummy_int
END DO
 CLOSE(50)

 
 
  OPEN(50,file='2MRS_mock5/2MRS_galaxies_mock5_final_simid.txt')
DO i=1,n_mock2(3)
READ(50,*)  dummy_int,mockid2(i,3),dummy_int
END DO
 CLOSE(50)

 
OPEN(50,file='2MRS_mock5/2MRS_allgalaxies_mock5_id_file.txt')
DO i=1,nall_mock2(3)
READ(50,*)  dummy_int,gid2(i,3),dummy_int
END DO
 CLOSE(50)
 
 
 
  
  OPEN(50,file='2MRS_mock8/2MRS_galaxies_mock8_final_simid.txt')
DO i=1,n_mock2(4)
READ(50,*)  dummy_int,mockid2(i,4),dummy_int
END DO
 CLOSE(50)

 
OPEN(50,file='2MRS_mock8/2MRS_allgalaxies_mock8_id_file.txt')
DO i=1,nall_mock2(4)
READ(50,*)  dummy_int,gid2(i,4),dummy_int
END DO
 CLOSE(50)

 
 
 WRITE(*,*) 'all IDs read in'
 
 
DO i=1,nall_mock(1)
gx(i,1)=gx(i,1)
gy(i,1)=gy(i,1)
gz(i,1)=gz(i,1)
END DO
 
DO i=1,nall_mock(2)
gx(i,2)=(shiftlength/h)-gx(i,2)
gy(i,2)=gy(i,2)
gz(i,2)=gz(i,2)
END DO

DO i=1,nall_mock(3)
gx(i,3)=(shiftlength/h)-gx(i,3)
gy(i,3)=(shiftlength/h)-gy(i,3)
gz(i,3)=gz(i,3)
END DO

DO i=1,nall_mock(4)
gx(i,4)=(shiftlength/h)-gx(i,4)
gy(i,4)=(shiftlength/h)-gy(i,4)
gz(i,4)=(shiftlength/h)-gz(i,4)
END DO


DO ii=1,9
DO i=1,nall_gal_max
active(i,ii)=.FALSE.
END DO 
END DO



DO iii=1,4
i=1
! WRITE(*,*) 'ok'
DO ii=1,n_mock(iii)
! WRITE(*,*) ii,n_mock(iii)
DO WHILE ((gid(i,iii)).NE.(mockid(ii,iii)))
i=i+1
END DO
! WRITE(*,*) i,ii,n_mock(iii)
active(i,iii)=.TRUE.
i=i+1

END DO
END DO

! WRITE(*,*) 'ok'


OPEN(50,file='SDSS_overlaps/SDSS_projection1.txt')
DO i=1,nall_mock(1)
IF (active(i,1)) THEN
WRITE(50,*)  gx(i,1),gy(i,1)
END IF
END DO
 CLOSE(50)

OPEN(50,file='SDSS_overlaps/SDSS_projection2.txt')
DO i=1,nall_mock(2)
IF (active(i,2)) THEN
WRITE(50,*)  gx(i,2),gy(i,2)
END IF
END DO
 CLOSE(50)

 OPEN(50,file='SDSS_overlaps/SDSS_projection3.txt')
DO i=1,nall_mock(3)
IF (active(i,3)) THEN
WRITE(50,*)  gx(i,3),gy(i,3)
END IF
END DO
 CLOSE(50)
 
 OPEN(50,file='SDSS_overlaps/SDSS_projection4.txt')
DO i=1,nall_mock(4)
IF (active(i,4)) THEN
WRITE(50,*)  gx(i,4),gy(i,4)
END IF
END DO
 CLOSE(50)
 
 
 
 
 
 
 
!  WRITE(*,*) 'ok'
 
 
 
 
 
 
DO i=1,nall_mock2(1)
gx2(i,1)=gx2(i,1)
gy2(i,1)=gy2(i,1)
gz2(i,1)=gz2(i,1)
END DO
 
DO i=1,nall_mock2(2)
gx2(i,2)=(shiftlength/h)-gx2(i,2)
gy2(i,2)=gy2(i,2)
gz2(i,2)=gz2(i,2)
END DO

DO i=1,nall_mock2(3)
gx2(i,3)=(shiftlength/h)-gx2(i,3)
gy2(i,3)=(shiftlength/h)-gy2(i,3)
gz2(i,3)=gz2(i,3)
END DO

DO i=1,nall_mock2(4)
gx2(i,4)=(shiftlength/h)-gx2(i,4)
gy2(i,4)=(shiftlength/h)-gy2(i,4)
gz2(i,4)=(shiftlength/h)-gz2(i,4)
END DO


DO ii=1,9
DO i=1,nall_gal_max2
active2(i,ii)=.FALSE.
END DO 
END DO

! DO ii=1,4
! DO i=1,nall_mock(ii)
! 
! END DO 
! END DO

DO iii=1,4
i=1
DO ii=1,n_mock2(iii)

DO WHILE ((gid2(i,iii)).NE.(mockid2(ii,iii)))
i=i+1
END DO

active2(i,iii)=.TRUE.
i=i+1

END DO
END DO




OPEN(50,file='2MRS_overlaps/2MRS_projection1.txt')
DO i=1,nall_mock2(1)
IF (active2(i,1)) THEN
WRITE(50,*)  gx2(i,1),gy2(i,1)
END IF
END DO
 CLOSE(50)

OPEN(50,file='2MRS_overlaps/2MRS_projection2.txt')
DO i=1,nall_mock2(2)
IF (active2(i,2)) THEN
WRITE(50,*)  gx2(i,2),gy2(i,2)
END IF
END DO
 CLOSE(50)

 OPEN(50,file='2MRS_overlaps/2MRS_projection3.txt')
DO i=1,nall_mock2(3)
IF (active2(i,3)) THEN
WRITE(50,*)  gx2(i,3),gy2(i,3)
END IF
END DO
 CLOSE(50)
 
 OPEN(50,file='2MRS_overlaps/2MRS_projection4.txt')
DO i=1,nall_mock2(4)
IF (active2(i,4)) THEN
WRITE(50,*)  gx2(i,4),gy2(i,4)
END IF
END DO
 CLOSE(50)
 
 WRITE(*,*) 'all projection output written'
 
allocate(minigx(1:n_gal_max,1:4))
allocate(minigy(1:n_gal_max,1:4))
allocate(minigz(1:n_gal_max,1:4))
allocate(miniactive(1:n_gal_max,1:9))
allocate(minigx2(1:n_gal_max2,1:4))
allocate(minigy2(1:n_gal_max2,1:4))
allocate(minigz2(1:n_gal_max2,1:4))
allocate(miniactive2(1:n_gal_max2,1:9))
 
DO iii=1,4
ii=0
DO i=1,nall_mock(iii)
IF (active(i,iii)) THEN
ii=ii+1
IF (ii>n_gal_max) THEN
WRITE(*,*) 'WARNING SDSS galaxy mismatch',iii
END IF
minigx(ii,iii)=gx(i,iii)
minigy(ii,iii)=gy(i,iii)
minigz(ii,iii)=gz(i,iii)
END IF
END DO
 END DO
 
 DO iii=1,4
ii=0
DO i=1,nall_mock2(iii)
IF (active2(i,iii)) THEN
ii=ii+1
IF (ii>n_gal_max2) THEN
WRITE(*,*) 'WARNING 2MRS galaxy mismatch',iii
END IF
minigx2(ii,iii)=gx2(i,iii)
minigy2(ii,iii)=gy2(i,iii)
minigz2(ii,iii)=gz2(i,iii)
END IF
END DO
 END DO
 
 
 DO iii=1,4
 DO i=1,n_gal_max
 IF (i<=n_mock(iii)) THEN
 miniactive(i,iii)=.TRUE.
 ELSE
 miniactive(i,iii)=.FALSE.
 END IF
 END DO 
 END DO
 
 DO iii=1,4
 DO i=1,n_gal_max2
 IF (i<=n_mock2(iii)) THEN
 miniactive2(i,iii)=.TRUE.
 ELSE
 miniactive2(i,iii)=.FALSE.
 END IF
 END DO 
 END DO
 
 
 
 

 
DO iii=1,4
DO i=1,n_gal_max
miniactive(i,iii+4)=miniactive(i,iii)
miniactive(i,9)=.FALSE.
END DO
END DO
 
 
 
DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
DO ii=1,n_mock(2)
IF (miniactive(ii,6)) THEN
! WRITE(*,*) i,ii
IF (mockid(i,1)==(mockid(ii,2))) THEN

miniactive(i,5)=.FALSE.
miniactive(ii,6)=.FALSE.
miniactive(i,9)=.TRUE.

END IF

 END IF
 END DO
 END IF
 END DO
 
 
  OPEN(50,file='SDSS_overlaps/SDSS_neighbours_1a.txt')
DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
WRITE(50,*)  minigx(i,1),minigy(i,1)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='SDSS_overlaps/SDSS_neighbours_1b.txt')
DO i=1,n_mock(2)
IF (miniactive(i,6)) THEN
WRITE(50,*)  minigx(i,2),minigy(i,2)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='SDSS_overlaps/SDSS_neighbours_1overlap.txt')
DO i=1,n_mock(1)
IF (miniactive(i,9)) THEN
WRITE(50,*)  minigx(i,1),minigy(i,1)
END IF
END DO
 CLOSE(50)
 
 
 
 DO i=0,1500
DO ii=0,1500
mapmap(i,ii)=0
map_a(i,ii)=0
map_b(i,ii)=0
map_ol(i,ii)=0
END DO
END DO


DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
help_x=NINT(minigx(i,1))
help_y=NINT(minigy(i,1))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_a(help_x,help_y)=1
END IF
END IF
END IF
END DO

DO i=1,n_mock(2)
IF (miniactive(i,6)) THEN
help_x=NINT(minigx(i,2))
help_y=NINT(minigy(i,2))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_b(help_x,help_y)=1
END IF
END IF
END IF
END DO


 DO i=0,750
DO ii=0,750
IF (map_a(i,ii)>0) THEN
mapmap(i,ii)=1
END IF
IF (map_b(i,ii)>0) THEN
mapmap(i,ii)=3
END IF
IF ((map_b(i,ii)>0).AND.(map_a(i,ii)>0)) THEN
mapmap(i,ii)=2
END IF
END DO
END DO


OPEN(61,file='SDSS_overlaps/SDSS_1map.txt')
DO i=0,750
DO ii=0,750
help_x2=DBLE(i)
help_y2=DBLE(ii)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 
 
 
 WRITE(*,*) 'SDSS 1. set done'
 
  
DO iii=1,4
DO i=1,n_gal_max
miniactive(i,iii+4)=miniactive(i,iii)
miniactive(i,9)=.FALSE.
END DO
END DO
 
 
 
DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
DO ii=1,n_mock(3)
IF (miniactive(ii,7)) THEN

IF (mockid(i,1)==(mockid(ii,3))) THEN

miniactive(i,5)=.FALSE.
miniactive(ii,7)=.FALSE.
miniactive(i,9)=.TRUE.

END IF

 END IF
 END DO
 END IF
 END DO
 
 
  OPEN(50,file='SDSS_overlaps/SDSS_neighbours_2a.txt')
DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
WRITE(50,*)  minigx(i,1),minigy(i,1)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='SDSS_overlaps/SDSS_neighbours_2b.txt')
DO i=1,n_mock(3)
IF (miniactive(i,7)) THEN
WRITE(50,*)  minigx(i,3),minigy(i,3)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='SDSS_overlaps/SDSS_neighbours_2overlap.txt')
DO i=1,n_mock(1)
IF (miniactive(i,9)) THEN
WRITE(50,*)  minigx(i,1),minigy(i,1)
END IF
END DO
 CLOSE(50)
 
 
 
 
 
 
 DO i=0,1500
DO ii=0,1500
mapmap(i,ii)=0
map_a(i,ii)=0
map_b(i,ii)=0
map_ol(i,ii)=0
END DO
END DO


DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
help_x=NINT(minigx(i,1))
help_y=NINT(minigy(i,1))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_a(help_x,help_y)=1
END IF
END IF
END IF
END DO

DO i=1,n_mock(3)
IF (miniactive(i,7)) THEN
help_x=NINT(minigx(i,3))
help_y=NINT(minigy(i,3))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_b(help_x,help_y)=1
END IF
END IF
END IF
END DO


 DO i=0,750
DO ii=0,750
IF (map_a(i,ii)>0) THEN
mapmap(i,ii)=1
END IF
IF (map_b(i,ii)>0) THEN
mapmap(i,ii)=3
END IF
IF ((map_b(i,ii)>0).AND.(map_a(i,ii)>0)) THEN
mapmap(i,ii)=2
END IF
END DO
END DO


OPEN(61,file='SDSS_overlaps/SDSS_2map.txt')
DO i=0,750
DO ii=0,750
help_x2=DBLE(i)
help_y2=DBLE(ii)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 
 
 
 
 
 
 
 
 
 
 
  WRITE(*,*) 'SDSS 2. set done'
 
 
  
DO iii=1,4
DO i=1,n_gal_max
miniactive(i,iii+4)=miniactive(i,iii)
miniactive(i,9)=.FALSE.
END DO
END DO
 
 
 
DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
DO ii=1,n_mock(4)
IF (miniactive(ii,8)) THEN

IF (mockid(i,1)==(mockid(ii,4))) THEN

miniactive(i,5)=.FALSE.
miniactive(ii,8)=.FALSE.
miniactive(i,9)=.TRUE.

END IF

 END IF
 END DO
 END IF
 END DO
 
 
  OPEN(50,file='SDSS_overlaps/SDSS_neighbours_3a.txt')
DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
WRITE(50,*)  minigx(i,1),minigy(i,1)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='SDSS_overlaps/SDSS_neighbours_3b.txt')
DO i=1,n_mock(4)
IF (miniactive(i,8)) THEN
WRITE(50,*)  minigx(i,4),minigy(i,4)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='SDSS_overlaps/SDSS_neighbours_3overlap.txt')
DO i=1,n_mock(1)
IF (miniactive(i,9)) THEN
WRITE(50,*)  minigx(i,1),minigy(i,1)
END IF
END DO
 CLOSE(50)
 
 
 DO i=0,1500
DO ii=0,1500
mapmap(i,ii)=0
map_a(i,ii)=0
map_b(i,ii)=0
map_ol(i,ii)=0
END DO
END DO


DO i=1,n_mock(1)
IF (miniactive(i,5)) THEN
help_x=NINT(minigx(i,1))
help_y=NINT(minigy(i,1))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_a(help_x,help_y)=1
END IF
END IF
END IF
END DO

DO i=1,n_mock(4)
IF (miniactive(i,8)) THEN
help_x=NINT(minigx(i,4))
help_y=NINT(minigy(i,4))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_b(help_x,help_y)=1
END IF
END IF
END IF
END DO


 DO i=0,750
DO ii=0,750
IF (map_a(i,ii)>0) THEN
mapmap(i,ii)=1
END IF
IF (map_b(i,ii)>0) THEN
mapmap(i,ii)=3
END IF
IF ((map_b(i,ii)>0).AND.(map_a(i,ii)>0)) THEN
mapmap(i,ii)=2
END IF
END DO
END DO


OPEN(61,file='SDSS_overlaps/SDSS_3map.txt')
DO i=0,750
DO ii=0,750
help_x2=DBLE(i)
help_y2=DBLE(ii)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 
 
 
 
 
  WRITE(*,*) 'SDSS 3. set done'
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
DO iii=1,4
DO i=1,n_gal_max2
miniactive2(i,iii+4)=miniactive2(i,iii)
miniactive2(i,9)=.FALSE.
END DO
END DO
 
 
 
DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
DO ii=1,n_mock2(2)
IF (miniactive2(ii,6)) THEN

IF (mockid2(i,1)==(mockid2(ii,2))) THEN

miniactive2(i,5)=.FALSE.
miniactive2(ii,6)=.FALSE.
miniactive2(i,9)=.TRUE.

END IF

 END IF
 END DO
 END IF
 END DO
 
 
  OPEN(50,file='2MRS_overlaps/2MRS_neighbours_1a.txt')
DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
WRITE(50,*)  minigx2(i,1),minigy2(i,1)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='2MRS_overlaps/2MRS_neighbours_1b.txt')
DO i=1,n_mock2(2)
IF (miniactive2(i,6)) THEN
WRITE(50,*)  minigx2(i,2),minigy2(i,2)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='2MRS_overlaps/2MRS_neighbours_1overlap.txt')
DO i=1,n_mock2(1)
IF (miniactive2(i,9)) THEN
WRITE(50,*)  minigx2(i,1),minigy2(i,1)
END IF
END DO
 CLOSE(50)
 
 
 DO i=0,1500
DO ii=0,1500
mapmap(i,ii)=0
map_a(i,ii)=0
map_b(i,ii)=0
map_ol(i,ii)=0
END DO
END DO


DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
help_x=NINT(minigx2(i,1))
help_y=NINT(minigy2(i,1))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_a(help_x,help_y)=1
END IF
END IF
END IF
END DO

DO i=1,n_mock2(2)
IF (miniactive2(i,6)) THEN
help_x=NINT(minigx2(i,2))
help_y=NINT(minigy2(i,2))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_b(help_x,help_y)=1
END IF
END IF
END IF
END DO


 DO i=0,750
DO ii=0,750
IF (map_a(i,ii)>0) THEN
mapmap(i,ii)=1
END IF
IF (map_b(i,ii)>0) THEN
mapmap(i,ii)=3
END IF
IF ((map_b(i,ii)>0).AND.(map_a(i,ii)>0)) THEN
mapmap(i,ii)=2
END IF
END DO
END DO


OPEN(61,file='2MRS_overlaps/2MRS_1map.txt')
DO i=0,750
DO ii=0,750
help_x2=DBLE(i)
help_y2=DBLE(ii)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
  WRITE(*,*) '2MRS 1. set done'
 
 
  
DO iii=1,4
DO i=1,n_gal_max2
miniactive2(i,iii+4)=miniactive2(i,iii)
miniactive2(i,9)=.FALSE.
END DO
END DO
 
 
 
DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
DO ii=1,n_mock2(3)
IF (miniactive2(ii,7)) THEN

IF (mockid2(i,1)==(mockid2(ii,3))) THEN

miniactive2(i,5)=.FALSE.
miniactive2(ii,7)=.FALSE.
miniactive2(i,9)=.TRUE.

END IF

 END IF
 END DO
 END IF
 END DO
 
 
  OPEN(50,file='2MRS_overlaps/2MRS_neighbours_2a.txt')
DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
WRITE(50,*)  minigx2(i,1),minigy2(i,1)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='2MRS_overlaps/2MRS_neighbours_2b.txt')
DO i=1,n_mock2(3)
IF (miniactive2(i,7)) THEN
WRITE(50,*)  minigx2(i,3),minigy2(i,3)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='2MRS_overlaps/2MRS_neighbours_2overlap.txt')
DO i=1,n_mock2(1)
IF (miniactive2(i,9)) THEN
WRITE(50,*)  minigx2(i,1),minigy2(i,1)
END IF
END DO
 CLOSE(50)
 
 
 
 DO i=0,1500
DO ii=0,1500
mapmap(i,ii)=0
map_a(i,ii)=0
map_b(i,ii)=0
map_ol(i,ii)=0
END DO
END DO


DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
help_x=NINT(minigx2(i,1))
help_y=NINT(minigy2(i,1))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_a(help_x,help_y)=1
END IF
END IF
END IF
END DO

DO i=1,n_mock2(3)
IF (miniactive2(i,7)) THEN
help_x=NINT(minigx2(i,3))
help_y=NINT(minigy2(i,3))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_b(help_x,help_y)=1
END IF
END IF
END IF
END DO


 DO i=0,750
DO ii=0,750
IF (map_a(i,ii)>0) THEN
mapmap(i,ii)=1
END IF
IF (map_b(i,ii)>0) THEN
mapmap(i,ii)=3
END IF
IF ((map_b(i,ii)>0).AND.(map_a(i,ii)>0)) THEN
mapmap(i,ii)=2
END IF
END DO
END DO


OPEN(61,file='2MRS_overlaps/2MRS_2map.txt')
DO i=0,750
DO ii=0,750
help_x2=DBLE(i)
help_y2=DBLE(ii)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 
   WRITE(*,*) '2MRS 2. set done'
 
  
DO iii=1,4
DO i=1,n_gal_max2
miniactive2(i,iii+4)=miniactive2(i,iii)
miniactive2(i,9)=.FALSE.
END DO
END DO
 
 
 
DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
DO ii=1,n_mock2(4)
IF (miniactive2(ii,8)) THEN

IF (mockid2(i,1)==(mockid2(ii,4))) THEN

miniactive2(i,5)=.FALSE.
miniactive2(ii,8)=.FALSE.
miniactive2(i,9)=.TRUE.

END IF

 END IF
 END DO
 END IF
 END DO
 
 
  OPEN(50,file='2MRS_overlaps/2MRS_neighbours_3a.txt')
DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
WRITE(50,*)  minigx2(i,1),minigy2(i,1)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='2MRS_overlaps/2MRS_neighbours_3b.txt')
DO i=1,n_mock2(4)
IF (miniactive2(i,8)) THEN
WRITE(50,*)  minigx2(i,4),minigy2(i,4)
END IF
END DO
 CLOSE(50)
 
   OPEN(50,file='2MRS_overlaps/2MRS_neighbours_3overlap.txt')
DO i=1,n_mock2(1)
IF (miniactive2(i,9)) THEN
WRITE(50,*)  minigx2(i,1),minigy2(i,1)
END IF
END DO
 CLOSE(50)

 
 DO i=0,1500
DO ii=0,1500
mapmap(i,ii)=0
map_a(i,ii)=0
map_b(i,ii)=0
map_ol(i,ii)=0
END DO
END DO


DO i=1,n_mock2(1)
IF (miniactive2(i,5)) THEN
help_x=NINT(minigx2(i,1))
help_y=NINT(minigy2(i,1))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_a(help_x,help_y)=1
END IF
END IF
END IF
END DO

DO i=1,n_mock2(4)
IF (miniactive2(i,8)) THEN
help_x=NINT(minigx2(i,4))
help_y=NINT(minigy2(i,4))
IF ((help_x>0).AND.(help_x<750)) THEN
IF ((help_y>0).AND.(help_y<750)) THEN
map_b(help_x,help_y)=1
END IF
END IF
END IF
END DO


 DO i=0,750
DO ii=0,750
IF (map_a(i,ii)>0) THEN
mapmap(i,ii)=1
END IF
IF (map_b(i,ii)>0) THEN
mapmap(i,ii)=3
END IF
IF ((map_b(i,ii)>0).AND.(map_a(i,ii)>0)) THEN
mapmap(i,ii)=2
END IF
END DO
END DO


OPEN(61,file='2MRS_overlaps/2MRS_3map.txt')
DO i=0,750
DO ii=0,750
help_x2=DBLE(i)
help_y2=DBLE(ii)
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

  
   WRITE(*,*) '2MRS 3. set done'
 
    WRITE(*,*) '-------------------------------------------------------------'
 
 
 

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'
WRITE(33,*) '============================================================'
WRITE(33,*) '    programme complete'
WRITE(33,*) '============================================================'

CLOSE(33)


END PROGRAM


