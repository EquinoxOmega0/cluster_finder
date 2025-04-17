PROGRAM recalibration
! declaration of variables
IMPLICIT NONE
! basic constants
integer, parameter :: dim_s=5
integer, parameter :: dim_m=8
integer :: io_err,n_2mrs_s,n_2mrs_m,n_sdss_s,n_sdss_m,i,ii,iii,hc1,hc2

double precision, allocatable :: mass_2mrs_s(:),lum_2mrs_s(:),lum2_2mrs_s(:),dist_2mrs_s(:)
double precision, allocatable :: mass_sdss_s(:),lum_sdss_s(:),lum2_sdss_s(:),dist_sdss_s(:)
double precision, allocatable :: mass_2mrs_m(:),lum_2mrs_m(:),lum2_2mrs_m(:),dist_2mrs_m(:)
double precision, allocatable :: sigma_2mrs_m(:),rad_2mrs_m(:),nvis_2mrs_m(:)
double precision, allocatable :: mass_sdss_m(:),lum_sdss_m(:),lum2_sdss_m(:),dist_sdss_m(:)
double precision, allocatable :: sigma_sdss_m(:),rad_sdss_m(:),nvis_sdss_m(:)
double precision, allocatable :: colour_sdss_m(:),colour_sdss_s(:),colour_2mrs_m(:),colour_2mrs_s(:)
double precision, allocatable :: lum3_2mrs_s(:),lum3_sdss_s(:),lum3_2mrs_m(:),lum3_sdss_m(:)
double precision, allocatable :: lum4_2mrs_s(:),lum4_sdss_s(:),lum4_2mrs_m(:),lum4_sdss_m(:)
double precision, allocatable :: lum5_2mrs_s(:),lum5_sdss_s(:),lum5_2mrs_m(:),lum5_sdss_m(:)
double precision, dimension(1:dim_s) :: Yvector_sdss_s,Yvector_2mrs_s
double precision, dimension(1:dim_s,1:dim_s) :: Amatrix_sdss_s,Amatrix_2mrs_s
double precision, dimension(1:dim_s,1:dim_s,1:dim_s) :: Bmatrices_sdss_s,Bmatrices_2mrs_s
double precision, dimension(1:dim_s) :: solution_sdss_s,solution_2mrs_s,errorbar_sdss_s,errorbar_2mrs_s
double precision, dimension(1:dim_m) :: Yvector_sdss_m,Yvector_2mrs_m
double precision, dimension(1:dim_m,1:dim_m) :: Amatrix_sdss_m,Amatrix_2mrs_m
double precision, dimension(1:dim_m,1:dim_m,1:dim_m) :: Bmatrices_sdss_m,Bmatrices_2mrs_m
double precision, dimension(1:dim_m) :: solution_sdss_m,solution_2mrs_m,errorbar_sdss_m,errorbar_2mrs_m
double precision, allocatable :: datalist_sdss_s(:,:),datalist_2mrs_s(:,:)
double precision, allocatable :: datalist_sdss_m(:,:),datalist_2mrs_m(:,:)
double precision :: rms_sdss_s,rms_2mrs_s,rms_sdss_m,rms_2mrs_m
double precision :: det_A
double precision, dimension(1:dim_s) :: det_B
double precision, dimension(1:dim_s,1:dim_s) :: help_matrix
double precision, dimension(1:(dim_s-1),1:(dim_s-1)) :: help_submatrix
double precision, dimension(1:dim_m) :: det_B2
double precision, dimension(1:dim_m,1:dim_m) :: help_matrix2
double precision, dimension(1:(dim_m-1),1:(dim_m-1)) :: help_submatrix2
double precision, allocatable :: deviation_sdss_m(:),deviation_2mrs_m(:)
double precision, allocatable :: deviation_sdss_s(:),deviation_2mrs_s(:)
double precision, allocatable :: otherside_sdss_m(:),otherside_2mrs_m(:)
double precision, allocatable :: otherside_sdss_s(:),otherside_2mrs_s(:)

double precision, dimension(0:500,0:500) :: mapmap
real :: help_x2,help_y2
integer :: help_x,help_y


double precision :: det2d,det3d,det4d,det5d,det6d,det7d,det8d,det9d,det10d,det11d



! ! define constants
! PI=ACOS(-1.D0) !pi



! get length of file
OPEN(50,file='mass/log_mass_dependences_multi_2MRS_combi.txt')
io_err=0
n_2mrs_m=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mrs_m=n_2mrs_m+1
END DO
 CLOSE(50)
n_2mrs_m=n_2mrs_m-1



! get length of file
OPEN(50,file='mass/log_mass_dependences_single_2MRS_combi.txt')
io_err=0
n_2mrs_s=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mrs_s=n_2mrs_s+1
END DO
 CLOSE(50)
n_2mrs_s=n_2mrs_s-1



! get length of file
OPEN(50,file='mass/log_mass_dependences_multi_SDSS_combi.txt')
io_err=0
n_sdss_m=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss_m=n_sdss_m+1
END DO
 CLOSE(50)
n_sdss_m=n_sdss_m-1



! get length of file
OPEN(50,file='mass/log_mass_dependences_single_SDSS_combi.txt')
io_err=0
n_sdss_s=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss_s=n_sdss_s+1
END DO
 CLOSE(50)
n_sdss_s=n_sdss_s-1


! allocate arrays


allocate(mass_2mrs_s(1:n_2mrs_s))
allocate(lum_2mrs_s(1:n_2mrs_s))
allocate(lum2_2mrs_s(1:n_2mrs_s))
allocate(dist_2mrs_s(1:n_2mrs_s))
allocate(colour_2mrs_s(1:n_2mrs_s))
allocate(lum3_2mrs_s(1:n_2mrs_s))
allocate(lum4_2mrs_s(1:n_2mrs_s))
allocate(lum5_2mrs_s(1:n_2mrs_s))

allocate(mass_sdss_s(1:n_sdss_s))
allocate(lum_sdss_s(1:n_sdss_s))
allocate(lum2_sdss_s(1:n_sdss_s))
allocate(dist_sdss_s(1:n_sdss_s))
allocate(colour_sdss_s(1:n_sdss_s))
allocate(lum3_sdss_s(1:n_sdss_s))
allocate(lum4_sdss_s(1:n_sdss_s))
allocate(lum5_sdss_s(1:n_sdss_s))

allocate(mass_2mrs_m(1:n_2mrs_m))
allocate(lum_2mrs_m(1:n_2mrs_m))
allocate(lum2_2mrs_m(1:n_2mrs_m))
allocate(dist_2mrs_m(1:n_2mrs_m))
allocate(sigma_2mrs_m(1:n_2mrs_m))
allocate(rad_2mrs_m(1:n_2mrs_m))
allocate(nvis_2mrs_m(1:n_2mrs_m))
allocate(colour_2mrs_m(1:n_2mrs_m))
allocate(lum3_2mrs_m(1:n_2mrs_m))
allocate(lum4_2mrs_m(1:n_2mrs_m))
allocate(lum5_2mrs_m(1:n_2mrs_m))

allocate(mass_sdss_m(1:n_sdss_m))
allocate(lum_sdss_m(1:n_sdss_m))
allocate(lum2_sdss_m(1:n_sdss_m))
allocate(dist_sdss_m(1:n_sdss_m))
allocate(sigma_sdss_m(1:n_sdss_m))
allocate(rad_sdss_m(1:n_sdss_m))
allocate(nvis_sdss_m(1:n_sdss_m))
allocate(colour_sdss_m(1:n_sdss_m))
allocate(lum3_sdss_m(1:n_sdss_m))
allocate(lum4_sdss_m(1:n_sdss_m))
allocate(lum5_sdss_m(1:n_sdss_m))

 
!  dim_s=6
!  dim_m=9
!  


allocate(datalist_sdss_s(1:(dim_s+1),1:n_sdss_s))
allocate(datalist_2mrs_s(1:(dim_s+1),1:n_2mrs_s))
allocate(datalist_sdss_m(1:(dim_m+1),1:n_sdss_m))
allocate(datalist_2mrs_m(1:(dim_m+1),1:n_2mrs_m))

allocate(deviation_sdss_m(1:n_sdss_m))
allocate(deviation_2mrs_m(1:n_2mrs_m))
allocate(deviation_sdss_s(1:n_sdss_s))
allocate(deviation_2mrs_s(1:n_2mrs_s))
allocate(otherside_sdss_m(1:n_sdss_m))
allocate(otherside_2mrs_m(1:n_2mrs_m))
allocate(otherside_sdss_s(1:n_sdss_s))
allocate(otherside_2mrs_s(1:n_2mrs_s))

! read file
OPEN(50,file='mass/log_mass_dependences_single_2MRS_combi.txt')
DO i=1,n_2mrs_s
READ(50,*) mass_2mrs_s(i),lum_2mrs_s(i),dist_2mrs_s(i),colour_2mrs_s(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_2mrs_s,' galaxies'


! read file
OPEN(50,file='mass/log_mass_dependences_single_SDSS_combi.txt')
DO i=1,n_sdss_s
READ(50,*) mass_sdss_s(i),lum_sdss_s(i),dist_sdss_s(i),colour_sdss_s(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_sdss_s,' galaxies'


! read file
OPEN(50,file='mass/log_mass_dependences_multi_2MRS_combi.txt')
DO i=1,n_2mrs_m
READ(50,*) mass_2mrs_m(i),lum_2mrs_m(i),dist_2mrs_m(i),colour_2mrs_m(i),&
sigma_2mrs_m(i),rad_2mrs_m(i),nvis_2mrs_m(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_2mrs_m,' galaxies'


! read file
OPEN(50,file='mass/log_mass_dependences_multi_SDSS_combi.txt')
DO i=1,n_sdss_m
READ(50,*) mass_sdss_m(i),lum_sdss_m(i),dist_sdss_m(i),colour_sdss_m(i),&
sigma_sdss_m(i),rad_sdss_m(i),nvis_sdss_m(i)
END DO
 CLOSE(50)
WRITE(*,*) 'read in data for ',n_sdss_m,' galaxies'


 

DO i=1,n_2mrs_s
lum2_2mrs_s(i)=lum_2mrs_s(i)**2
END DO

DO i=1,n_sdss_s
lum2_sdss_s(i)=lum_sdss_s(i)**2
END DO

DO i=1,n_2mrs_m
lum2_2mrs_m(i)=lum_2mrs_m(i)**2
END DO

DO i=1,n_sdss_m
lum2_sdss_m(i)=lum_sdss_m(i)**2
END DO
 
 
 DO i=1,n_2mrs_s
lum3_2mrs_s(i)=lum_2mrs_s(i)**3
END DO

DO i=1,n_sdss_s
lum3_sdss_s(i)=lum_sdss_s(i)**3
END DO

DO i=1,n_2mrs_m
lum3_2mrs_m(i)=lum_2mrs_m(i)**3
END DO

DO i=1,n_sdss_m
lum3_sdss_m(i)=lum_sdss_m(i)**3
END DO
 
DO i=1,n_2mrs_s
lum4_2mrs_s(i)=lum_2mrs_s(i)**4
END DO

DO i=1,n_sdss_s
lum4_sdss_s(i)=lum_sdss_s(i)**4
END DO

DO i=1,n_2mrs_m
lum4_2mrs_m(i)=lum_2mrs_m(i)**4
END DO

DO i=1,n_sdss_m
lum4_sdss_m(i)=lum_sdss_m(i)**4
END DO
! 
DO i=1,n_2mrs_s
lum5_2mrs_s(i)=lum_2mrs_s(i)**5
END DO

DO i=1,n_sdss_s
lum5_sdss_s(i)=lum_sdss_s(i)**5
END DO

DO i=1,n_2mrs_m
lum5_2mrs_m(i)=lum_2mrs_m(i)**5
END DO

DO i=1,n_sdss_m
lum5_sdss_m(i)=lum_sdss_m(i)**5
END DO


! 
! 
!  DO i=1,n_sdss_s
!  datalist_sdss_s(1,i)=mass_sdss_s(i)
!  datalist_sdss_s(2,i)=lum_sdss_s(i)
!  datalist_sdss_s(3,i)=lum2_sdss_s(i)  
!  datalist_sdss_s(4,i)=lum3_sdss_s(i)  
!  datalist_sdss_s(5,i)=lum4_sdss_s(i) 
!  datalist_sdss_s(6,i)=lum5_sdss_s(i)
!  datalist_sdss_s(7,i)=dist_sdss_s(i)
!  datalist_sdss_s(8,i)=colour_sdss_s(i) 
!  datalist_sdss_s(9,i)=1.D0
!  END DO
!  
!  DO i=1,n_2mrs_s
!  datalist_2mrs_s(1,i)=mass_2mrs_s(i)
!  datalist_2mrs_s(2,i)=lum_2mrs_s(i) 
!  datalist_2mrs_s(3,i)=lum2_2mrs_s(i) 
!  datalist_2mrs_s(4,i)=lum3_2mrs_s(i)
!  datalist_2mrs_s(5,i)=lum4_2mrs_s(i)
!  datalist_2mrs_s(6,i)=lum5_2mrs_s(i) 
!  datalist_2mrs_s(7,i)=dist_2mrs_s(i)
!  datalist_2mrs_s(8,i)=colour_2mrs_s(i) 
!  datalist_2mrs_s(9,i)=1.D0
!  END DO
!  
! 
!  
!  DO i=1,n_sdss_m
!  datalist_sdss_m(1,i)=mass_sdss_m(i)
!  datalist_sdss_m(2,i)=lum_sdss_m(i) 
!  datalist_sdss_m(3,i)=lum2_sdss_m(i)
!  datalist_sdss_m(4,i)=lum3_sdss_m(i)
!  datalist_sdss_m(5,i)=lum4_sdss_m(i)
!  datalist_sdss_m(6,i)=lum5_sdss_m(i)
!  datalist_sdss_m(7,i)=dist_sdss_m(i)
!  datalist_sdss_m(8,i)=sigma_sdss_m(i) 
!  datalist_sdss_m(9,i)=rad_sdss_m(i) 
!  datalist_sdss_m(10,i)=nvis_sdss_m(i)
!  datalist_sdss_m(11,i)=colour_sdss_m(i) 
!  datalist_sdss_m(12,i)=1.D0
!  END DO
!  
!  DO i=1,n_2mrs_m
!  datalist_2mrs_m(1,i)=mass_2mrs_m(i)
!  datalist_2mrs_m(2,i)=lum_2mrs_m(i)
!  datalist_2mrs_m(3,i)=lum2_2mrs_m(i)
!  datalist_2mrs_m(4,i)=lum3_2mrs_m(i)
!  datalist_2mrs_m(5,i)=lum4_2mrs_m(i)
!  datalist_2mrs_m(6,i)=lum5_2mrs_m(i)
!  datalist_2mrs_m(7,i)=dist_2mrs_m(i)
!  datalist_2mrs_m(8,i)=sigma_2mrs_m(i)
!  datalist_2mrs_m(9,i)=rad_2mrs_m(i)
!  datalist_2mrs_m(10,i)=nvis_2mrs_m(i)
!  datalist_2mrs_m(11,i)=colour_2mrs_m(i)  
!  datalist_2mrs_m(12,i)=1.D0
!  END DO
!  
!  
 
 
 

 DO i=1,n_sdss_s
 datalist_sdss_s(1,i)=mass_sdss_s(i)
 datalist_sdss_s(2,i)=lum_sdss_s(i)
 datalist_sdss_s(3,i)=lum2_sdss_s(i)  
 datalist_sdss_s(4,i)=lum3_sdss_s(i)  
 datalist_sdss_s(5,i)=dist_sdss_s(i)
!  datalist_sdss_s(6,i)=colour_sdss_s(i) 
 datalist_sdss_s(6,i)=1.D0
 END DO
 
 DO i=1,n_2mrs_s
 datalist_2mrs_s(1,i)=mass_2mrs_s(i)
 datalist_2mrs_s(2,i)=lum_2mrs_s(i) 
 datalist_2mrs_s(3,i)=lum2_2mrs_s(i) 
 datalist_2mrs_s(4,i)=lum3_2mrs_s(i)
 datalist_2mrs_s(5,i)=dist_2mrs_s(i)
!  datalist_2mrs_s(6,i)=colour_2mrs_s(i) 
 datalist_2mrs_s(6,i)=1.D0
 END DO
 

  
 
 DO i=1,n_sdss_m
 datalist_sdss_m(1,i)=mass_sdss_m(i)
 datalist_sdss_m(2,i)=lum_sdss_m(i) 
 datalist_sdss_m(3,i)=lum2_sdss_m(i)
 datalist_sdss_m(4,i)=lum3_sdss_m(i)
 datalist_sdss_m(5,i)=dist_sdss_m(i)
!  datalist_sdss_m(6,i)=colour_sdss_m(i) 
 datalist_sdss_m(6,i)=sigma_sdss_m(i) 
 datalist_sdss_m(7,i)=rad_sdss_m(i) 
 datalist_sdss_m(8,i)=nvis_sdss_m(i)
 datalist_sdss_m(9,i)=1.D0
 END DO
 
 DO i=1,n_2mrs_m
 datalist_2mrs_m(1,i)=mass_2mrs_m(i)
 datalist_2mrs_m(2,i)=lum_2mrs_m(i)
 datalist_2mrs_m(3,i)=lum2_2mrs_m(i)
 datalist_2mrs_m(4,i)=lum3_2mrs_m(i)
 datalist_2mrs_m(5,i)=dist_2mrs_m(i)
!  datalist_2mrs_m(6,i)=colour_2mrs_m(i)  
 datalist_2mrs_m(6,i)=sigma_2mrs_m(i)
 datalist_2mrs_m(7,i)=rad_2mrs_m(i)
 datalist_2mrs_m(8,i)=nvis_2mrs_m(i)
 datalist_2mrs_m(9,i)=1.D0
 END DO
 
 
 
 
 
 DO i=1,dim_s
 Yvector_sdss_s(i)=0.D0
 DO ii=1,dim_s
 Amatrix_sdss_s(i,ii)=0.D0
 DO iii=1,dim_s
 Bmatrices_sdss_s(i,ii,iii)=0.D0 
 END DO
 END DO
 END DO
 


 DO i=1,n_sdss_s
 
 DO ii=1,dim_s
 DO iii=1,dim_s
 Amatrix_sdss_s(ii,iii)=Amatrix_sdss_s(ii,iii)+datalist_sdss_s(ii+1,i)*datalist_sdss_s(iii+1,i)
 END DO
 END DO
 
 DO ii=1,dim_s
 Yvector_sdss_s(ii)=Yvector_sdss_s(ii)+datalist_sdss_s(1,i)*datalist_sdss_s(ii+1,i)
 END DO
 
 END DO
 
 
det_A=det5d(Amatrix_sdss_s)       !                    determinate

DO iii=1,dim_s
DO i=1,dim_s
DO ii=1,dim_s
Bmatrices_sdss_s(i,ii,iii)=Amatrix_sdss_s(i,ii)
END DO
END DO
END DO

DO iii=1,dim_s
DO i=1,dim_s
Bmatrices_sdss_s(iii,i,iii)=Yvector_sdss_s(i)
END DO
END DO

DO iii=1,dim_s
DO i=1,dim_s
DO ii=1,dim_s
help_matrix(i,ii)=Bmatrices_sdss_s(i,ii,iii)
END DO
END DO
det_B(iii)=det5d(help_matrix)                           ! determinate
END DO


DO iii=1,dim_s
solution_sdss_s(iii)=det_B(iii)/det_A
END DO
! get solutions



!calculate right side of fundamental plane
DO i=1,n_sdss_s
otherside_sdss_s(i)=0.D0
DO ii=1,dim_s
otherside_sdss_s(i)=otherside_sdss_s(i)+solution_sdss_s(ii)*datalist_sdss_s(ii+1,i)
END DO 
END DO

! propagation of error
!root mean square
rms_sdss_s=0.D0


DO i=1,n_sdss_s
deviation_sdss_s(i)=otherside_sdss_s(i)-datalist_sdss_s(1,i)
rms_sdss_s=rms_sdss_s+(deviation_sdss_s(i)**2)

END DO
rms_sdss_s=SQRT(rms_sdss_s/DBLE(n_sdss_s))



! error of fitting parameters


DO i=1,dim_s

hc1=0
DO ii=1,dim_s
IF (ii.NE.i) THEN
hc1=hc1+1
hc2=0
DO iii=1,dim_s
IF (iii.NE.i) THEN
hc2=hc2+1
help_submatrix(hc1,hc2)=Amatrix_sdss_s(ii,iii)
END IF
END DO
END IF
END DO


errorbar_sdss_s(i)=det4d(help_submatrix)/det_A            !subdeterminate
! IF (errorbar_sdss_s(i)<0.D0) THEN
! WRITE(*,*) errorbar_sdss_s(i)
! END IF
errorbar_sdss_s(i)=SQRT(errorbar_sdss_s(i))
END DO

! output fitparameter results
! OPEN(60,file='mass/plane_sdss_s.txt')
! DO i=1,n_sdss_s
! WRITE(60,*) datalist_sdss_s(1,i),otherside_sdss_s(i)
! END DO
! CLOSE(60)



 
! output fitparameter results
OPEN(60,file='mass/fitparameters_SDSS_s.txt')
WRITE(60,*) solution_sdss_s(1:dim_s)
WRITE(60,*) errorbar_sdss_s(1:dim_s)
WRITE(60,*) rms_sdss_s
 CLOSE(60)


 
WRITE(*,*) 'SDSS single fitted'












 
 DO i=1,dim_s
 Yvector_2mrs_s(i)=0.D0
 DO ii=1,dim_s
 Amatrix_2mrs_s(i,ii)=0.D0
 DO iii=1,dim_s
 Bmatrices_2mrs_s(i,ii,iii)=0.D0 
 END DO
 END DO
 END DO
 


 DO i=1,n_2mrs_s
 
 DO ii=1,dim_s
 DO iii=1,dim_s
 Amatrix_2mrs_s(ii,iii)=Amatrix_2mrs_s(ii,iii)+datalist_2mrs_s(ii+1,i)*datalist_2mrs_s(iii+1,i)
 END DO
 END DO
 
 DO ii=1,dim_s
 Yvector_2mrs_s(ii)=Yvector_2mrs_s(ii)+datalist_2mrs_s(1,i)*datalist_2mrs_s(ii+1,i)
 END DO
 
 END DO
 
 
 
det_A=det5d(Amatrix_2mrs_s)                          !determinate

DO iii=1,dim_s
DO i=1,dim_s
DO ii=1,dim_s
Bmatrices_2mrs_s(i,ii,iii)=Amatrix_2mrs_s(i,ii)
END DO
END DO
END DO

DO iii=1,dim_s
DO i=1,dim_s
Bmatrices_2mrs_s(iii,i,iii)=Yvector_2mrs_s(i)
END DO
END DO

DO iii=1,dim_s
DO i=1,dim_s
DO ii=1,dim_s
help_matrix(i,ii)=Bmatrices_2mrs_s(i,ii,iii)
END DO
END DO
det_B(iii)=det5d(help_matrix)                          !   determinate
END DO


DO iii=1,dim_s
solution_2mrs_s(iii)=det_B(iii)/det_A
END DO
! get solutions



!calculate right side of fundamental plane
DO i=1,n_2mrs_s
otherside_2mrs_s(i)=0.D0
DO ii=1,dim_s
otherside_2mrs_s(i)=otherside_2mrs_s(i)+solution_2mrs_s(ii)*datalist_2mrs_s(ii+1,i)
END DO 
END DO

! propagation of error
!root mean square
rms_2mrs_s=0.D0


DO i=1,n_2mrs_s
deviation_2mrs_s(i)=otherside_2mrs_s(i)-datalist_2mrs_s(1,i)
rms_2mrs_s=rms_2mrs_s+(deviation_2mrs_s(i)**2)

END DO
rms_2mrs_s=SQRT(rms_2mrs_s/DBLE(n_2mrs_s))



! error of fitting parameters


DO i=1,dim_s

hc1=0
DO ii=1,dim_s
IF (ii.NE.i) THEN
hc1=hc1+1
hc2=0
DO iii=1,dim_s
IF (iii.NE.i) THEN
hc2=hc2+1
help_submatrix(hc1,hc2)=Amatrix_2mrs_s(ii,iii)
END IF
END DO
END IF
END DO


errorbar_2mrs_s(i)=det4d(help_submatrix)/det_A                   !subdeterminate
errorbar_2mrs_s(i)=SQRT(errorbar_2mrs_s(i))
END DO

! output fitparameter results
! OPEN(60,file='mass/plane_2mrs_s.txt')
! DO i=1,n_2mrs_s
! WRITE(60,*) datalist_2mrs_s(1,i),otherside_2mrs_s(i)
! END DO
! CLOSE(60)



 
! output fitparameter results
OPEN(60,file='mass/fitparameters_2MRS_s.txt')
WRITE(60,*) solution_2mrs_s(1:dim_s)
WRITE(60,*) errorbar_2mrs_s(1:dim_s)
WRITE(60,*) rms_2mrs_s
 CLOSE(60)


 
WRITE(*,*) '2MRS single fitted'









 
 
 
 DO i=1,dim_m
 Yvector_sdss_m(i)=0.D0
 DO ii=1,dim_m
 Amatrix_sdss_m(i,ii)=0.D0
 DO iii=1,dim_m
 Bmatrices_sdss_m(i,ii,iii)=0.D0 
 END DO
 END DO
 END DO
 


 DO i=1,n_sdss_m
 
 DO ii=1,dim_m
 DO iii=1,dim_m
 Amatrix_sdss_m(ii,iii)=Amatrix_sdss_m(ii,iii)+datalist_sdss_m(ii+1,i)*datalist_sdss_m(iii+1,i)
 END DO
 END DO
 
 DO ii=1,dim_m
 Yvector_sdss_m(ii)=Yvector_sdss_m(ii)+datalist_sdss_m(1,i)*datalist_sdss_m(ii+1,i)
 END DO
 
 END DO
 
 
 
det_A=det8d(Amatrix_sdss_m)                           !determinate

DO iii=1,dim_m
DO i=1,dim_m
DO ii=1,dim_m
Bmatrices_sdss_m(i,ii,iii)=Amatrix_sdss_m(i,ii)
END DO
END DO
END DO

DO iii=1,dim_m
DO i=1,dim_m
Bmatrices_sdss_m(iii,i,iii)=Yvector_sdss_m(i)
END DO
END DO

DO iii=1,dim_m
DO i=1,dim_m
DO ii=1,dim_m
help_matrix2(i,ii)=Bmatrices_sdss_m(i,ii,iii)
END DO
END DO
det_B2(iii)=det8d(help_matrix2)                                       !determinate
END DO


DO iii=1,dim_m
solution_sdss_m(iii)=det_B2(iii)/det_A
END DO
! get solutions



!calculate right side of fundamental plane
DO i=1,n_sdss_m
otherside_sdss_m(i)=0.D0
DO ii=1,dim_m
otherside_sdss_m(i)=otherside_sdss_m(i)+solution_sdss_m(ii)*datalist_sdss_m(ii+1,i)
END DO 
END DO

! propagation of error
!root mean square
rms_sdss_m=0.D0


DO i=1,n_sdss_m
deviation_sdss_m(i)=otherside_sdss_m(i)-datalist_sdss_m(1,i)
rms_sdss_m=rms_sdss_m+(deviation_sdss_m(i)**2)

END DO
rms_sdss_m=SQRT(rms_sdss_m/DBLE(n_sdss_m))



! error of fitting parameters


DO i=1,dim_m

hc1=0
DO ii=1,dim_m
IF (ii.NE.i) THEN
hc1=hc1+1
hc2=0
DO iii=1,dim_m
IF (iii.NE.i) THEN
hc2=hc2+1
help_submatrix2(hc1,hc2)=Amatrix_sdss_m(ii,iii)
END IF
END DO
END IF
END DO


errorbar_sdss_m(i)=det7d(help_submatrix2)/det_A                              !subdeterminate
errorbar_sdss_m(i)=SQRT(errorbar_sdss_m(i))
END DO

! ! output fitparameter results
! ! OPEN(60,file='mass/plane_sdss_m.txt')
! ! DO i=1,n_sdss_m
! ! WRITE(60,*) datalist_sdss_m(1,i),otherside_sdss_m(i)
! ! END DO
! ! CLOSE(60)
! 


 
! output fitparameter results
OPEN(60,file='mass/fitparameters_SDSS_m.txt')
WRITE(60,*) solution_sdss_m(1:dim_m)
WRITE(60,*) errorbar_sdss_m(1:dim_m)
WRITE(60,*) rms_sdss_m
 CLOSE(60)


 
WRITE(*,*) 'SDSS multi fitted'





! 
! 
! 
! 
! 
! 

 
 DO i=1,dim_m
 Yvector_2mrs_m(i)=0.D0
 DO ii=1,dim_m
 Amatrix_2mrs_m(i,ii)=0.D0
 DO iii=1,dim_m
 Bmatrices_2mrs_m(i,ii,iii)=0.D0 
 END DO
 END DO
 END DO
 


 DO i=1,n_2mrs_m
 
 DO ii=1,dim_m
 DO iii=1,dim_m
 Amatrix_2mrs_m(ii,iii)=Amatrix_2mrs_m(ii,iii)+datalist_2mrs_m(ii+1,i)*datalist_2mrs_m(iii+1,i)
 END DO
 END DO
 
 DO ii=1,dim_m
 Yvector_2mrs_m(ii)=Yvector_2mrs_m(ii)+datalist_2mrs_m(1,i)*datalist_2mrs_m(ii+1,i)
 END DO
 
 END DO
 
 
 
det_A=det8d(Amatrix_2mrs_m)                              !determinate

DO iii=1,dim_m
DO i=1,dim_m
DO ii=1,dim_m
Bmatrices_2mrs_m(i,ii,iii)=Amatrix_2mrs_m(i,ii)
END DO
END DO
END DO

DO iii=1,dim_m
DO i=1,dim_m
Bmatrices_2mrs_m(iii,i,iii)=Yvector_2mrs_m(i)
END DO
END DO

DO iii=1,dim_m
DO i=1,dim_m
DO ii=1,dim_m
help_matrix2(i,ii)=Bmatrices_2mrs_m(i,ii,iii)
END DO
END DO
det_B2(iii)=det8d(help_matrix2)                            !determinate
END DO


DO iii=1,dim_m
solution_2mrs_m(iii)=det_B2(iii)/det_A
END DO
! get solutions



!calculate right side of fundamental plane
DO i=1,n_2mrs_m
otherside_2mrs_m(i)=0.D0
DO ii=1,dim_m
otherside_2mrs_m(i)=otherside_2mrs_m(i)+solution_2mrs_m(ii)*datalist_2mrs_m(ii+1,i)
END DO 
END DO

! propagation of error
!root mean square
rms_2mrs_m=0.D0


DO i=1,n_2mrs_m
deviation_2mrs_m(i)=otherside_2mrs_m(i)-datalist_2mrs_m(1,i)
rms_2mrs_m=rms_2mrs_m+(deviation_2mrs_m(i)**2)

END DO
rms_2mrs_m=SQRT(rms_2mrs_m/DBLE(n_2mrs_m))



! error of fitting parameters


DO i=1,dim_m

hc1=0
DO ii=1,dim_m
IF (ii.NE.i) THEN
hc1=hc1+1
hc2=0
DO iii=1,dim_m
IF (iii.NE.i) THEN
hc2=hc2+1
help_submatrix2(hc1,hc2)=Amatrix_2mrs_m(ii,iii)
END IF
END DO
END IF
END DO


errorbar_2mrs_m(i)=det7d(help_submatrix2)/det_A                  !subdeterminate
errorbar_2mrs_m(i)=SQRT(errorbar_2mrs_m(i))
END DO

! output fitparameter results
! OPEN(60,file='mass/plane_2mrs_m.txt')
! DO i=1,n_2mrs_m
! WRITE(60,*) datalist_2mrs_m(1,i),otherside_2mrs_m(i)
! END DO
! CLOSE(60)



 
! output fitparameter results
OPEN(60,file='mass/fitparameters_2MRS_m.txt')
WRITE(60,*) solution_2mrs_m(1:dim_m)
WRITE(60,*) errorbar_2mrs_m(1:dim_m)
WRITE(60,*) rms_2mrs_m
 CLOSE(60)


 
WRITE(*,*) '2MRS multi fitted'






















DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_m.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_m.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(1,i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+10.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
 
 
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_m.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_s.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_m.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(2,i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_s.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+7.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

 
 
 
 
 
 
 
 
 
 
 
 
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(5,i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_m.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(5,i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(5,i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_m.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(5,i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_s.txt')
DO i=0,250
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

 
 
 
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(6,i)+0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_colour_2mrs_m.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)-0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_s
help_x=NINT((datalist_2mrs_s(6,i)+0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_colour_2mrs_s.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)-0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 
DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(6,i)+0.5D0)*50.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_colour_sdss_m.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)-0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_s
help_x=NINT((datalist_sdss_s(6,i)+0.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_colour_sdss_s.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)-0.5D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)

 

 
 
 
 

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(7,i)-1.4D0)*50.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_sigma_2mrs_m.txt')
DO i=0,100
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.4D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(7,i)-1.4D0)*50.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_sigma_sdss_m.txt')
DO i=0,100
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.4D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 

 
 

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(8,i))*50.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_rad_2mrs_m.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(8,i))*50.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF
END DO

OPEN(61,file='mass/res_rad_sdss_m.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 
 
 

 
  

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_2mrs_m
help_x=NINT((datalist_2mrs_m(9,i))*50.D0)
help_y=NINT((deviation_2mrs_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<100.D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+0.05D0
END IF
END IF
END IF
END DO

OPEN(61,file='mass/res_nvis_2mrs_m.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_sdss_m
help_x=NINT((datalist_sdss_m(9,i))*50.D0)
help_y=NINT((deviation_sdss_m(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
IF (mapmap(help_x,help_y)<100.D0) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
ELSE
mapmap(help_x,help_y)=mapmap(help_x,help_y)+0.05D0
END IF
END IF
END IF
END DO

OPEN(61,file='mass/res_nvis_sdss_m.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 OPEN(61,file='mass/all_residuals_sdss_s.txt')
DO i=1,n_sdss_s
WRITE(61,*) deviation_sdss_s(i),datalist_sdss_s(2:dim_s,i)
END DO
 CLOSE(61)

 
  OPEN(61,file='mass/all_residuals_sdss_m.txt')
DO i=1,n_sdss_m
WRITE(61,*) deviation_sdss_m(i),datalist_sdss_m(2:dim_m,i)
END DO
 CLOSE(61)
  
 
 OPEN(61,file='mass/all_residuals_2mrs_s.txt')
DO i=1,n_2mrs_s
WRITE(61,*) deviation_2mrs_s(i),datalist_2mrs_s(2:dim_s,i)
END DO
 CLOSE(61)

 
  OPEN(61,file='mass/all_residuals_2mrs_m.txt')
DO i=1,n_2mrs_m
WRITE(61,*) deviation_2mrs_m(i),datalist_2mrs_m(2:dim_m,i)
END DO
 CLOSE(61)
 
 
 
 
 
 




END PROGRAM 











double precision FUNCTION det11d(matrix)

double precision, dimension(1:11,1:11) :: matrix
double precision, dimension(1:10,1:10) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det10d


helpsum=0.D0

DO ii=1,11
 c=0
DO iii=1,10

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,11
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det10d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det10d(submatrix)*matrix(ii,1)
END IF
END DO


 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det11d'
 END IF 
 

 det11d=helpsum

! 
RETURN
END FUNCTION det11d








double precision FUNCTION det10d(matrix)

double precision, dimension(1:10,1:10) :: matrix
double precision, dimension(1:9,1:9) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det9d


helpsum=0.D0

DO ii=1,10
 c=0
DO iii=1,9

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,10
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det9d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det9d(submatrix)*matrix(ii,1)
END IF
END DO

 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det10d'
 END IF 
 


 det10d=helpsum

! 
RETURN
END FUNCTION det10d












double precision FUNCTION det9d(matrix)

double precision, dimension(1:9,1:9) :: matrix
double precision, dimension(1:8,1:8) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det8d


helpsum=0.D0

DO ii=1,9
 c=0
DO iii=1,8

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,9
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det8d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det8d(submatrix)*matrix(ii,1)
END IF
END DO


 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det9d'
 END IF 
 

 det9d=helpsum

! 
RETURN
END FUNCTION det9d











double precision FUNCTION det8d(matrix)

double precision, dimension(1:8,1:8) :: matrix
double precision, dimension(1:7,1:7) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det7d


helpsum=0.D0

DO ii=1,8
 c=0
DO iii=1,7

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,8
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det7d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det7d(submatrix)*matrix(ii,1)
END IF
END DO


 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det8d'
 END IF 
 

 det8d=helpsum

! 
RETURN
END FUNCTION det8d








double precision FUNCTION det7d(matrix)

double precision, dimension(1:7,1:7) :: matrix
double precision, dimension(1:6,1:6) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det6d


helpsum=0.D0

DO ii=1,7
 c=0
DO iii=1,6

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,7
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det6d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det6d(submatrix)*matrix(ii,1)
END IF
END DO


 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det7d'
 END IF 
 

 det7d=helpsum

! 
RETURN
END FUNCTION det7d






double precision FUNCTION det6d(matrix)

double precision, dimension(1:6,1:6) :: matrix
double precision, dimension(1:5,1:5) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det5d


helpsum=0.D0

DO ii=1,6
 c=0
DO iii=1,5

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,6
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det5d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det5d(submatrix)*matrix(ii,1)
END IF
END DO


 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det6d'
 END IF 
 

 det6d=helpsum

! 
RETURN
END FUNCTION det6d





double precision FUNCTION det5d(matrix)

double precision, dimension(1:5,1:5) :: matrix
double precision, dimension(1:4,1:4) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det4d


helpsum=0.D0

DO ii=1,5
 c=0
DO iii=1,4

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,5
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det4d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det4d(submatrix)*matrix(ii,1)
END IF
END DO

 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det5d'
 END IF 
 


 det5d=helpsum

! 
RETURN
END FUNCTION det5d







double precision FUNCTION det4d(matrix)

double precision, dimension(1:4,1:4) :: matrix
double precision, dimension(1:3,1:3) :: submatrix
double precision :: helpsum
integer :: i,ii,iii,c
double precision :: det3d


helpsum=0.D0

DO ii=1,4
 c=0
DO iii=1,3

 IF (iii==ii) THEN
 c=1
 END IF
 
DO i=2,4
submatrix(iii,i-1)=matrix((iii+c),i)
END DO
END DO

IF (MOD(ii,2)==1) THEN
helpsum=helpsum+det3d(submatrix)*matrix(ii,1)
ELSE
helpsum=helpsum-det3d(submatrix)*matrix(ii,1)
END IF
END DO


 IF (IsNaN(helpsum)) THEN
 WRITE(*,*) 'error in det4d'
 END IF 
 

 det4d=helpsum

! 
RETURN
END FUNCTION det4d










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

 
 
 