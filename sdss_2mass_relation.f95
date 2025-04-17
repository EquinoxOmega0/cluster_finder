PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,h,Omega_m,Omega_l
integer :: io_err,n,i,ii,iii,n_sdss,n_2mass
double precision :: d_J,dummy_coeffient,d_Ks,e_J,e_Ks,f_J,f_Ks,delta_sum_J,delta_sum_Ks,dummy_id
double precision, allocatable :: gx(:),gy(:),gz(:),mag_r(:),mag_g(:),mag_i(:),mag_z(:),mag_J(:),mag_H(:),mag_Ks(:)
double precision, allocatable :: model_J(:),model_Ks(:),delta_J(:),delta_Ks(:)
double precision, allocatable :: otherside_J(:),otherside_Ks(:),deviation_J(:),deviation_Ks(:)

double precision, allocatable :: model_H(:),delta_H(:),otherside_H(:),deviation_H(:)

double precision, allocatable ::sdss_ra(:),sdss_dec(:),sdss_z(:),sdss_gl(:),sdss_gb(:),z_cor(:)
double precision, allocatable ::  cModelMag_g(:),cModelMag_r(:),cModelMag_i(:)
double precision, allocatable ::  extinction_g(:),extinction_r(:),extinction_i(:)
double precision, allocatable :: cModelMagErr_g(:),cModelMagErr_r(:),cModelMagErr_i(:)
double precision, allocatable :: twomass_ra(:),twomass_dec(:)
double precision, allocatable :: twomass_gl(:),twomass_gb(:)
double precision, allocatable :: Ktmag(:),Jtmag(:),Kmagerr(:),Jmagerr(:)
double precision, allocatable :: E_BV(:),twomass_cz(:),twomass_ecz(:) 

double precision, dimension(0:5,0:3) :: K_coeff_r,K_coeff_i,K_coeff_z,K_coeff_J,K_coeff_Ks
double precision, dimension(0:7,0:3) :: K_coeff_g,K_coeff_H
double precision, allocatable :: K_g(:),K_r(:),K_i(:)

 double precision, allocatable :: magapp_g(:),magapp_r(:),magapp_i(:)
 double precision, allocatable :: K_J(:),K_Ks(:),z_2mass(:),z_sdss(:)
 double precision, allocatable :: magapp_J(:),magapp_Ks(:)

 logical, allocatable :: twomass_active(:),sdss_active(:)
 
double precision, dimension(1:3,1:3) :: Amatrix_J,Bmatrix1_J,Amatrix_Ks,Bmatrix1_Ks
double precision, dimension(1:3,1:3) :: Bmatrix2_J,Bmatrix3_J,Bmatrix2_Ks,Bmatrix3_Ks
double precision, dimension(1:3) :: Yvector_J,Yvector_Ks
double precision :: detA_J,detA_Ks,detB1_J,detB1_Ks,detB2_J,detB2_Ks,detB3_J,detB3_Ks

double precision :: a_J,b_J,c_J,a_Ks,b_Ks,c_Ks,rms_J,rms_Ks
double precision :: a_err_J,b_err_J,c_err_J,a_err_Ks,b_err_Ks,c_err_Ks
double precision :: v_cmb,l_cmb,b_cmb,z_cmb,gal_b,gal_l,z_cmb_x,z_cmb_y,z_cmb_z,z_gal_x,z_gal_y,z_gal_z
double precision :: fivearcsec,angular_sep

double precision ::  a_J_bilir,a_H_bilir,a_Ks_bilir,b_J_bilir,b_H_bilir,b_Ks_bilir,c_J_bilir,c_H_bilir,c_Ks_bilir
double precision ::  a_J_mill,a_H_mill,a_Ks_mill,b_J_mill,b_H_mill,b_Ks_mill,c_J_mill,c_H_mill,c_Ks_mill

 double precision, allocatable :: K_H(:),magapp_H(:),Htmag(:),Hmagerr(:)
 double precision, allocatable :: leftside_J(:),leftside_H(:),leftside_Ks(:)
 double precision, allocatable :: rightside_J(:),rightside_H(:),rightside_Ks(:)
 
double precision, dimension(1:3,1:3) :: Amatrix_H,Bmatrix1_H,Bmatrix2_H,Bmatrix3_H
double precision, dimension(1:3) :: Yvector_H
double precision :: detA_H,detB1_H,detB2_H,detB3_H

double precision :: a_H,b_H,c_H,rms_H,d_H,e_H,f_H,delta_sum_H
double precision :: a_err_H,b_err_H,c_err_H,delta_z

 double precision, allocatable :: m_J(:),m_H(:),m_Ks(:),m_g(:),m_r(:),m_i(:)
 integer :: hcount

integer, dimension(0:500,0:400) :: mapmap
real :: help_x2,help_y2
integer :: help_x,help_y

 
 
 
WRITE(*,*) '============================================================'
WRITE(*,*) '    programme SDSS-2MASS RELATION started'
WRITE(*,*) '============================================================'


! define constants
PI=ACOS(-1.D0)
fivearcsec=10.D0/3600.D0*PI/180.D0

Omega_m=0.272D0
Omega_l=0.728D0

q0=Omega_m/2.D0-Omega_l
light=3.D5
H0=70.4D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2



! get length of file
OPEN(50,file='SDSSmag_2MASS_mag.txt')
io_err=0
n=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n=n+1
END DO
 CLOSE(50)
n=n-2 ! remove header




allocate(gx(1:n))
allocate(gy(1:n))
allocate(gz(1:n))
allocate(mag_r(1:n))
allocate(mag_g(1:n))
allocate(mag_i(1:n))
allocate(mag_z(1:n))
allocate(mag_J(1:n))
allocate(mag_H(1:n))
allocate(mag_Ks(1:n))


allocate(model_J(1:n))
allocate(model_Ks(1:n))
allocate(delta_J(1:n))
allocate(delta_Ks(1:n))
allocate(otherside_J(1:n))
allocate(otherside_Ks(1:n))
allocate(deviation_J(1:n))
allocate(deviation_Ks(1:n))

allocate(model_H(1:n))
allocate(delta_H(1:n))
allocate(otherside_H(1:n))
allocate(deviation_H(1:n))









! read file from milleniums simulation
OPEN(50,file='SDSSmag_2MASS_mag.txt')
READ(50,*)
DO i=1,n

READ(50,*) gx(i),gy(i),gz(i),mag_r(i),mag_g(i),mag_i(i),mag_z(i),mag_J(i),mag_H(i),mag_Ks(i)


END DO

 CLOSE(50)




WRITE(*,*) n,'galaxies'



 
 OPEN(52,file='2mass_sdss_transformation.txt')
READ(52,*) 
READ(52,*) d_J,d_H,d_Ks
READ(52,*) e_J,e_H,e_Ks
READ(52,*) f_J,f_H,f_Ks
 CLOSE(52)
 
 a_J_bilir=d_J
 a_H_bilir=d_H
 a_Ks_bilir=d_Ks
 b_J_bilir=e_J
 b_H_bilir=e_H
 b_Ks_bilir=e_Ks
 c_J_bilir=f_J
 c_H_bilir=f_H
 c_Ks_bilir=f_Ks
 
 delta_sum_J=0.D0
 delta_sum_Ks=0.D0
 delta_sum_H=0.D0
 
 DO i=1,n
 
model_J(i)=mag_g(i)-(d_J*(mag_g(i)-mag_r(i))+e_J*(mag_r(i)-mag_i(i))+f_J)
model_Ks(i)=mag_g(i)-(d_Ks*(mag_g(i)-mag_r(i))+e_Ks*(mag_r(i)-mag_i(i))+f_Ks)
model_H(i)=mag_g(i)-(d_H*(mag_g(i)-mag_r(i))+e_H*(mag_r(i)-mag_i(i))+f_H)


delta_J(i)=mag_J(i)-model_J(i)
delta_Ks(i)=mag_Ks(i)-model_Ks(i)
delta_sum_J=delta_sum_J+delta_J(i)**2
delta_sum_Ks=delta_sum_Ks+delta_Ks(i)**2

delta_H(i)=mag_H(i)-model_H(i)
delta_sum_H=delta_sum_H+delta_H(i)**2

END DO

delta_sum_J=SQRT(delta_sum_J/DBLE(n))
delta_sum_Ks=SQRT(delta_sum_Ks/DBLE(n))
delta_sum_H=SQRT(delta_sum_H/DBLE(n))

WRITE(*,*) delta_sum_J,delta_sum_H,delta_sum_Ks



OPEN(50,file='mill_delta_J.txt')
DO i=1,n
WRITE(50,*) mag_J(i),delta_J(i)
END DO
 CLOSE(50)
OPEN(50,file='mill_delta_Ks.txt')
DO i=1,n
WRITE(50,*) mag_Ks(i),delta_Ks(i)
END DO
 CLOSE(50)
 OPEN(50,file='mill_delta_H.txt')
DO i=1,n
WRITE(50,*) mag_H(i),delta_H(i)
END DO
 CLOSE(50)
 
 
 
 
 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,n

help_x=NINT((-mag_J(i)-12.D0)*10.D0)
help_y=NINT((delta_J(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_J_bilirmill_map.txt')
DO i=0,150
DO ii=0,150
help_x2=-(DBLE(i)/10.D0)-12.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,n

help_x=NINT((-mag_Ks(i)-12.D0)*10.D0)
help_y=NINT((delta_Ks(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_Ks_bilirmill_map.txt')
DO i=0,150
DO ii=0,150
help_x2=-(DBLE(i)/10.D0)-12.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,n

help_x=NINT((-mag_H(i)-12.D0)*10.D0)
help_y=NINT((delta_H(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_H_bilirmill_map.txt')
DO i=0,150
DO ii=0,150
help_x2=-(DBLE(i)/10.D0)-12.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 

!preparing matrizes
DO i=1,3 
Yvector_J(i)=0.D0
Yvector_Ks(i)=0.D0
Yvector_H(i)=0.D0
DO ii=1,3
Amatrix_J(i,ii)=0.D0
Amatrix_Ks(i,ii)=0.D0
Amatrix_H(i,ii)=0.D0
Bmatrix1_J(i,ii)=0.D0
Bmatrix2_J(i,ii)=0.D0
Bmatrix3_J(i,ii)=0.D0
Bmatrix1_Ks(i,ii)=0.D0
Bmatrix2_Ks(i,ii)=0.D0
Bmatrix3_Ks(i,ii)=0.D0
Bmatrix1_H(i,ii)=0.D0
Bmatrix2_H(i,ii)=0.D0
Bmatrix3_H(i,ii)=0.D0
END DO
END DO

!calculating basic matrix und vector
DO i=1,n


Amatrix_J(1,1)=Amatrix_J(1,1)+(mag_g(i)-mag_r(i))*(mag_g(i)-mag_r(i))
Amatrix_J(2,2)=Amatrix_J(2,2)+(mag_r(i)-mag_i(i))*(mag_r(i)-mag_i(i))
Amatrix_J(3,3)=Amatrix_J(3,3)+1.D0
Amatrix_J(1,2)=Amatrix_J(1,2)+(mag_g(i)-mag_r(i))*(mag_r(i)-mag_i(i))
Amatrix_J(2,1)=Amatrix_J(2,1)+(mag_g(i)-mag_r(i))*(mag_r(i)-mag_i(i))
Amatrix_J(1,3)=Amatrix_J(1,3)+(mag_g(i)-mag_r(i))
Amatrix_J(3,1)=Amatrix_J(3,1)+(mag_g(i)-mag_r(i))
Amatrix_J(2,3)=Amatrix_J(2,3)+(mag_r(i)-mag_i(i))
Amatrix_J(3,2)=Amatrix_J(3,2)+(mag_r(i)-mag_i(i))
Yvector_J(1)=Yvector_J(1)+(mag_g(i)-mag_J(i))*(mag_g(i)-mag_r(i))
Yvector_J(2)=Yvector_J(2)+(mag_g(i)-mag_J(i))*(mag_r(i)-mag_i(i))
Yvector_J(3)=Yvector_J(3)+(mag_g(i)-mag_J(i))

Amatrix_Ks(1,1)=Amatrix_Ks(1,1)+(mag_g(i)-mag_r(i))*(mag_g(i)-mag_r(i))
Amatrix_Ks(2,2)=Amatrix_Ks(2,2)+(mag_r(i)-mag_i(i))*(mag_r(i)-mag_i(i))
Amatrix_Ks(3,3)=Amatrix_Ks(3,3)+1.D0
Amatrix_Ks(1,2)=Amatrix_Ks(1,2)+(mag_g(i)-mag_r(i))*(mag_r(i)-mag_i(i))
Amatrix_Ks(2,1)=Amatrix_Ks(2,1)+(mag_g(i)-mag_r(i))*(mag_r(i)-mag_i(i))
Amatrix_Ks(1,3)=Amatrix_Ks(1,3)+(mag_g(i)-mag_r(i))
Amatrix_Ks(3,1)=Amatrix_Ks(3,1)+(mag_g(i)-mag_r(i))
Amatrix_Ks(2,3)=Amatrix_Ks(2,3)+(mag_r(i)-mag_i(i))
Amatrix_Ks(3,2)=Amatrix_Ks(3,2)+(mag_r(i)-mag_i(i))
Yvector_Ks(1)=Yvector_Ks(1)+(mag_g(i)-mag_Ks(i))*(mag_g(i)-mag_r(i))
Yvector_Ks(2)=Yvector_Ks(2)+(mag_g(i)-mag_Ks(i))*(mag_r(i)-mag_i(i))
Yvector_Ks(3)=Yvector_Ks(3)+(mag_g(i)-mag_Ks(i))

Amatrix_H(1,1)=Amatrix_H(1,1)+(mag_g(i)-mag_r(i))*(mag_g(i)-mag_r(i))
Amatrix_H(2,2)=Amatrix_H(2,2)+(mag_r(i)-mag_i(i))*(mag_r(i)-mag_i(i))
Amatrix_H(3,3)=Amatrix_H(3,3)+1.D0
Amatrix_H(1,2)=Amatrix_H(1,2)+(mag_g(i)-mag_r(i))*(mag_r(i)-mag_i(i))
Amatrix_H(2,1)=Amatrix_H(2,1)+(mag_g(i)-mag_r(i))*(mag_r(i)-mag_i(i))
Amatrix_H(1,3)=Amatrix_H(1,3)+(mag_g(i)-mag_r(i))
Amatrix_H(3,1)=Amatrix_H(3,1)+(mag_g(i)-mag_r(i))
Amatrix_H(2,3)=Amatrix_H(2,3)+(mag_r(i)-mag_i(i))
Amatrix_H(3,2)=Amatrix_H(3,2)+(mag_r(i)-mag_i(i))
Yvector_H(1)=Yvector_H(1)+(mag_g(i)-mag_H(i))*(mag_g(i)-mag_r(i))
Yvector_H(2)=Yvector_H(2)+(mag_g(i)-mag_H(i))*(mag_r(i)-mag_i(i))
Yvector_H(3)=Yvector_H(3)+(mag_g(i)-mag_H(i))

END DO

!get determinat of Amatrix

detA_J=Amatrix_J(1,1)*Amatrix_J(2,2)*Amatrix_J(3,3)
detA_J=detA_J+Amatrix_J(2,1)*Amatrix_J(3,2)*Amatrix_J(1,3)
detA_J=detA_J+Amatrix_J(1,2)*Amatrix_J(2,3)*Amatrix_J(3,1)
detA_J=detA_J-Amatrix_J(1,3)*Amatrix_J(2,2)*Amatrix_J(3,1)
detA_J=detA_J-Amatrix_J(2,1)*Amatrix_J(1,2)*Amatrix_J(3,3)
detA_J=detA_J-Amatrix_J(2,3)*Amatrix_J(3,2)*Amatrix_J(1,1)

detA_Ks=Amatrix_Ks(1,1)*Amatrix_Ks(2,2)*Amatrix_Ks(3,3)
detA_Ks=detA_Ks+Amatrix_Ks(2,1)*Amatrix_Ks(3,2)*Amatrix_Ks(1,3)
detA_Ks=detA_Ks+Amatrix_Ks(1,2)*Amatrix_Ks(2,3)*Amatrix_Ks(3,1)
detA_Ks=detA_Ks-Amatrix_Ks(1,3)*Amatrix_Ks(2,2)*Amatrix_Ks(3,1)
detA_Ks=detA_Ks-Amatrix_Ks(2,1)*Amatrix_Ks(1,2)*Amatrix_Ks(3,3)
detA_Ks=detA_Ks-Amatrix_Ks(2,3)*Amatrix_Ks(3,2)*Amatrix_Ks(1,1)

detA_H=Amatrix_H(1,1)*Amatrix_H(2,2)*Amatrix_H(3,3)
detA_H=detA_H+Amatrix_H(2,1)*Amatrix_H(3,2)*Amatrix_H(1,3)
detA_H=detA_H+Amatrix_H(1,2)*Amatrix_H(2,3)*Amatrix_H(3,1)
detA_H=detA_H-Amatrix_H(1,3)*Amatrix_H(2,2)*Amatrix_H(3,1)
detA_H=detA_H-Amatrix_H(2,1)*Amatrix_H(1,2)*Amatrix_H(3,3)
detA_H=detA_H-Amatrix_H(2,3)*Amatrix_H(3,2)*Amatrix_H(1,1)

!calculate B matrizes
DO i=1,3
DO ii=1,3

Bmatrix1_J(i,ii)=Amatrix_J(i,ii)
Bmatrix2_J(i,ii)=Amatrix_J(i,ii)
Bmatrix3_J(i,ii)=Amatrix_J(i,ii)

Bmatrix1_Ks(i,ii)=Amatrix_Ks(i,ii)
Bmatrix2_Ks(i,ii)=Amatrix_Ks(i,ii)
Bmatrix3_Ks(i,ii)=Amatrix_Ks(i,ii)

Bmatrix1_H(i,ii)=Amatrix_H(i,ii)
Bmatrix2_H(i,ii)=Amatrix_H(i,ii)
Bmatrix3_H(i,ii)=Amatrix_H(i,ii)

END DO
END DO

DO i=1,3

Bmatrix1_J(1,i)=Yvector_J(i)
Bmatrix2_J(2,i)=Yvector_J(i)
Bmatrix3_J(3,i)=Yvector_J(i)

Bmatrix1_Ks(1,i)=Yvector_Ks(i)
Bmatrix2_Ks(2,i)=Yvector_Ks(i)
Bmatrix3_Ks(3,i)=Yvector_Ks(i)

Bmatrix1_H(1,i)=Yvector_H(i)
Bmatrix2_H(2,i)=Yvector_H(i)
Bmatrix3_H(3,i)=Yvector_H(i)

END DO





! get determinat of Bmatrizes

detB1_J=Bmatrix1_J(1,1)*Bmatrix1_J(2,2)*Bmatrix1_J(3,3)
detB1_J=detB1_J+Bmatrix1_J(2,1)*Bmatrix1_J(3,2)*Bmatrix1_J(1,3)
detB1_J=detB1_J+Bmatrix1_J(1,2)*Bmatrix1_J(2,3)*Bmatrix1_J(3,1)
detB1_J=detB1_J-Bmatrix1_J(1,3)*Bmatrix1_J(2,2)*Bmatrix1_J(3,1)
detB1_J=detB1_J-Bmatrix1_J(2,1)*Bmatrix1_J(1,2)*Bmatrix1_J(3,3)
detB1_J=detB1_J-Bmatrix1_J(2,3)*Bmatrix1_J(3,2)*Bmatrix1_J(1,1)
detB2_J=Bmatrix2_J(1,1)*Bmatrix2_J(2,2)*Bmatrix2_J(3,3)
detB2_J=detB2_J+Bmatrix2_J(2,1)*Bmatrix2_J(3,2)*Bmatrix2_J(1,3)
detB2_J=detB2_J+Bmatrix2_J(1,2)*Bmatrix2_J(2,3)*Bmatrix2_J(3,1)
detB2_J=detB2_J-Bmatrix2_J(1,3)*Bmatrix2_J(2,2)*Bmatrix2_J(3,1)
detB2_J=detB2_J-Bmatrix2_J(2,1)*Bmatrix2_J(1,2)*Bmatrix2_J(3,3)
detB2_J=detB2_J-Bmatrix2_J(2,3)*Bmatrix2_J(3,2)*Bmatrix2_J(1,1)
detB3_J=Bmatrix3_J(1,1)*Bmatrix3_J(2,2)*Bmatrix3_J(3,3)
detB3_J=detB3_J+Bmatrix3_J(2,1)*Bmatrix3_J(3,2)*Bmatrix3_J(1,3)
detB3_J=detB3_J+Bmatrix3_J(1,2)*Bmatrix3_J(2,3)*Bmatrix3_J(3,1)
detB3_J=detB3_J-Bmatrix3_J(1,3)*Bmatrix3_J(2,2)*Bmatrix3_J(3,1)
detB3_J=detB3_J-Bmatrix3_J(2,1)*Bmatrix3_J(1,2)*Bmatrix3_J(3,3)
detB3_J=detB3_J-Bmatrix3_J(2,3)*Bmatrix3_J(3,2)*Bmatrix3_J(1,1)

detB1_Ks=Bmatrix1_Ks(1,1)*Bmatrix1_Ks(2,2)*Bmatrix1_Ks(3,3)
detB1_Ks=detB1_Ks+Bmatrix1_Ks(2,1)*Bmatrix1_Ks(3,2)*Bmatrix1_Ks(1,3)
detB1_Ks=detB1_Ks+Bmatrix1_Ks(1,2)*Bmatrix1_Ks(2,3)*Bmatrix1_Ks(3,1)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(1,3)*Bmatrix1_Ks(2,2)*Bmatrix1_Ks(3,1)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(2,1)*Bmatrix1_Ks(1,2)*Bmatrix1_Ks(3,3)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(2,3)*Bmatrix1_Ks(3,2)*Bmatrix1_Ks(1,1)
detB2_Ks=Bmatrix2_Ks(1,1)*Bmatrix2_Ks(2,2)*Bmatrix2_Ks(3,3)
detB2_Ks=detB2_Ks+Bmatrix2_Ks(2,1)*Bmatrix2_Ks(3,2)*Bmatrix2_Ks(1,3)
detB2_Ks=detB2_Ks+Bmatrix2_Ks(1,2)*Bmatrix2_Ks(2,3)*Bmatrix2_Ks(3,1)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(1,3)*Bmatrix2_Ks(2,2)*Bmatrix2_Ks(3,1)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(2,1)*Bmatrix2_Ks(1,2)*Bmatrix2_Ks(3,3)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(2,3)*Bmatrix2_Ks(3,2)*Bmatrix2_Ks(1,1)
detB3_Ks=Bmatrix3_Ks(1,1)*Bmatrix3_Ks(2,2)*Bmatrix3_Ks(3,3)
detB3_Ks=detB3_Ks+Bmatrix3_Ks(2,1)*Bmatrix3_Ks(3,2)*Bmatrix3_Ks(1,3)
detB3_Ks=detB3_Ks+Bmatrix3_Ks(1,2)*Bmatrix3_Ks(2,3)*Bmatrix3_Ks(3,1)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(1,3)*Bmatrix3_Ks(2,2)*Bmatrix3_Ks(3,1)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(2,1)*Bmatrix3_Ks(1,2)*Bmatrix3_Ks(3,3)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(2,3)*Bmatrix3_Ks(3,2)*Bmatrix3_Ks(1,1)

detB1_H=Bmatrix1_H(1,1)*Bmatrix1_H(2,2)*Bmatrix1_H(3,3)
detB1_H=detB1_H+Bmatrix1_H(2,1)*Bmatrix1_H(3,2)*Bmatrix1_H(1,3)
detB1_H=detB1_H+Bmatrix1_H(1,2)*Bmatrix1_H(2,3)*Bmatrix1_H(3,1)
detB1_H=detB1_H-Bmatrix1_H(1,3)*Bmatrix1_H(2,2)*Bmatrix1_H(3,1)
detB1_H=detB1_H-Bmatrix1_H(2,1)*Bmatrix1_H(1,2)*Bmatrix1_H(3,3)
detB1_H=detB1_H-Bmatrix1_H(2,3)*Bmatrix1_H(3,2)*Bmatrix1_H(1,1)
detB2_H=Bmatrix2_H(1,1)*Bmatrix2_H(2,2)*Bmatrix2_H(3,3)
detB2_H=detB2_H+Bmatrix2_H(2,1)*Bmatrix2_H(3,2)*Bmatrix2_H(1,3)
detB2_H=detB2_H+Bmatrix2_H(1,2)*Bmatrix2_H(2,3)*Bmatrix2_H(3,1)
detB2_H=detB2_H-Bmatrix2_H(1,3)*Bmatrix2_H(2,2)*Bmatrix2_H(3,1)
detB2_H=detB2_H-Bmatrix2_H(2,1)*Bmatrix2_H(1,2)*Bmatrix2_H(3,3)
detB2_H=detB2_H-Bmatrix2_H(2,3)*Bmatrix2_H(3,2)*Bmatrix2_H(1,1)
detB3_H=Bmatrix3_H(1,1)*Bmatrix3_H(2,2)*Bmatrix3_H(3,3)
detB3_H=detB3_H+Bmatrix3_H(2,1)*Bmatrix3_H(3,2)*Bmatrix3_H(1,3)
detB3_H=detB3_H+Bmatrix3_H(1,2)*Bmatrix3_H(2,3)*Bmatrix3_H(3,1)
detB3_H=detB3_H-Bmatrix3_H(1,3)*Bmatrix3_H(2,2)*Bmatrix3_H(3,1)
detB3_H=detB3_H-Bmatrix3_H(2,1)*Bmatrix3_H(1,2)*Bmatrix3_H(3,3)
detB3_H=detB3_H-Bmatrix3_H(2,3)*Bmatrix3_H(3,2)*Bmatrix3_H(1,1)

! get solutions

 a_J=detB1_J/detA_J
 b_J=detB2_J/detA_J
 c_J=detB3_J/detA_J

 a_Ks=detB1_Ks/detA_Ks
 b_Ks=detB2_Ks/detA_Ks
 c_Ks=detB3_Ks/detA_Ks

 a_H=detB1_H/detA_H
 b_H=detB2_H/detA_H
 c_H=detB3_H/detA_H




DO i=1,n


otherside_J(i)=a_J*(mag_g(i)-mag_r(i))+b_J*(mag_r(i)-mag_i(i))+c_J
otherside_Ks(i)=a_Ks*(mag_g(i)-mag_r(i))+b_Ks*(mag_r(i)-mag_i(i))+c_Ks
otherside_H(i)=a_H*(mag_g(i)-mag_r(i))+b_H*(mag_r(i)-mag_i(i))+c_H

END DO

! propagation of error
!root mean square
rms_J=0.D0
rms_Ks=0.D0
rms_H=0.D0

DO i=1,n


deviation_J(i)=otherside_J(i)-(mag_g(i)-mag_J(i))
deviation_Ks(i)=otherside_Ks(i)-(mag_g(i)-mag_Ks(i))
deviation_H(i)=otherside_H(i)-(mag_g(i)-mag_H(i))

rms_J=rms_J+((deviation_J(i))**2)
rms_Ks=rms_Ks+((deviation_Ks(i))**2)
rms_H=rms_H+((deviation_H(i))**2)

END DO


rms_J=SQRT(rms_J/DBLE(n))
rms_Ks=SQRT(rms_Ks/DBLE(n))
rms_H=SQRT(rms_H/DBLE(n))



! error of fitting parameters


 a_err_J=(Amatrix_J(2,2)*Amatrix_J(3,3)-Amatrix_J(2,3)*Amatrix_J(3,2))/detA_J
 a_err_J=SQRT(a_err_J)
 b_err_J=(Amatrix_J(1,1)*Amatrix_J(3,3)-Amatrix_J(1,3)*Amatrix_J(3,1))/detA_J
 b_err_J=SQRT(b_err_J)
 c_err_J=(Amatrix_J(1,1)*Amatrix_J(2,2)-Amatrix_J(1,2)*Amatrix_J(2,1))/detA_J
 c_err_J=SQRT(c_err_J)

 a_err_Ks=(Amatrix_Ks(2,2)*Amatrix_Ks(3,3)-Amatrix_Ks(2,3)*Amatrix_Ks(3,2))/detA_Ks
 a_err_Ks=SQRT(a_err_Ks)
 b_err_Ks=(Amatrix_Ks(1,1)*Amatrix_Ks(3,3)-Amatrix_Ks(1,3)*Amatrix_Ks(3,1))/detA_Ks
 b_err_Ks=SQRT(b_err_Ks)
 c_err_Ks=(Amatrix_Ks(1,1)*Amatrix_Ks(2,2)-Amatrix_Ks(1,2)*Amatrix_Ks(2,1))/detA_Ks
 c_err_Ks=SQRT(c_err_Ks)

 a_err_H=(Amatrix_H(2,2)*Amatrix_H(3,3)-Amatrix_H(2,3)*Amatrix_H(3,2))/detA_H
 a_err_H=SQRT(a_err_H)
 b_err_H=(Amatrix_H(1,1)*Amatrix_H(3,3)-Amatrix_H(1,3)*Amatrix_H(3,1))/detA_H
 b_err_H=SQRT(b_err_H)
 c_err_H=(Amatrix_H(1,1)*Amatrix_H(2,2)-Amatrix_H(1,2)*Amatrix_H(2,1))/detA_H
 c_err_H=SQRT(c_err_H)
 
 
 OPEN(50,file='dev_J.txt')
DO i=1,n
WRITE(50,*) mag_J(i),deviation_J(i)
END DO
 CLOSE(50)
OPEN(50,file='dev_Ks.txt')
DO i=1,n
WRITE(50,*) mag_Ks(i),deviation_Ks(i)
END DO
 CLOSE(50)
 OPEN(50,file='dev_H.txt')
DO i=1,n
WRITE(50,*) mag_H(i),deviation_H(i)
END DO
 CLOSE(50)
 
 
      WRITE(*,*) '------------------'    
 WRITE(*,*) 'J'
  WRITE(*,*) 'd:',a_J,'±',a_err_J
  WRITE(*,*) 'e:',b_J,'±',b_err_J
    WRITE(*,*) 'f:',c_J,'±',c_err_J
     WRITE(*,*) 'RMS:',rms_J
     WRITE(*,*) '------------------'    
      WRITE(*,*) 'H'
  WRITE(*,*) 'd:',a_H,'±',a_err_H
  WRITE(*,*) 'e:',b_H,'±',b_err_H
    WRITE(*,*) 'f:',c_H,'±',c_err_H
   WRITE(*,*) 'RMS:',rms_H
         WRITE(*,*) '------------------'    
 WRITE(*,*) 'Ks'
  WRITE(*,*) 'd:',a_Ks,'±',a_err_Ks
  WRITE(*,*) 'e:',b_Ks,'±',b_err_Ks
    WRITE(*,*) 'f:',c_Ks,'±',c_err_Ks
   WRITE(*,*) 'RMS:',rms_Ks
         WRITE(*,*) '------------------'    
 
 dummy_coeffient=0.D0
  OPEN(52,file='2mass_sdss_transformation2.txt')
WRITE(52,*) 'g-J   g-H   g-Ks'
WRITE(52,*) a_J,a_H,a_Ks
WRITE(52,*) b_J,b_H,b_Ks
WRITE(52,*) c_J,c_H,c_Ks
WRITE(52,*) dummy_coeffient,dummy_coeffient,dummy_coeffient
WRITE(52,*) rms_J,rms_H,rms_Ks
WRITE(52,*) 
WRITE(52,*) 'fit on Milleniums-Simulation'
 CLOSE(52)
 
 a_J_mill=a_J
 a_H_mill=a_H
 a_Ks_mill=a_Ks
 b_J_mill=b_J
 b_H_mill=b_H
 b_Ks_mill=b_Ks
 c_J_mill=c_J
 c_H_mill=c_H
 c_Ks_mill=c_Ks
 
 
 
 
 deallocate(otherside_J)
deallocate(otherside_Ks)
deallocate(deviation_J)
deallocate(deviation_Ks)
  deallocate(otherside_H)
deallocate(deviation_H)
 
 
 ! get length of file
OPEN(50,file='2MRS.txt')
io_err=0
n_2mass=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_2mass=n_2mass+1
END DO
 CLOSE(50)
n_2mass=n_2mass-3 ! remove header

! get length of file
OPEN(50,file='SDSS_DR12.txt')
io_err=0
n_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_sdss=n_sdss+1
END DO
 CLOSE(50)
n_sdss=n_sdss-2 ! remove header







allocate(sdss_ra(1:n_sdss))
allocate(sdss_dec(1:n_sdss))
allocate(sdss_gl(1:n_sdss))
allocate(sdss_gb(1:n_sdss))
allocate(sdss_z(1:n_sdss))
allocate(cModelMag_g(1:n_sdss))
allocate(cModelMag_r(1:n_sdss))
allocate(cModelMag_i(1:n_sdss))
allocate(cModelMagErr_g(1:n_sdss))
allocate(cModelMagErr_r(1:n_sdss))
allocate(cModelMagErr_i(1:n_sdss))
allocate(extinction_g(1:n_sdss))
allocate(extinction_r(1:n_sdss))
allocate(extinction_i(1:n_sdss))
allocate(sdss_active(1:n_sdss))


allocate(twomass_ra(1:n_2mass))
allocate(twomass_dec(1:n_2mass))
allocate(twomass_gl(1:n_2mass))
allocate(twomass_gb(1:n_2mass))
allocate(Ktmag(1:n_2mass))
allocate(Jtmag(1:n_2mass))
allocate(Kmagerr(1:n_2mass))
allocate(Jmagerr(1:n_2mass))
allocate(E_BV(1:n_2mass))
allocate(twomass_cz(1:n_2mass))
allocate(twomass_ecz(1:n_2mass))
allocate(twomass_active(1:n_2mass))
 
allocate(Htmag(1:n_2mass))
allocate(Hmagerr(1:n_2mass))
 
allocate(K_g(1:n_sdss))
allocate(K_r(1:n_sdss))
allocate(K_i(1:n_sdss))
allocate(z_sdss(1:n_sdss))
allocate(magapp_g(1:n_sdss))
allocate(magapp_r(1:n_sdss))
allocate(magapp_i(1:n_sdss))

allocate(K_J(1:n_2mass))
allocate(K_H(1:n_2mass))
allocate(K_Ks(1:n_2mass))
allocate(z_2mass(1:n_2mass))
allocate(magapp_J(1:n_2mass))
allocate(magapp_Ks(1:n_2mass))
allocate(magapp_H(1:n_2mass))

allocate(m_J(1:n_2mass))
allocate(m_H(1:n_2mass))
allocate(m_Ks(1:n_2mass))
allocate(m_g(1:n_2mass))
allocate(m_r(1:n_2mass))
allocate(m_i(1:n_2mass))
allocate(otherside_J(1:n_2mass))
allocate(otherside_Ks(1:n_2mass))
allocate(deviation_J(1:n_2mass))
allocate(deviation_Ks(1:n_2mass))
allocate(otherside_H(1:n_2mass))
allocate(deviation_H(1:n_2mass))
 
allocate(leftside_J(1:n_2mass))
allocate(leftside_H(1:n_2mass))
allocate(leftside_Ks(1:n_2mass))
allocate(rightside_J(1:n_2mass))
allocate(rightside_H(1:n_2mass))
allocate(rightside_Ks(1:n_2mass))

! read file observational data
OPEN(50,file='SDSS_DR12.txt')
READ(50,*)
DO i=1,n_sdss

READ(50,*) dummy_id,sdss_ra(i),sdss_dec(i),sdss_gb(i),sdss_gl(i),sdss_z(i),&
 cModelMag_g(i),cModelMag_r(i),cModelMag_i(i),cModelMagErr_g(i),cModelMagErr_r(i),cModelMagErr_i(i),&
 extinction_g(i),extinction_r(i),extinction_i(i)

 

sdss_active(i)=.TRUE.

END DO

 CLOSE(50)

 
 DO i=1,n_sdss
  IF (cModelMag_g(i)<5.D0) THEN
sdss_active(i)=.FALSE.
  END IF
    IF (cModelMag_r(i)<5.D0) THEN
sdss_active(i)=.FALSE.
  END IF
    IF (cModelMag_i(i)<5.D0) THEN
sdss_active(i)=.FALSE.
  END IF
  IF (cModelMag_g(i)>25.D0) THEN
sdss_active(i)=.FALSE.
  END IF
    IF (cModelMag_r(i)>25.D0) THEN
sdss_active(i)=.FALSE.
  END IF
    IF (cModelMag_i(i)>25.D0) THEN
sdss_active(i)=.FALSE. 
  END IF
 END DO

WRITE(*,*) n_sdss,'galaxies from SDSS DR12'



! read file
OPEN(50,file='2MRS_H.txt')

READ(50,*)
READ(50,*)

DO i=1,n_2mass 
READ(50,*) dummy_id,dummy_id,twomass_ra(i),twomass_dec(i),&
twomass_gl(i),twomass_gb(i),Ktmag(i),Kmagerr(i),&
Htmag(i),Hmagerr(i),Jtmag(i),Jmagerr(i),&
E_BV(i),twomass_cz(i),twomass_ecz(i)

twomass_active(i)=.TRUE.

END DO

 CLOSE(50)
WRITE(*,*) n_2mass,'galaxies from 2MRS'


 
 
 !get parameters for K correction  
OPEN(52,file='K_correction_g.txt')
DO i=0,7
READ(52,*) K_coeff_g(i,0:3)
END DO
 CLOSE(52)
OPEN(52,file='K_correction_r.txt')
DO i=0,5
READ(52,*) K_coeff_r(i,0:3)
END DO
 CLOSE(52)
OPEN(52,file='K_correction_i.txt')
DO i=0,5
READ(52,*) K_coeff_i(i,0:3)
END DO
 CLOSE(52)

OPEN(52,file='K_correction_J.txt')
DO i=0,5
READ(52,*) K_coeff_J(i,0:3)
END DO
 CLOSE(52)

 OPEN(52,file='K_correction_H.txt')
DO i=0,7
READ(52,*) K_coeff_H(i,0:3)
END DO
 CLOSE(52)
 
OPEN(52,file='K_correction_Ks.txt')
DO i=0,5
READ(52,*) K_coeff_Ks(i,0:3)
END DO
 CLOSE(52)
 
 
 
 
 
!prepare for calculation of correction for CMB motion
! values for the motion of the sun relativ to the CMB
v_cmb=369.D0
l_cmb=263.99D0
b_cmb=48.26D0

! calculate sun's motion in z space relativ to the CMB
z_cmb=v_cmb/light
gal_b=b_cmb*PI/180.D0
gal_l=l_cmb*PI/180.D0
z_cmb_x=z_cmb*COS(gal_b)*COS(gal_l)
z_cmb_y=z_cmb*COS(gal_b)*SIN(gal_l)
z_cmb_z=z_cmb*SIN(gal_b)

DO i=1,n_sdss
gal_b=sdss_gb(i)*PI/180.D0
gal_l=sdss_gl(i)*PI/180.D0
z_gal_x=sdss_z(i)*COS(gal_b)*COS(gal_l)
z_gal_y=sdss_z(i)*COS(gal_b)*SIN(gal_l)
z_gal_z=sdss_z(i)*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
z_sdss(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


DO i=1,n_2mass
gal_b=twomass_gb(i)*PI/180.D0
gal_l=twomass_gl(i)*PI/180.D0
z_gal_x=twomass_cz(i)/light*COS(gal_b)*COS(gal_l)
z_gal_y=twomass_cz(i)/light*COS(gal_b)*SIN(gal_l)
z_gal_z=twomass_cz(i)/light*SIN(gal_b)
z_gal_x=((1.D0+z_gal_x)*(1.D0+z_cmb_x))-1.D0
z_gal_y=((1.D0+z_gal_y)*(1.D0+z_cmb_y))-1.D0
z_gal_z=((1.D0+z_gal_z)*(1.D0+z_cmb_z))-1.D0
z_2mass(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
END DO


! DO i=1,n_sdss
! gal_b=sdss_gb(i)*PI/180.D0
! gal_l=sdss_gl(i)*PI/180.D0
! z_gal_x=sdss_z(i)*COS(gal_b)*COS(gal_l)
! z_gal_y=sdss_z(i)*COS(gal_b)*SIN(gal_l)
! z_gal_z=sdss_z(i)*SIN(gal_b)
! z_gal_x=z_gal_x+z_cmb_x
! z_gal_y=z_gal_y+z_cmb_y
! z_gal_z=z_gal_z+z_cmb_z
! z_sdss(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
! END DO
! 
! 
! DO i=1,n_2mass
! gal_b=twomass_gb(i)*PI/180.D0
! gal_l=twomass_gl(i)*PI/180.D0
! z_gal_x=twomass_cz(i)/light*COS(gal_b)*COS(gal_l)
! z_gal_y=twomass_cz(i)/light*COS(gal_b)*SIN(gal_l)
! z_gal_z=twomass_cz(i)/light*SIN(gal_b)
! z_gal_x=z_gal_x+z_cmb_x
! z_gal_y=z_gal_y+z_cmb_y
! z_gal_z=z_gal_z+z_cmb_z
! z_2mass(i)=z_gal_x*COS(gal_b)*COS(gal_l)+z_gal_y*COS(gal_b)*SIN(gal_l)+z_gal_z*SIN(gal_b)
! END DO
! 

 
 
 
!correct for galactic extinction
DO i=1,n_sdss
magapp_g(i)=cModelMag_g(i)-extinction_g(i)
magapp_r(i)=cModelMag_r(i)-extinction_r(i)
magapp_i(i)=cModelMag_i(i)-extinction_i(i)
END DO


DO i=1,n_2mass
magapp_J(i)=Jtmag(i)
magapp_H(i)=Htmag(i)
magapp_Ks(i)=Ktmag(i)
END DO

!K-correction
DO i=1,n_sdss

K_g(i)=0.D0
DO ii=0,7
DO iii=0,3
K_g(i)=K_g(i)+K_coeff_g(ii,iii)*(sdss_z(i)**ii)*((magapp_g(i)-magapp_r(i))**iii)
END DO
END DO

K_r(i)=0.D0
DO ii=0,5
DO iii=0,3
K_r(i)=K_r(i)+K_coeff_r(ii,iii)*(sdss_z(i)**ii)*((magapp_g(i)-magapp_r(i))**iii)
END DO
END DO

K_i(i)=0.D0
DO ii=0,5
DO iii=0,3
K_i(i)=K_i(i)+K_coeff_i(ii,iii)*(sdss_z(i)**ii)*((magapp_g(i)-magapp_i(i))**iii)
END DO
END DO

END DO 



!K-correction
DO i=1,n_2mass

K_J(i)=0.D0
DO ii=0,5
DO iii=0,3
K_J(i)=K_J(i)+K_coeff_J(ii,iii)*((twomass_cz(i)/light)**ii)*((magapp_J(i)-magapp_Ks(i))**iii)
END DO
END DO

K_H(i)=0.D0
DO ii=0,7
DO iii=0,3
K_H(i)=K_H(i)+K_coeff_H(ii,iii)*((twomass_cz(i)/light)**ii)*((magapp_H(i)-magapp_Ks(i))**iii)
END DO
END DO

K_Ks(i)=0.D0
DO ii=0,5
DO iii=0,3
K_Ks(i)=K_Ks(i)+K_coeff_Ks(ii,iii)*((twomass_cz(i)/light)**ii)*((magapp_J(i)-magapp_Ks(i))**iii)
END DO
END DO

END DO


! corrected apparent magnitudes
DO i=1,n_sdss

magapp_g(i)=magapp_g(i)-K_g(i) 
magapp_r(i)=magapp_r(i)-K_r(i) 
magapp_i(i)=magapp_i(i)-K_i(i) 



END DO



! corrected apparent magnitudes
DO i=1,n_2mass

magapp_J(i)=magapp_J(i)-K_J(i) 
magapp_H(i)=magapp_H(i)-K_H(i) 
magapp_Ks(i)=magapp_Ks(i)-K_Ks(i) 

END DO

 
 
 
 
 !Crossmatch objects 
 hcount=0
 DO i=1,n_sdss
 IF (sdss_active(i)) THEN
!  WRITE(*,*) i
 DO ii=1,n_2mass
 IF (twomass_active(ii)) THEN 
 
 
 
 
 
delta_z=z_2mass(ii)-z_sdss(i)
IF (delta_z<0.D0) THEN
delta_z=-delta_z
END IF
IF (delta_z<1.D-3) THEN


angular_sep=COS(sdss_dec(i)*PI/180.D0)*COS(twomass_dec(ii)*PI/180.D0)*COS((twomass_ra(ii)-sdss_ra(i))*PI/180.D0)
angular_sep=angular_sep+SIN(sdss_dec(i)*PI/180.D0)*SIN(twomass_dec(ii)*PI/180.D0)
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


IF (angular_sep<fivearcsec) THEN
 
 
 
 
 
!  IF ((z_sdss(i)-z_2mass(ii))<(1.D-2)) THEN
!  IF ((z_sdss(i)-z_2mass(ii))>(-1.D-2)) THEN 
!  
! angular_sep=COS(sdss_dec(i)*PI/180.D0)*COS(twomass_dec(ii)*PI/180.D0)*COS((twomass_ra(ii)-sdss_ra(i))*PI/180.D0)
! angular_sep=angular_sep+SIN(sdss_dec(i)*PI/180.D0)*SIN(twomass_dec(ii)*PI/180.D0)
! angular_sep=ACOS(angular_sep)*180.D0/PI
! 
! IF (angular_sep<0.D0) THEN
! angular_sep=-angular_sep
! END IF
!  
!  
!  IF (angular_sep<(fivearcsec)) THEN
!  
  hcount=hcount+1
 sdss_active(i)=.FALSE.
 twomass_active(ii)=.FALSE.
 
m_J(hcount)=magapp_J(ii)
m_H(hcount)=magapp_H(ii)
m_Ks(hcount)=magapp_Ks(ii)
m_g(hcount)=magapp_g(i)
m_r(hcount)=magapp_r(i)
m_i(hcount)=magapp_i(i)

 
!  END IF
 
 END IF
 END IF
 
 END IF
 END DO
 END IF
 END DO
 
 

 
 
 
 
 
 
  WRITE(*,*) '----------------------------------------------------------------'
 WRITE(*,*) hcount,'galaxies found that are in both catalogues'
 
 
 
 
 
 
 
 
 
 
 
 

!preparing matrizes
DO i=1,3 
Yvector_J(i)=0.D0
Yvector_H(i)=0.D0
Yvector_Ks(i)=0.D0
DO ii=1,3
Amatrix_J(i,ii)=0.D0
Amatrix_H(i,ii)=0.D0
Amatrix_Ks(i,ii)=0.D0
Bmatrix1_J(i,ii)=0.D0
Bmatrix2_J(i,ii)=0.D0
Bmatrix3_J(i,ii)=0.D0
Bmatrix1_Ks(i,ii)=0.D0
Bmatrix2_Ks(i,ii)=0.D0
Bmatrix3_Ks(i,ii)=0.D0
Bmatrix1_H(i,ii)=0.D0
Bmatrix2_H(i,ii)=0.D0
Bmatrix3_H(i,ii)=0.D0
END DO
END DO

!calculating basic matrix und vector
DO i=1,hcount


Amatrix_J(1,1)=Amatrix_J(1,1)+(m_g(i)-m_r(i))*(m_g(i)-m_r(i))
Amatrix_J(2,2)=Amatrix_J(2,2)+(m_r(i)-m_i(i))*(m_r(i)-m_i(i))
Amatrix_J(3,3)=Amatrix_J(3,3)+1.D0
Amatrix_J(1,2)=Amatrix_J(1,2)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_J(2,1)=Amatrix_J(2,1)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_J(1,3)=Amatrix_J(1,3)+(m_g(i)-m_r(i))
Amatrix_J(3,1)=Amatrix_J(3,1)+(m_g(i)-m_r(i))
Amatrix_J(2,3)=Amatrix_J(2,3)+(m_r(i)-m_i(i))
Amatrix_J(3,2)=Amatrix_J(3,2)+(m_r(i)-m_i(i))
Yvector_J(1)=Yvector_J(1)+(m_g(i)-m_J(i))*(m_g(i)-m_r(i))
Yvector_J(2)=Yvector_J(2)+(m_g(i)-m_J(i))*(m_r(i)-m_i(i))
Yvector_J(3)=Yvector_J(3)+(m_g(i)-m_J(i))

Amatrix_Ks(1,1)=Amatrix_Ks(1,1)+(m_g(i)-m_r(i))*(m_g(i)-m_r(i))
Amatrix_Ks(2,2)=Amatrix_Ks(2,2)+(m_r(i)-m_i(i))*(m_r(i)-m_i(i))
Amatrix_Ks(3,3)=Amatrix_Ks(3,3)+1.D0
Amatrix_Ks(1,2)=Amatrix_Ks(1,2)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_Ks(2,1)=Amatrix_Ks(2,1)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_Ks(1,3)=Amatrix_Ks(1,3)+(m_g(i)-m_r(i))
Amatrix_Ks(3,1)=Amatrix_Ks(3,1)+(m_g(i)-m_r(i))
Amatrix_Ks(2,3)=Amatrix_Ks(2,3)+(m_r(i)-m_i(i))
Amatrix_Ks(3,2)=Amatrix_Ks(3,2)+(m_r(i)-m_i(i))
Yvector_Ks(1)=Yvector_Ks(1)+(m_g(i)-m_Ks(i))*(m_g(i)-m_r(i))
Yvector_Ks(2)=Yvector_Ks(2)+(m_g(i)-m_Ks(i))*(m_r(i)-m_i(i))
Yvector_Ks(3)=Yvector_Ks(3)+(m_g(i)-m_Ks(i))


Amatrix_H(1,1)=Amatrix_H(1,1)+(m_g(i)-m_r(i))*(m_g(i)-m_r(i))
Amatrix_H(2,2)=Amatrix_H(2,2)+(m_r(i)-m_i(i))*(m_r(i)-m_i(i))
Amatrix_H(3,3)=Amatrix_H(3,3)+1.D0
Amatrix_H(1,2)=Amatrix_H(1,2)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_H(2,1)=Amatrix_H(2,1)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_H(1,3)=Amatrix_H(1,3)+(m_g(i)-m_r(i))
Amatrix_H(3,1)=Amatrix_H(3,1)+(m_g(i)-m_r(i))
Amatrix_H(2,3)=Amatrix_H(2,3)+(m_r(i)-m_i(i))
Amatrix_H(3,2)=Amatrix_H(3,2)+(m_r(i)-m_i(i))
Yvector_H(1)=Yvector_H(1)+(m_g(i)-m_H(i))*(m_g(i)-m_r(i))
Yvector_H(2)=Yvector_H(2)+(m_g(i)-m_H(i))*(m_r(i)-m_i(i))
Yvector_H(3)=Yvector_H(3)+(m_g(i)-m_H(i))

END DO

!get determinat of Amatrix

detA_J=Amatrix_J(1,1)*Amatrix_J(2,2)*Amatrix_J(3,3)
detA_J=detA_J+Amatrix_J(2,1)*Amatrix_J(3,2)*Amatrix_J(1,3)
detA_J=detA_J+Amatrix_J(1,2)*Amatrix_J(2,3)*Amatrix_J(3,1)
detA_J=detA_J-Amatrix_J(1,3)*Amatrix_J(2,2)*Amatrix_J(3,1)
detA_J=detA_J-Amatrix_J(2,1)*Amatrix_J(1,2)*Amatrix_J(3,3)
detA_J=detA_J-Amatrix_J(2,3)*Amatrix_J(3,2)*Amatrix_J(1,1)

detA_Ks=Amatrix_Ks(1,1)*Amatrix_Ks(2,2)*Amatrix_Ks(3,3)
detA_Ks=detA_Ks+Amatrix_Ks(2,1)*Amatrix_Ks(3,2)*Amatrix_Ks(1,3)
detA_Ks=detA_Ks+Amatrix_Ks(1,2)*Amatrix_Ks(2,3)*Amatrix_Ks(3,1)
detA_Ks=detA_Ks-Amatrix_Ks(1,3)*Amatrix_Ks(2,2)*Amatrix_Ks(3,1)
detA_Ks=detA_Ks-Amatrix_Ks(2,1)*Amatrix_Ks(1,2)*Amatrix_Ks(3,3)
detA_Ks=detA_Ks-Amatrix_Ks(2,3)*Amatrix_Ks(3,2)*Amatrix_Ks(1,1)

detA_H=Amatrix_H(1,1)*Amatrix_H(2,2)*Amatrix_H(3,3)
detA_H=detA_H+Amatrix_H(2,1)*Amatrix_H(3,2)*Amatrix_H(1,3)
detA_H=detA_H+Amatrix_H(1,2)*Amatrix_H(2,3)*Amatrix_H(3,1)
detA_H=detA_H-Amatrix_H(1,3)*Amatrix_H(2,2)*Amatrix_H(3,1)
detA_H=detA_H-Amatrix_H(2,1)*Amatrix_H(1,2)*Amatrix_H(3,3)
detA_H=detA_H-Amatrix_H(2,3)*Amatrix_H(3,2)*Amatrix_H(1,1)

!calculate B matrizes
DO i=1,3
DO ii=1,3

Bmatrix1_J(i,ii)=Amatrix_J(i,ii)
Bmatrix2_J(i,ii)=Amatrix_J(i,ii)
Bmatrix3_J(i,ii)=Amatrix_J(i,ii)

Bmatrix1_Ks(i,ii)=Amatrix_Ks(i,ii)
Bmatrix2_Ks(i,ii)=Amatrix_Ks(i,ii)
Bmatrix3_Ks(i,ii)=Amatrix_Ks(i,ii)

Bmatrix1_H(i,ii)=Amatrix_H(i,ii)
Bmatrix2_H(i,ii)=Amatrix_H(i,ii)
Bmatrix3_H(i,ii)=Amatrix_H(i,ii)

END DO
END DO

DO i=1,3

Bmatrix1_J(1,i)=Yvector_J(i)
Bmatrix2_J(2,i)=Yvector_J(i)
Bmatrix3_J(3,i)=Yvector_J(i)

Bmatrix1_Ks(1,i)=Yvector_Ks(i)
Bmatrix2_Ks(2,i)=Yvector_Ks(i)
Bmatrix3_Ks(3,i)=Yvector_Ks(i)

Bmatrix1_H(1,i)=Yvector_H(i)
Bmatrix2_H(2,i)=Yvector_H(i)
Bmatrix3_H(3,i)=Yvector_H(i)

END DO





! get determinat of Bmatrizes

detB1_J=Bmatrix1_J(1,1)*Bmatrix1_J(2,2)*Bmatrix1_J(3,3)
detB1_J=detB1_J+Bmatrix1_J(2,1)*Bmatrix1_J(3,2)*Bmatrix1_J(1,3)
detB1_J=detB1_J+Bmatrix1_J(1,2)*Bmatrix1_J(2,3)*Bmatrix1_J(3,1)
detB1_J=detB1_J-Bmatrix1_J(1,3)*Bmatrix1_J(2,2)*Bmatrix1_J(3,1)
detB1_J=detB1_J-Bmatrix1_J(2,1)*Bmatrix1_J(1,2)*Bmatrix1_J(3,3)
detB1_J=detB1_J-Bmatrix1_J(2,3)*Bmatrix1_J(3,2)*Bmatrix1_J(1,1)
detB2_J=Bmatrix2_J(1,1)*Bmatrix2_J(2,2)*Bmatrix2_J(3,3)
detB2_J=detB2_J+Bmatrix2_J(2,1)*Bmatrix2_J(3,2)*Bmatrix2_J(1,3)
detB2_J=detB2_J+Bmatrix2_J(1,2)*Bmatrix2_J(2,3)*Bmatrix2_J(3,1)
detB2_J=detB2_J-Bmatrix2_J(1,3)*Bmatrix2_J(2,2)*Bmatrix2_J(3,1)
detB2_J=detB2_J-Bmatrix2_J(2,1)*Bmatrix2_J(1,2)*Bmatrix2_J(3,3)
detB2_J=detB2_J-Bmatrix2_J(2,3)*Bmatrix2_J(3,2)*Bmatrix2_J(1,1)
detB3_J=Bmatrix3_J(1,1)*Bmatrix3_J(2,2)*Bmatrix3_J(3,3)
detB3_J=detB3_J+Bmatrix3_J(2,1)*Bmatrix3_J(3,2)*Bmatrix3_J(1,3)
detB3_J=detB3_J+Bmatrix3_J(1,2)*Bmatrix3_J(2,3)*Bmatrix3_J(3,1)
detB3_J=detB3_J-Bmatrix3_J(1,3)*Bmatrix3_J(2,2)*Bmatrix3_J(3,1)
detB3_J=detB3_J-Bmatrix3_J(2,1)*Bmatrix3_J(1,2)*Bmatrix3_J(3,3)
detB3_J=detB3_J-Bmatrix3_J(2,3)*Bmatrix3_J(3,2)*Bmatrix3_J(1,1)

detB1_Ks=Bmatrix1_Ks(1,1)*Bmatrix1_Ks(2,2)*Bmatrix1_Ks(3,3)
detB1_Ks=detB1_Ks+Bmatrix1_Ks(2,1)*Bmatrix1_Ks(3,2)*Bmatrix1_Ks(1,3)
detB1_Ks=detB1_Ks+Bmatrix1_Ks(1,2)*Bmatrix1_Ks(2,3)*Bmatrix1_Ks(3,1)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(1,3)*Bmatrix1_Ks(2,2)*Bmatrix1_Ks(3,1)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(2,1)*Bmatrix1_Ks(1,2)*Bmatrix1_Ks(3,3)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(2,3)*Bmatrix1_Ks(3,2)*Bmatrix1_Ks(1,1)
detB2_Ks=Bmatrix2_Ks(1,1)*Bmatrix2_Ks(2,2)*Bmatrix2_Ks(3,3)
detB2_Ks=detB2_Ks+Bmatrix2_Ks(2,1)*Bmatrix2_Ks(3,2)*Bmatrix2_Ks(1,3)
detB2_Ks=detB2_Ks+Bmatrix2_Ks(1,2)*Bmatrix2_Ks(2,3)*Bmatrix2_Ks(3,1)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(1,3)*Bmatrix2_Ks(2,2)*Bmatrix2_Ks(3,1)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(2,1)*Bmatrix2_Ks(1,2)*Bmatrix2_Ks(3,3)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(2,3)*Bmatrix2_Ks(3,2)*Bmatrix2_Ks(1,1)
detB3_Ks=Bmatrix3_Ks(1,1)*Bmatrix3_Ks(2,2)*Bmatrix3_Ks(3,3)
detB3_Ks=detB3_Ks+Bmatrix3_Ks(2,1)*Bmatrix3_Ks(3,2)*Bmatrix3_Ks(1,3)
detB3_Ks=detB3_Ks+Bmatrix3_Ks(1,2)*Bmatrix3_Ks(2,3)*Bmatrix3_Ks(3,1)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(1,3)*Bmatrix3_Ks(2,2)*Bmatrix3_Ks(3,1)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(2,1)*Bmatrix3_Ks(1,2)*Bmatrix3_Ks(3,3)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(2,3)*Bmatrix3_Ks(3,2)*Bmatrix3_Ks(1,1)


detB1_H=Bmatrix1_H(1,1)*Bmatrix1_H(2,2)*Bmatrix1_H(3,3)
detB1_H=detB1_H+Bmatrix1_H(2,1)*Bmatrix1_H(3,2)*Bmatrix1_H(1,3)
detB1_H=detB1_H+Bmatrix1_H(1,2)*Bmatrix1_H(2,3)*Bmatrix1_H(3,1)
detB1_H=detB1_H-Bmatrix1_H(1,3)*Bmatrix1_H(2,2)*Bmatrix1_H(3,1)
detB1_H=detB1_H-Bmatrix1_H(2,1)*Bmatrix1_H(1,2)*Bmatrix1_H(3,3)
detB1_H=detB1_H-Bmatrix1_H(2,3)*Bmatrix1_H(3,2)*Bmatrix1_H(1,1)
detB2_H=Bmatrix2_H(1,1)*Bmatrix2_H(2,2)*Bmatrix2_H(3,3)
detB2_H=detB2_H+Bmatrix2_H(2,1)*Bmatrix2_H(3,2)*Bmatrix2_H(1,3)
detB2_H=detB2_H+Bmatrix2_H(1,2)*Bmatrix2_H(2,3)*Bmatrix2_H(3,1)
detB2_H=detB2_H-Bmatrix2_H(1,3)*Bmatrix2_H(2,2)*Bmatrix2_H(3,1)
detB2_H=detB2_H-Bmatrix2_H(2,1)*Bmatrix2_H(1,2)*Bmatrix2_H(3,3)
detB2_H=detB2_H-Bmatrix2_H(2,3)*Bmatrix2_H(3,2)*Bmatrix2_H(1,1)
detB3_H=Bmatrix3_H(1,1)*Bmatrix3_H(2,2)*Bmatrix3_H(3,3)
detB3_H=detB3_H+Bmatrix3_H(2,1)*Bmatrix3_H(3,2)*Bmatrix3_H(1,3)
detB3_H=detB3_H+Bmatrix3_H(1,2)*Bmatrix3_H(2,3)*Bmatrix3_H(3,1)
detB3_H=detB3_H-Bmatrix3_H(1,3)*Bmatrix3_H(2,2)*Bmatrix3_H(3,1)
detB3_H=detB3_H-Bmatrix3_H(2,1)*Bmatrix3_H(1,2)*Bmatrix3_H(3,3)
detB3_H=detB3_H-Bmatrix3_H(2,3)*Bmatrix3_H(3,2)*Bmatrix3_H(1,1)

! get solutions

 a_J=detB1_J/detA_J
 b_J=detB2_J/detA_J
 c_J=detB3_J/detA_J

 a_Ks=detB1_Ks/detA_Ks
 b_Ks=detB2_Ks/detA_Ks
 c_Ks=detB3_Ks/detA_Ks


 a_H=detB1_H/detA_H
 b_H=detB2_H/detA_H
 c_H=detB3_H/detA_H



!calculate right side of fundamental plane
DO i=1,hcount


otherside_J(i)=a_J*(m_g(i)-m_r(i))+b_J*(m_r(i)-m_i(i))+c_J
otherside_Ks(i)=a_Ks*(m_g(i)-m_r(i))+b_Ks*(m_r(i)-m_i(i))+c_Ks
otherside_H(i)=a_H*(m_g(i)-m_r(i))+b_H*(m_r(i)-m_i(i))+c_H

rightside_J(i)=otherside_J(i)
rightside_Ks(i)=otherside_Ks(i)
rightside_H(i)=otherside_H(i)

leftside_J(i)=(m_g(i)-m_J(i))
leftside_Ks(i)=(m_g(i)-m_Ks(i))
leftside_H(i)=(m_g(i)-m_H(i))

END DO

! propagation of error
!root mean square
rms_J=0.D0
rms_Ks=0.D0
rms_H=0.D0

DO i=1,hcount


deviation_J(i)=otherside_J(i)-(m_g(i)-m_J(i))
deviation_Ks(i)=otherside_Ks(i)-(m_g(i)-m_Ks(i))
deviation_H(i)=otherside_H(i)-(m_g(i)-m_H(i))

rms_J=rms_J+((deviation_J(i))**2)
rms_Ks=rms_Ks+((deviation_Ks(i))**2)
rms_H=rms_H+((deviation_H(i))**2)


END DO


rms_J=SQRT(rms_J/DBLE(hcount))
rms_Ks=SQRT(rms_Ks/DBLE(hcount))
rms_H=SQRT(rms_H/DBLE(hcount))



! error of fitting parameters


 a_err_J=(Amatrix_J(2,2)*Amatrix_J(3,3)-Amatrix_J(2,3)*Amatrix_J(3,2))/detA_J
 a_err_J=SQRT(a_err_J)
 b_err_J=(Amatrix_J(1,1)*Amatrix_J(3,3)-Amatrix_J(1,3)*Amatrix_J(3,1))/detA_J
 b_err_J=SQRT(b_err_J)
 c_err_J=(Amatrix_J(1,1)*Amatrix_J(2,2)-Amatrix_J(1,2)*Amatrix_J(2,1))/detA_J
 c_err_J=SQRT(c_err_J)

 a_err_Ks=(Amatrix_Ks(2,2)*Amatrix_Ks(3,3)-Amatrix_Ks(2,3)*Amatrix_Ks(3,2))/detA_Ks
 a_err_Ks=SQRT(a_err_Ks)
 b_err_Ks=(Amatrix_Ks(1,1)*Amatrix_Ks(3,3)-Amatrix_Ks(1,3)*Amatrix_Ks(3,1))/detA_Ks
 b_err_Ks=SQRT(b_err_Ks)
 c_err_Ks=(Amatrix_Ks(1,1)*Amatrix_Ks(2,2)-Amatrix_Ks(1,2)*Amatrix_Ks(2,1))/detA_Ks
 c_err_Ks=SQRT(c_err_Ks)

 a_err_H=(Amatrix_H(2,2)*Amatrix_H(3,3)-Amatrix_H(2,3)*Amatrix_H(3,2))/detA_H
 a_err_H=SQRT(a_err_H)
 b_err_H=(Amatrix_H(1,1)*Amatrix_H(3,3)-Amatrix_H(1,3)*Amatrix_H(3,1))/detA_H
 b_err_H=SQRT(b_err_H)
 c_err_H=(Amatrix_H(1,1)*Amatrix_H(2,2)-Amatrix_H(1,2)*Amatrix_H(2,1))/detA_H
 c_err_H=SQRT(c_err_H)

 

 a_err_J=a_err_J/SQRT(DBLE(hcount))
 b_err_J=b_err_J/SQRT(DBLE(hcount))
 c_err_J=c_err_J/SQRT(DBLE(hcount))
 
 a_err_Ks=a_err_Ks/SQRT(DBLE(hcount))
 b_err_Ks=b_err_Ks/SQRT(DBLE(hcount))
 c_err_Ks=c_err_Ks/SQRT(DBLE(hcount))
 
 a_err_H=a_err_H/SQRT(DBLE(hcount))
 b_err_H=b_err_H/SQRT(DBLE(hcount))
 c_err_H=c_err_H/SQRT(DBLE(hcount))


 
 
 
 OPEN(50,file='dev_J_obs.txt')
DO i=1,hcount
WRITE(50,*) m_J(i),deviation_J(i)
END DO
 CLOSE(50)
  OPEN(50,file='dev_H_obs.txt')
DO i=1,hcount
WRITE(50,*) m_H(i),deviation_H(i)
END DO
 CLOSE(50)
OPEN(50,file='dev_Ks_obs.txt')
DO i=1,hcount
WRITE(50,*) m_Ks(i),deviation_Ks(i)
END DO
 CLOSE(50)
 
      WRITE(*,*) '------------------'    
 WRITE(*,*) 'J'
  WRITE(*,*) 'd:',a_J,'±',a_err_J
  WRITE(*,*) 'e:',b_J,'±',b_err_J
    WRITE(*,*) 'f:',c_J,'±',c_err_J
     WRITE(*,*) 'RMS:',rms_J
     WRITE(*,*) '------------------'    
      WRITE(*,*) 'H'
  WRITE(*,*) 'd:',a_H,'±',a_err_H
  WRITE(*,*) 'e:',b_H,'±',b_err_H
    WRITE(*,*) 'f:',c_H,'±',c_err_H
     WRITE(*,*) 'RMS:',rms_H
     WRITE(*,*) '------------------'    
 WRITE(*,*) 'Ks'
  WRITE(*,*) 'd:',a_Ks,'±',a_err_Ks
  WRITE(*,*) 'e:',b_Ks,'±',b_err_Ks
    WRITE(*,*) 'f:',c_Ks,'±',c_err_Ks
   WRITE(*,*) 'RMS:',rms_Ks
         WRITE(*,*) '------------------'    
 
 dummy_coeffient=0.D0
  OPEN(52,file='2mass_sdss_transformation3_1it.txt')
WRITE(52,*) 'g-J   g-H   g-Ks'
WRITE(52,*) a_J,a_H,a_Ks
WRITE(52,*) b_J,b_H,b_Ks
WRITE(52,*) c_J,c_H,c_Ks
WRITE(52,*) dummy_coeffient,dummy_coeffient,dummy_coeffient
WRITE(52,*) rms_J,rms_H,rms_Ks
WRITE(52,*) 
WRITE(52,*) 'fit on observed data ... 1 iteration'
 CLOSE(52)
 
 


ii=0
i=0

DO WHILE (i<hcount)
i=i+1


IF ((abs(deviation_J(i))>(3.D0*rms_J)).OR.(abs(deviation_J(i))>(3.D0*rms_J)).OR.(abs(deviation_J(i))>(3.D0*rms_J))) THEN

m_J(i)=m_J(hcount)
m_H(i)=m_H(hcount)
m_Ks(i)=m_Ks(hcount)
m_g(i)=m_g(hcount)
m_r(i)=m_r(hcount)
m_i(i)=m_i(hcount)
deviation_J(i)=deviation_J(hcount)
deviation_Ks(i)=deviation_Ks(hcount)
deviation_H(i)=deviation_H(hcount)

hcount=hcount-1
i=i-1
END IF


END DO



 
 
  WRITE(*,*) '----------------------------------------------------------------'
 WRITE(*,*) 'using ',hcount,'galaxies after 3-sigma clipping'
 
 
 
 
 
 
 
 
 
 
 
 

!preparing matrizes
DO i=1,3 
Yvector_J(i)=0.D0
Yvector_H(i)=0.D0
Yvector_Ks(i)=0.D0
DO ii=1,3
Amatrix_J(i,ii)=0.D0
Amatrix_H(i,ii)=0.D0
Amatrix_Ks(i,ii)=0.D0
Bmatrix1_J(i,ii)=0.D0
Bmatrix2_J(i,ii)=0.D0
Bmatrix3_J(i,ii)=0.D0
Bmatrix1_Ks(i,ii)=0.D0
Bmatrix2_Ks(i,ii)=0.D0
Bmatrix3_Ks(i,ii)=0.D0
Bmatrix1_H(i,ii)=0.D0
Bmatrix2_H(i,ii)=0.D0
Bmatrix3_H(i,ii)=0.D0
END DO
END DO

!calculating basic matrix und vector
DO i=1,hcount


Amatrix_J(1,1)=Amatrix_J(1,1)+(m_g(i)-m_r(i))*(m_g(i)-m_r(i))
Amatrix_J(2,2)=Amatrix_J(2,2)+(m_r(i)-m_i(i))*(m_r(i)-m_i(i))
Amatrix_J(3,3)=Amatrix_J(3,3)+1.D0
Amatrix_J(1,2)=Amatrix_J(1,2)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_J(2,1)=Amatrix_J(2,1)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_J(1,3)=Amatrix_J(1,3)+(m_g(i)-m_r(i))
Amatrix_J(3,1)=Amatrix_J(3,1)+(m_g(i)-m_r(i))
Amatrix_J(2,3)=Amatrix_J(2,3)+(m_r(i)-m_i(i))
Amatrix_J(3,2)=Amatrix_J(3,2)+(m_r(i)-m_i(i))
Yvector_J(1)=Yvector_J(1)+(m_g(i)-m_J(i))*(m_g(i)-m_r(i))
Yvector_J(2)=Yvector_J(2)+(m_g(i)-m_J(i))*(m_r(i)-m_i(i))
Yvector_J(3)=Yvector_J(3)+(m_g(i)-m_J(i))

Amatrix_Ks(1,1)=Amatrix_Ks(1,1)+(m_g(i)-m_r(i))*(m_g(i)-m_r(i))
Amatrix_Ks(2,2)=Amatrix_Ks(2,2)+(m_r(i)-m_i(i))*(m_r(i)-m_i(i))
Amatrix_Ks(3,3)=Amatrix_Ks(3,3)+1.D0
Amatrix_Ks(1,2)=Amatrix_Ks(1,2)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_Ks(2,1)=Amatrix_Ks(2,1)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_Ks(1,3)=Amatrix_Ks(1,3)+(m_g(i)-m_r(i))
Amatrix_Ks(3,1)=Amatrix_Ks(3,1)+(m_g(i)-m_r(i))
Amatrix_Ks(2,3)=Amatrix_Ks(2,3)+(m_r(i)-m_i(i))
Amatrix_Ks(3,2)=Amatrix_Ks(3,2)+(m_r(i)-m_i(i))
Yvector_Ks(1)=Yvector_Ks(1)+(m_g(i)-m_Ks(i))*(m_g(i)-m_r(i))
Yvector_Ks(2)=Yvector_Ks(2)+(m_g(i)-m_Ks(i))*(m_r(i)-m_i(i))
Yvector_Ks(3)=Yvector_Ks(3)+(m_g(i)-m_Ks(i))


Amatrix_H(1,1)=Amatrix_H(1,1)+(m_g(i)-m_r(i))*(m_g(i)-m_r(i))
Amatrix_H(2,2)=Amatrix_H(2,2)+(m_r(i)-m_i(i))*(m_r(i)-m_i(i))
Amatrix_H(3,3)=Amatrix_H(3,3)+1.D0
Amatrix_H(1,2)=Amatrix_H(1,2)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_H(2,1)=Amatrix_H(2,1)+(m_g(i)-m_r(i))*(m_r(i)-m_i(i))
Amatrix_H(1,3)=Amatrix_H(1,3)+(m_g(i)-m_r(i))
Amatrix_H(3,1)=Amatrix_H(3,1)+(m_g(i)-m_r(i))
Amatrix_H(2,3)=Amatrix_H(2,3)+(m_r(i)-m_i(i))
Amatrix_H(3,2)=Amatrix_H(3,2)+(m_r(i)-m_i(i))
Yvector_H(1)=Yvector_H(1)+(m_g(i)-m_H(i))*(m_g(i)-m_r(i))
Yvector_H(2)=Yvector_H(2)+(m_g(i)-m_H(i))*(m_r(i)-m_i(i))
Yvector_H(3)=Yvector_H(3)+(m_g(i)-m_H(i))

END DO

!get determinat of Amatrix

detA_J=Amatrix_J(1,1)*Amatrix_J(2,2)*Amatrix_J(3,3)
detA_J=detA_J+Amatrix_J(2,1)*Amatrix_J(3,2)*Amatrix_J(1,3)
detA_J=detA_J+Amatrix_J(1,2)*Amatrix_J(2,3)*Amatrix_J(3,1)
detA_J=detA_J-Amatrix_J(1,3)*Amatrix_J(2,2)*Amatrix_J(3,1)
detA_J=detA_J-Amatrix_J(2,1)*Amatrix_J(1,2)*Amatrix_J(3,3)
detA_J=detA_J-Amatrix_J(2,3)*Amatrix_J(3,2)*Amatrix_J(1,1)

detA_Ks=Amatrix_Ks(1,1)*Amatrix_Ks(2,2)*Amatrix_Ks(3,3)
detA_Ks=detA_Ks+Amatrix_Ks(2,1)*Amatrix_Ks(3,2)*Amatrix_Ks(1,3)
detA_Ks=detA_Ks+Amatrix_Ks(1,2)*Amatrix_Ks(2,3)*Amatrix_Ks(3,1)
detA_Ks=detA_Ks-Amatrix_Ks(1,3)*Amatrix_Ks(2,2)*Amatrix_Ks(3,1)
detA_Ks=detA_Ks-Amatrix_Ks(2,1)*Amatrix_Ks(1,2)*Amatrix_Ks(3,3)
detA_Ks=detA_Ks-Amatrix_Ks(2,3)*Amatrix_Ks(3,2)*Amatrix_Ks(1,1)

detA_H=Amatrix_H(1,1)*Amatrix_H(2,2)*Amatrix_H(3,3)
detA_H=detA_H+Amatrix_H(2,1)*Amatrix_H(3,2)*Amatrix_H(1,3)
detA_H=detA_H+Amatrix_H(1,2)*Amatrix_H(2,3)*Amatrix_H(3,1)
detA_H=detA_H-Amatrix_H(1,3)*Amatrix_H(2,2)*Amatrix_H(3,1)
detA_H=detA_H-Amatrix_H(2,1)*Amatrix_H(1,2)*Amatrix_H(3,3)
detA_H=detA_H-Amatrix_H(2,3)*Amatrix_H(3,2)*Amatrix_H(1,1)

!calculate B matrizes
DO i=1,3
DO ii=1,3

Bmatrix1_J(i,ii)=Amatrix_J(i,ii)
Bmatrix2_J(i,ii)=Amatrix_J(i,ii)
Bmatrix3_J(i,ii)=Amatrix_J(i,ii)

Bmatrix1_Ks(i,ii)=Amatrix_Ks(i,ii)
Bmatrix2_Ks(i,ii)=Amatrix_Ks(i,ii)
Bmatrix3_Ks(i,ii)=Amatrix_Ks(i,ii)

Bmatrix1_H(i,ii)=Amatrix_H(i,ii)
Bmatrix2_H(i,ii)=Amatrix_H(i,ii)
Bmatrix3_H(i,ii)=Amatrix_H(i,ii)

END DO
END DO

DO i=1,3

Bmatrix1_J(1,i)=Yvector_J(i)
Bmatrix2_J(2,i)=Yvector_J(i)
Bmatrix3_J(3,i)=Yvector_J(i)

Bmatrix1_Ks(1,i)=Yvector_Ks(i)
Bmatrix2_Ks(2,i)=Yvector_Ks(i)
Bmatrix3_Ks(3,i)=Yvector_Ks(i)

Bmatrix1_H(1,i)=Yvector_H(i)
Bmatrix2_H(2,i)=Yvector_H(i)
Bmatrix3_H(3,i)=Yvector_H(i)

END DO





! get determinat of Bmatrizes

detB1_J=Bmatrix1_J(1,1)*Bmatrix1_J(2,2)*Bmatrix1_J(3,3)
detB1_J=detB1_J+Bmatrix1_J(2,1)*Bmatrix1_J(3,2)*Bmatrix1_J(1,3)
detB1_J=detB1_J+Bmatrix1_J(1,2)*Bmatrix1_J(2,3)*Bmatrix1_J(3,1)
detB1_J=detB1_J-Bmatrix1_J(1,3)*Bmatrix1_J(2,2)*Bmatrix1_J(3,1)
detB1_J=detB1_J-Bmatrix1_J(2,1)*Bmatrix1_J(1,2)*Bmatrix1_J(3,3)
detB1_J=detB1_J-Bmatrix1_J(2,3)*Bmatrix1_J(3,2)*Bmatrix1_J(1,1)
detB2_J=Bmatrix2_J(1,1)*Bmatrix2_J(2,2)*Bmatrix2_J(3,3)
detB2_J=detB2_J+Bmatrix2_J(2,1)*Bmatrix2_J(3,2)*Bmatrix2_J(1,3)
detB2_J=detB2_J+Bmatrix2_J(1,2)*Bmatrix2_J(2,3)*Bmatrix2_J(3,1)
detB2_J=detB2_J-Bmatrix2_J(1,3)*Bmatrix2_J(2,2)*Bmatrix2_J(3,1)
detB2_J=detB2_J-Bmatrix2_J(2,1)*Bmatrix2_J(1,2)*Bmatrix2_J(3,3)
detB2_J=detB2_J-Bmatrix2_J(2,3)*Bmatrix2_J(3,2)*Bmatrix2_J(1,1)
detB3_J=Bmatrix3_J(1,1)*Bmatrix3_J(2,2)*Bmatrix3_J(3,3)
detB3_J=detB3_J+Bmatrix3_J(2,1)*Bmatrix3_J(3,2)*Bmatrix3_J(1,3)
detB3_J=detB3_J+Bmatrix3_J(1,2)*Bmatrix3_J(2,3)*Bmatrix3_J(3,1)
detB3_J=detB3_J-Bmatrix3_J(1,3)*Bmatrix3_J(2,2)*Bmatrix3_J(3,1)
detB3_J=detB3_J-Bmatrix3_J(2,1)*Bmatrix3_J(1,2)*Bmatrix3_J(3,3)
detB3_J=detB3_J-Bmatrix3_J(2,3)*Bmatrix3_J(3,2)*Bmatrix3_J(1,1)

detB1_Ks=Bmatrix1_Ks(1,1)*Bmatrix1_Ks(2,2)*Bmatrix1_Ks(3,3)
detB1_Ks=detB1_Ks+Bmatrix1_Ks(2,1)*Bmatrix1_Ks(3,2)*Bmatrix1_Ks(1,3)
detB1_Ks=detB1_Ks+Bmatrix1_Ks(1,2)*Bmatrix1_Ks(2,3)*Bmatrix1_Ks(3,1)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(1,3)*Bmatrix1_Ks(2,2)*Bmatrix1_Ks(3,1)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(2,1)*Bmatrix1_Ks(1,2)*Bmatrix1_Ks(3,3)
detB1_Ks=detB1_Ks-Bmatrix1_Ks(2,3)*Bmatrix1_Ks(3,2)*Bmatrix1_Ks(1,1)
detB2_Ks=Bmatrix2_Ks(1,1)*Bmatrix2_Ks(2,2)*Bmatrix2_Ks(3,3)
detB2_Ks=detB2_Ks+Bmatrix2_Ks(2,1)*Bmatrix2_Ks(3,2)*Bmatrix2_Ks(1,3)
detB2_Ks=detB2_Ks+Bmatrix2_Ks(1,2)*Bmatrix2_Ks(2,3)*Bmatrix2_Ks(3,1)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(1,3)*Bmatrix2_Ks(2,2)*Bmatrix2_Ks(3,1)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(2,1)*Bmatrix2_Ks(1,2)*Bmatrix2_Ks(3,3)
detB2_Ks=detB2_Ks-Bmatrix2_Ks(2,3)*Bmatrix2_Ks(3,2)*Bmatrix2_Ks(1,1)
detB3_Ks=Bmatrix3_Ks(1,1)*Bmatrix3_Ks(2,2)*Bmatrix3_Ks(3,3)
detB3_Ks=detB3_Ks+Bmatrix3_Ks(2,1)*Bmatrix3_Ks(3,2)*Bmatrix3_Ks(1,3)
detB3_Ks=detB3_Ks+Bmatrix3_Ks(1,2)*Bmatrix3_Ks(2,3)*Bmatrix3_Ks(3,1)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(1,3)*Bmatrix3_Ks(2,2)*Bmatrix3_Ks(3,1)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(2,1)*Bmatrix3_Ks(1,2)*Bmatrix3_Ks(3,3)
detB3_Ks=detB3_Ks-Bmatrix3_Ks(2,3)*Bmatrix3_Ks(3,2)*Bmatrix3_Ks(1,1)


detB1_H=Bmatrix1_H(1,1)*Bmatrix1_H(2,2)*Bmatrix1_H(3,3)
detB1_H=detB1_H+Bmatrix1_H(2,1)*Bmatrix1_H(3,2)*Bmatrix1_H(1,3)
detB1_H=detB1_H+Bmatrix1_H(1,2)*Bmatrix1_H(2,3)*Bmatrix1_H(3,1)
detB1_H=detB1_H-Bmatrix1_H(1,3)*Bmatrix1_H(2,2)*Bmatrix1_H(3,1)
detB1_H=detB1_H-Bmatrix1_H(2,1)*Bmatrix1_H(1,2)*Bmatrix1_H(3,3)
detB1_H=detB1_H-Bmatrix1_H(2,3)*Bmatrix1_H(3,2)*Bmatrix1_H(1,1)
detB2_H=Bmatrix2_H(1,1)*Bmatrix2_H(2,2)*Bmatrix2_H(3,3)
detB2_H=detB2_H+Bmatrix2_H(2,1)*Bmatrix2_H(3,2)*Bmatrix2_H(1,3)
detB2_H=detB2_H+Bmatrix2_H(1,2)*Bmatrix2_H(2,3)*Bmatrix2_H(3,1)
detB2_H=detB2_H-Bmatrix2_H(1,3)*Bmatrix2_H(2,2)*Bmatrix2_H(3,1)
detB2_H=detB2_H-Bmatrix2_H(2,1)*Bmatrix2_H(1,2)*Bmatrix2_H(3,3)
detB2_H=detB2_H-Bmatrix2_H(2,3)*Bmatrix2_H(3,2)*Bmatrix2_H(1,1)
detB3_H=Bmatrix3_H(1,1)*Bmatrix3_H(2,2)*Bmatrix3_H(3,3)
detB3_H=detB3_H+Bmatrix3_H(2,1)*Bmatrix3_H(3,2)*Bmatrix3_H(1,3)
detB3_H=detB3_H+Bmatrix3_H(1,2)*Bmatrix3_H(2,3)*Bmatrix3_H(3,1)
detB3_H=detB3_H-Bmatrix3_H(1,3)*Bmatrix3_H(2,2)*Bmatrix3_H(3,1)
detB3_H=detB3_H-Bmatrix3_H(2,1)*Bmatrix3_H(1,2)*Bmatrix3_H(3,3)
detB3_H=detB3_H-Bmatrix3_H(2,3)*Bmatrix3_H(3,2)*Bmatrix3_H(1,1)

! get solutions

 a_J=detB1_J/detA_J
 b_J=detB2_J/detA_J
 c_J=detB3_J/detA_J

 a_Ks=detB1_Ks/detA_Ks
 b_Ks=detB2_Ks/detA_Ks
 c_Ks=detB3_Ks/detA_Ks


 a_H=detB1_H/detA_H
 b_H=detB2_H/detA_H
 c_H=detB3_H/detA_H



!calculate right side of fundamental plane
DO i=1,hcount


otherside_J(i)=a_J*(m_g(i)-m_r(i))+b_J*(m_r(i)-m_i(i))+c_J
otherside_Ks(i)=a_Ks*(m_g(i)-m_r(i))+b_Ks*(m_r(i)-m_i(i))+c_Ks
otherside_H(i)=a_H*(m_g(i)-m_r(i))+b_H*(m_r(i)-m_i(i))+c_H

rightside_J(i)=otherside_J(i)
rightside_Ks(i)=otherside_Ks(i)
rightside_H(i)=otherside_H(i)

leftside_J(i)=(m_g(i)-m_J(i))
leftside_Ks(i)=(m_g(i)-m_Ks(i))
leftside_H(i)=(m_g(i)-m_H(i))

END DO

! propagation of error
!root mean square
rms_J=0.D0
rms_Ks=0.D0
rms_H=0.D0

DO i=1,hcount


deviation_J(i)=otherside_J(i)-(m_g(i)-m_J(i))
deviation_Ks(i)=otherside_Ks(i)-(m_g(i)-m_Ks(i))
deviation_H(i)=otherside_H(i)-(m_g(i)-m_H(i))

rms_J=rms_J+((deviation_J(i))**2)
rms_Ks=rms_Ks+((deviation_Ks(i))**2)
rms_H=rms_H+((deviation_H(i))**2)


END DO


rms_J=SQRT(rms_J/DBLE(hcount))
rms_Ks=SQRT(rms_Ks/DBLE(hcount))
rms_H=SQRT(rms_H/DBLE(hcount))



! error of fitting parameters


 a_err_J=(Amatrix_J(2,2)*Amatrix_J(3,3)-Amatrix_J(2,3)*Amatrix_J(3,2))/detA_J
 a_err_J=SQRT(a_err_J)
 b_err_J=(Amatrix_J(1,1)*Amatrix_J(3,3)-Amatrix_J(1,3)*Amatrix_J(3,1))/detA_J
 b_err_J=SQRT(b_err_J)
 c_err_J=(Amatrix_J(1,1)*Amatrix_J(2,2)-Amatrix_J(1,2)*Amatrix_J(2,1))/detA_J
 c_err_J=SQRT(c_err_J)

 a_err_Ks=(Amatrix_Ks(2,2)*Amatrix_Ks(3,3)-Amatrix_Ks(2,3)*Amatrix_Ks(3,2))/detA_Ks
 a_err_Ks=SQRT(a_err_Ks)
 b_err_Ks=(Amatrix_Ks(1,1)*Amatrix_Ks(3,3)-Amatrix_Ks(1,3)*Amatrix_Ks(3,1))/detA_Ks
 b_err_Ks=SQRT(b_err_Ks)
 c_err_Ks=(Amatrix_Ks(1,1)*Amatrix_Ks(2,2)-Amatrix_Ks(1,2)*Amatrix_Ks(2,1))/detA_Ks
 c_err_Ks=SQRT(c_err_Ks)

 a_err_H=(Amatrix_H(2,2)*Amatrix_H(3,3)-Amatrix_H(2,3)*Amatrix_H(3,2))/detA_H
 a_err_H=SQRT(a_err_H)
 b_err_H=(Amatrix_H(1,1)*Amatrix_H(3,3)-Amatrix_H(1,3)*Amatrix_H(3,1))/detA_H
 b_err_H=SQRT(b_err_H)
 c_err_H=(Amatrix_H(1,1)*Amatrix_H(2,2)-Amatrix_H(1,2)*Amatrix_H(2,1))/detA_H
 c_err_H=SQRT(c_err_H)

 a_err_J=a_err_J/SQRT(DBLE(hcount))
 b_err_J=b_err_J/SQRT(DBLE(hcount))
 c_err_J=c_err_J/SQRT(DBLE(hcount))
 
 a_err_Ks=a_err_Ks/SQRT(DBLE(hcount))
 b_err_Ks=b_err_Ks/SQRT(DBLE(hcount))
 c_err_Ks=c_err_Ks/SQRT(DBLE(hcount))
 
 a_err_H=a_err_H/SQRT(DBLE(hcount))
 b_err_H=b_err_H/SQRT(DBLE(hcount))
 c_err_H=c_err_H/SQRT(DBLE(hcount))
 
 OPEN(50,file='dev_J_obs.txt')
DO i=1,hcount
WRITE(50,*) m_J(i),deviation_J(i)
END DO
 CLOSE(50)
  OPEN(50,file='dev_H_obs.txt')
DO i=1,hcount
WRITE(50,*) m_H(i),deviation_H(i)
END DO
 CLOSE(50)
OPEN(50,file='dev_Ks_obs.txt')
DO i=1,hcount
WRITE(50,*) m_Ks(i),deviation_Ks(i)
END DO
 CLOSE(50)
 
      WRITE(*,*) '------------------'    
 WRITE(*,*) 'J'
  WRITE(*,*) 'd:',a_J,'±',a_err_J
  WRITE(*,*) 'e:',b_J,'±',b_err_J
    WRITE(*,*) 'f:',c_J,'±',c_err_J
     WRITE(*,*) 'RMS:',rms_J
     WRITE(*,*) '------------------'    
      WRITE(*,*) 'H'
  WRITE(*,*) 'd:',a_H,'±',a_err_H
  WRITE(*,*) 'e:',b_H,'±',b_err_H
    WRITE(*,*) 'f:',c_H,'±',c_err_H
     WRITE(*,*) 'RMS:',rms_H
     WRITE(*,*) '------------------'    
 WRITE(*,*) 'Ks'
  WRITE(*,*) 'd:',a_Ks,'±',a_err_Ks
  WRITE(*,*) 'e:',b_Ks,'±',b_err_Ks
    WRITE(*,*) 'f:',c_Ks,'±',c_err_Ks
   WRITE(*,*) 'RMS:',rms_Ks
         WRITE(*,*) '------------------'    
 
 dummy_coeffient=0.D0
  OPEN(52,file='2mass_sdss_transformation3.txt')
WRITE(52,*) 'g-J   g-H   g-Ks'
WRITE(52,*) a_J,a_H,a_Ks
WRITE(52,*) b_J,b_H,b_Ks
WRITE(52,*) c_J,c_H,c_Ks
WRITE(52,*) dummy_coeffient,dummy_coeffient,dummy_coeffient
WRITE(52,*) rms_J,rms_H,rms_Ks
WRITE(52,*) a_err_J,a_err_H,a_err_Ks
WRITE(52,*) b_err_J,b_err_H,b_err_Ks
WRITE(52,*) c_err_J,c_err_H,c_err_Ks
WRITE(52,*) 
WRITE(52,*) 'fit on observed data ... 2. iteration'
 CLOSE(52)
 
 






deallocate(model_J)
deallocate(model_Ks)
deallocate(delta_J)
deallocate(delta_Ks)

deallocate(model_H)
deallocate(delta_H)


allocate(model_J(1:hcount))
allocate(model_Ks(1:hcount))
allocate(delta_J(1:hcount))
allocate(delta_Ks(1:hcount)) 
 
allocate(model_H(1:hcount))
allocate(delta_H(1:hcount))
  
 delta_sum_J=0.D0
 delta_sum_Ks=0.D0
  delta_sum_H=0.D0
 
 DO i=1,hcount
 
model_J(i)=m_g(i)-(a_J*(m_g(i)-m_r(i))+b_J*(m_r(i)-m_i(i))+c_J)
model_Ks(i)=m_g(i)-(a_Ks*(m_g(i)-m_r(i))+b_Ks*(m_r(i)-m_i(i))+c_Ks)
model_H(i)=m_g(i)-(a_H*(m_g(i)-m_r(i))+b_H*(m_r(i)-m_i(i))+c_H)

! WRITE(*,*) model_J(i),m_J(i)
delta_J(i)=m_J(i)-model_J(i)
delta_Ks(i)=m_Ks(i)-model_Ks(i)
delta_sum_J=delta_sum_J+delta_J(i)**2
delta_sum_Ks=delta_sum_Ks+delta_Ks(i)**2

delta_H(i)=m_H(i)-model_H(i)
delta_sum_H=delta_sum_H+delta_H(i)**2

END DO

 delta_sum_J=SQRT(delta_sum_J/DBLE(hcount))
delta_sum_Ks=SQRT(delta_sum_Ks/DBLE(hcount))
delta_sum_H=SQRT(delta_sum_H/DBLE(hcount))













! WRITE(*,*) delta_sum_J,delta_sum_H,delta_sum_Ks


OPEN(50,file='edgeplane_J.txt')
DO i=1,hcount
WRITE(50,*) rightside_J(i),leftside_J(i)
END DO
 CLOSE(50)
OPEN(50,file='edgeplane_Ks.txt')
DO i=1,hcount
WRITE(50,*) rightside_Ks(i),leftside_Ks(i)
END DO
 CLOSE(50)
 OPEN(50,file='edgeplane_H.txt')
DO i=1,hcount
WRITE(50,*) rightside_H(i),leftside_H(i)
END DO
 CLOSE(50)
 
 
 
 
 
 
 
 
 
 
 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((leftside_J(i)-1.5D0)*25.D0)
help_y=NINT((rightside_J(i)-1.5D0)*25.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<75)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='edgeon_J.txt')
DO i=0,100
DO ii=0,75
help_x2=(DBLE(i)/25.D0)+1.5D0
help_y2=(DBLE(ii)/25.D0)+1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((leftside_Ks(i)-1.5D0)*25.D0)
help_y=NINT((rightside_Ks(i)-1.5D0)*25.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<75)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='edgeon_Ks.txt')
DO i=0,100
DO ii=0,75
help_x2=(DBLE(i)/25.D0)+1.5D0
help_y2=(DBLE(ii)/25.D0)+1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((leftside_H(i)-1.5D0)*25.D0)
help_y=NINT((rightside_H(i)-1.5D0)*25.D0)
IF ((help_x>0).AND.(help_x<100)) THEN
IF ((help_y>0).AND.(help_y<75)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='edgeon_H.txt')
DO i=0,100
DO ii=0,75
help_x2=(DBLE(i)/25.D0)+1.5D0
help_y2=(DBLE(ii)/25.D0)+1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


OPEN(50,file='delta_J_real.txt')
DO i=1,hcount
WRITE(50,*) m_J(i),delta_J(i)
END DO
 CLOSE(50)
OPEN(50,file='delta_Ks_real.txt')
DO i=1,hcount
WRITE(50,*) m_Ks(i),delta_Ks(i)
END DO
 CLOSE(50)
 OPEN(50,file='delta_H_real.txt')
DO i=1,hcount
WRITE(50,*) m_H(i),delta_H(i)
END DO
 CLOSE(50)
 
 
 
 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_J(i)-8.D0)*25.D0)
help_y=NINT((delta_J(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_J_real_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_Ks(i)-8.D0)*25.D0)
help_y=NINT((delta_Ks(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_Ks_real_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_H(i)-8.D0)*25.D0)
help_y=NINT((delta_H(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_H_real_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
 
 
 
 
 
 
 
  DO i=1,hcount
model_J(i)=m_g(i)-(a_J_bilir*(m_g(i)-m_r(i))+b_J_bilir*(m_r(i)-m_i(i))+c_J_bilir)
model_Ks(i)=m_g(i)-(a_Ks_bilir*(m_g(i)-m_r(i))+b_Ks_bilir*(m_r(i)-m_i(i))+c_Ks_bilir)
model_H(i)=m_g(i)-(a_H_bilir*(m_g(i)-m_r(i))+b_H_bilir*(m_r(i)-m_i(i))+c_H_bilir)

 delta_J(i)=m_J(i)-model_J(i)
delta_Ks(i)=m_Ks(i)-model_Ks(i)
delta_sum_J=delta_sum_J+delta_J(i)**2
delta_sum_Ks=delta_sum_Ks+delta_Ks(i)**2

delta_H(i)=m_H(i)-model_H(i)
delta_sum_H=delta_sum_H+delta_H(i)**2

END DO
 
  delta_sum_J=SQRT(delta_sum_J/DBLE(hcount))
delta_sum_Ks=SQRT(delta_sum_Ks/DBLE(hcount))
delta_sum_H=SQRT(delta_sum_H/DBLE(hcount))
 
OPEN(50,file='delta_J_bilir.txt')
DO i=1,hcount
WRITE(50,*) m_J(i),delta_J(i)
END DO
 CLOSE(50)
OPEN(50,file='delta_Ks_bilir.txt')
DO i=1,hcount
WRITE(50,*) m_Ks(i),delta_Ks(i)
END DO
 CLOSE(50)
 OPEN(50,file='delta_H_bilir.txt')
DO i=1,hcount
WRITE(50,*) m_H(i),delta_H(i)
END DO
 CLOSE(50)
 

 
 
 
 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_J(i)-8.D0)*25.D0)
help_y=NINT((delta_J(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_J_bilir_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_Ks(i)-8.D0)*25.D0)
help_y=NINT((delta_Ks(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_Ks_bilir_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_H(i)-8.D0)*25.D0)
help_y=NINT((delta_H(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_H_bilir_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
  
  DO i=1,hcount
model_J(i)=m_g(i)-(a_J_mill*(m_g(i)-m_r(i))+b_J_mill*(m_r(i)-m_i(i))+c_J_mill)
model_Ks(i)=m_g(i)-(a_Ks_mill*(m_g(i)-m_r(i))+b_Ks_mill*(m_r(i)-m_i(i))+c_Ks_mill)
model_H(i)=m_g(i)-(a_H_mill*(m_g(i)-m_r(i))+b_H_mill*(m_r(i)-m_i(i))+c_H_mill)

 delta_J(i)=m_J(i)-model_J(i)
delta_Ks(i)=m_Ks(i)-model_Ks(i)
delta_sum_J=delta_sum_J+delta_J(i)**2
delta_sum_Ks=delta_sum_Ks+delta_Ks(i)**2

delta_H(i)=m_H(i)-model_H(i)
delta_sum_H=delta_sum_H+delta_H(i)**2

END DO
 
 delta_sum_J=SQRT(delta_sum_J/DBLE(hcount))
delta_sum_Ks=SQRT(delta_sum_Ks/DBLE(hcount))
delta_sum_H=SQRT(delta_sum_H/DBLE(hcount))
 
OPEN(50,file='delta_J_mill.txt')
DO i=1,hcount
WRITE(50,*) m_J(i),delta_J(i)
END DO
 CLOSE(50)
OPEN(50,file='delta_Ks_mill.txt')
DO i=1,hcount
WRITE(50,*) m_Ks(i),delta_Ks(i)
END DO
 CLOSE(50)
 OPEN(50,file='delta_H_mill.txt')
DO i=1,hcount
WRITE(50,*) m_H(i),delta_H(i)
END DO
 CLOSE(50)
 
 
 
 
DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_J(i)-8.D0)*25.D0)
help_y=NINT((delta_J(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_J_mill_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_Ks(i)-8.D0)*25.D0)
help_y=NINT((delta_Ks(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_Ks_mill_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
 
 
 
 DO i=0,500
DO ii=0,400
mapmap(i,ii)=0
END DO
END DO


DO i=1,hcount

help_x=NINT((m_H(i)-8.D0)*25.D0)
help_y=NINT((delta_H(i)+3.D0)*25.D0)
IF ((help_x>0).AND.(help_x<125)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1
END IF
END IF


END DO


OPEN(61,file='delta_H_mill_map.txt')
DO i=0,125
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+8.D0
help_y2=(DBLE(ii)/25.D0)-3.D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)
 
WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


