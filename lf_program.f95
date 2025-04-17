PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
double precision :: PI,H0,q0,light,G,tophat,rho_crit
integer :: io_err,i,ii,iii,hint,counter,basicmode,old_fof_n,fof_n,n_dreieck,nmax
character(200), dimension(0:16) :: filename,directory,full_filename
character(200) :: appendix

double precision :: cV_mill,m_expected,Omega_m,Omega_l,h,part_exp,divisor,b_factor
double precision :: d_red,d_ang_dist,angular_sep,delta_z,hreal,av_z,av_l_a,av_l_z,l_next,hl_z,hl_a,l_old
double precision :: mag_max,mag_min,fit_max,fit_min

double precision, allocatable ::  ra(:,:),dec(:,:),z(:,:),mag_g(:,:),mag_r(:,:),volume_weight(:,:)
logical, allocatable :: active(:,:)
integer, dimension(0:16) :: n

double precision :: mag_limit,D_limit,t1,t2,t3,z_limit,V_limit,V_sum,binhelp
double precision :: sdss_cover,twomass_cover,mock_cover,Vmax

double precision, dimension(0:110,0:18) :: bin
double precision, dimension(0:30) :: binhelp2,bin_help_power,binhelp2power2
double precision, dimension(0:30,1:3) ::logmag,predicted_n,nbin,lnn

double precision, dimension(1:3,1:3) :: Amatrix,Bmatrix1
double precision, dimension(1:3,1:3) :: Bmatrix2,Bmatrix3
double precision, dimension(1:3) :: Yvector
double precision :: detA,detB1,detB2,detB3
double precision :: saturation_u,saturation_g,saturation_r,saturation_i,saturation_z
double precision :: Vminsat,Dmin

double precision, dimension(1:3)  :: a,b,c,rms,alpha,Mstar,Phistar

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme LUMINOSITY FUNCTION started'
WRITE(*,*) '============================================================'



! define constants
PI=ACOS(-1.D0)



Omega_m=0.25D0
Omega_l=0.75D0

sdss_cover=9274.D0/41253.D0
twomass_cover=0.91D0
mock_cover=1.D0/8.D0


q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
H0=73.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)
 cV_mill=(330.D6/h)**3
m_expected=cV_mill*rho_crit*Omega_m
part_exp=(2160.D0**3)*((330.D0/500.D0)**3)


z_limit=0.11D0
divisor=(SQRT(1.D0+2.D0*q0*z_limit)+1.D0+q0*z_limit)
D_limit=light/H0*z_limit*(1.D0+((z_limit*(1.D0-q0))/divisor))

OPEN(50,file='dist_limit.txt')
WRITE(50,*) 'dist_limit=',LOG10(D_limit*1.D6)
 CLOSE(50)

D_limit=D_limit*((1.D0+z_limit)**(-1.D0))
Vmax=4.D0/3.D0*PI*(D_limit**3)*mock_cover

 OPEN(52,file='saturation.txt')
READ(52,*) saturation_u
READ(52,*) saturation_g
READ(52,*) saturation_r
READ(52,*) saturation_i
READ(52,*) saturation_z
 CLOSE(52)
 



! Vmax=(350.D0-10.D0/h)**3
! mag_max=-30.D0
! mag_min=-15.D0

OPEN(50,file='inputfilelist_2MRS.txt')
DO i=0,16
READ(50,*) directory(i),filename(i)
full_filename(i)=TRIM(adjustl(directory(i)))//'/'//filename(i)
END DO
READ(50,*)  mag_max,mag_min
READ(50,*) appendix
READ(50,*) mag_limit
READ(50,*) fit_max,fit_min
appendix=TRIM(adjustl(appendix))
 CLOSE(50)
 

DO ii=0,16
! WRITE(*,*) full_filename(ii)
! get length of file
OPEN(50,file=full_filename(ii))

io_err=0
n(ii)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n(ii)=n(ii)+1
END DO
 CLOSE(50)
n(ii)=n(ii)-1

WRITE(*,*) n(ii),'galaxies will be used from '//directory(ii)

END DO



nmax=0
DO ii=0,16
IF (n(ii)>nmax) THEN
nmax=n(ii)
END IF
END DO


allocate(ra(1:nmax,0:16))
allocate(dec(1:nmax,0:16))
allocate(z(1:nmax,0:16))
allocate(mag_g(1:nmax,0:16))
allocate(mag_r(1:nmax,0:16))
allocate(volume_weight(1:nmax,0:16))
allocate(active(1:nmax,0:16))


DO ii=0,16
! read file
OPEN(50,file=full_filename(ii))
DO i=1,n(ii)
READ(50,*) ra(i,ii),dec(i,ii),z(i,ii),mag_g(i,ii),mag_r(i,ii)
END DO
CLOSE(50)

END DO


DO ii=0,16
V_sum=0.D0

DO i=1,n(ii)
active(i,ii)=.TRUE.
volume_weight(i,ii)=0.D0
IF (mag_r(i,ii)<mag_max) THEN
active(i,ii)=.FALSE.
END IF

IF (mag_r(i,ii)>mag_min) THEN
active(i,ii)=.FALSE.
END IF


IF (active(i,ii)) THEN
D_limit=-0.2D0*mag_r(i,ii)+((mag_limit+5.D0)/5.D0)
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

Vminsat=0.D0

IF (TRIM(appendix)=='SDSS') THEN
Dmin=(-mag_r(i,ii)+saturation_r+5.D0)/5.D0
Dmin=10.D0**Dmin
Dmin=Dmin/1.D6

t1=(light**2)*q0-(light**2)
t2=(light**4)*(q0**2)-2.D0*(light**4)*q0+(light**4)+2.D0*(light**3)*Dmin*H0*(q0**2)
t2=t2-4.D0*(light**3)*Dmin*H0*q0+2.D0*(light**3)*Dmin*H0
t2=SQRT(t2)
t3=light*Dmin*H0*q0

z_limit=1/(light**2)*(t1+t2+t3)
Dmin=Dmin*((1.D0+z_limit)**(-1.D0))!comoving diameter distance

Vminsat=(Dmin**3)

IF (Vminsat<0.D0) THEN
Vminsat=0.D0
END IF

END IF

V_limit=4.D0/3.D0*PI*((D_limit**3)-Vminsat)




IF (ii==0) THEN

IF (TRIM(appendix)=='SDSS') THEN
V_limit=V_limit*sdss_cover
ELSE
V_limit=V_limit*twomass_cover
END IF

END IF

IF ((ii>=1).AND.(ii<=8)) THEN
V_limit=V_limit*mock_cover
END IF




! WRITE(*,*) V_limit,D_limit,z_limit
volume_weight(i,ii)=1.D0/V_limit
! V_sum=V_sum+1.D0/V_limit



END IF
END DO

IF (ii>8) THEN
DO i=1,n(ii)

 IF (ra(i,ii)<0.D0) THEN
 active(i,ii)=.FALSE. 
 END IF
 IF (ra(i,ii)>90.D0) THEN
 active(i,ii)=.FALSE. 
 END IF
 IF (dec(i,ii)<0.D0) THEN
 active(i,ii)=.FALSE. 
 END IF
 IF (z(i,ii)>0.11D0) THEN
 active(i,ii)=.FALSE. 
 END IF


IF (active(i,ii)) THEN

volume_weight(i,ii)=1.D0/Vmax

END IF
END DO
END IF

END DO





DO ii=0,16

DO i=0,110
bin(i,ii)=0
END DO

DO i=1,n(ii)


DO iii=0,30
binhelp=(DBLE(iii)/2.D0)+mag_max
IF ((mag_r(i,ii)>(binhelp-(1.D0/4.D0))).AND.(mag_r(i,ii)<(binhelp+(1.D0/4.D0)))) THEN
IF (active(i,ii)) THEN
bin(iii,ii)=bin(iii,ii)+volume_weight(i,ii)  
END IF
END IF
END DO


END DO
END DO

OPEN(70,file='corrected_luminosity_function_'//TRIM(appendix)//'.txt')
DO i=0,30
binhelp=(DBLE(i)/2.D0)+mag_max
WRITE(70,*) binhelp,bin(i,0:16)
END DO
 CLOSE (70)



 DO ii=17,18
DO i=0,30
bin(i,ii)=0
END DO
 END DO
 
 
  DO ii=1,8
  DO i=0,30
 bin(i,17)= bin(i,17)+ bin(i,ii)
 END DO
 END DO
 
 
   DO ii=9,16
  DO i=0,30
 bin(i,18)= bin(i,18)+ bin(i,ii)
 END DO
 END DO

 DO ii=17,18
DO i=0,30
bin(i,ii)=bin(i,ii)/8.D0
END DO
 END DO
 
 
 OPEN(70,file='average_luminosity_function_'//TRIM(appendix)//'.txt')
DO i=0,30
binhelp=(DBLE(i)/2.D0)+mag_max
binhelp2(i)=(DBLE(i)/2.D0)+mag_max
bin_help_power(i)=10.D0**(-0.4D0*binhelp2(i))
WRITE(70,*) binhelp,bin(i,0),bin(i,17),bin(i,18)
END DO
 CLOSE (70)
 
 
 DO i=0,30
 
 logmag(i,1)=-1001.D0
 logmag(i,2)=-1001.D0
 logmag(i,3)=-1001.D0

 IF (binhelp2(i)>fit_max) THEN
  IF (binhelp2(i)<fit_min) THEN
 IF (bin(i,0)>0.D0) THEN
 logmag(i,1)=LOG10(bin(i,0))
 END IF
 
 IF (bin(i,17)>0.D0) THEN
 logmag(i,2)=LOG10(bin(i,17))
 END IF
 
 IF (bin(i,18)>0.D0) THEN
 logmag(i,3)=LOG10(bin(i,18))
 END IF
 END IF
  END IF
 END DO
 
 
 DO iii=1,3

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
DO i=0,30

IF (logmag(i,iii)>-1000.D0) THEN
Amatrix(1,1)=Amatrix(1,1)+bin_help_power(i)*bin_help_power(i)
Amatrix(2,2)=Amatrix(2,2)+binhelp2(i)*binhelp2(i)
Amatrix(3,3)=Amatrix(3,3)+1.D0
Amatrix(1,2)=Amatrix(1,2)+bin_help_power(i)*binhelp2(i)
Amatrix(2,1)=Amatrix(2,1)+bin_help_power(i)*binhelp2(i)
Amatrix(1,3)=Amatrix(1,3)+bin_help_power(i)
Amatrix(3,1)=Amatrix(3,1)+bin_help_power(i)
Amatrix(2,3)=Amatrix(2,3)+binhelp2(i)
Amatrix(3,2)=Amatrix(3,2)+binhelp2(i)
Yvector(1)=Yvector(1)+logmag(i,iii)*bin_help_power(i)
Yvector(2)=Yvector(2)+logmag(i,iii)*binhelp2(i)
Yvector(3)=Yvector(3)+logmag(i,iii)
END IF

END DO

!get determinat of Amatrix

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


! get solutions

 a(iii)=detB1/detA
 b(iii)=detB2/detA
 c(iii)=detB3/detA

 
alpha(iii)=(b(iii)/(-0.4D0))-1.D0
Mstar(iii)=LOG10(-LOG(10.D0)/10.D0*a(iii))/0.4D0
Phistar(iii)=(1.D0/(0.4D0*LOG(10.D0)))*(10.D0**(c(iii)-(alpha(iii)+1.D0)*0.4D0*Mstar(iii)))
 
 
 
 
 END DO
 
 
 DO iii=1,3
 hint=0
 rms(iii)=0.D0
 DO i=0,30
 predicted_n(i,iii)=0.D0
 IF (binhelp2(i)>fit_max) THEN
 IF (binhelp2(i)<fit_min) THEN
 predicted_n(i,iii)=0.4D0*LOG(10.D0)*Phistar(iii)*((10.D0**(-0.4D0*(binhelp2(i)-Mstar(iii))))**(alpha(iii)+1.D0))
 predicted_n(i,iii)=predicted_n(i,iii)*EXP(-(10.D0**(-0.4D0*(binhelp2(i)-Mstar(iii))))) 
 hint=hint+1
 rms(iii)=rms(iii)+(predicted_n(i,iii)-(10.D0**logmag(i,iii)))**2
 END IF
 END IF
 END DO
 rms(iii)=SQRT(rms(iii)/DBLE(hint))
 END DO
 
 
  OPEN(70,file='fit_luminosity_function_'//TRIM(appendix)//'.txt')
WRITE(70,*) a(1:3)
WRITE(70,*) b(1:3)
WRITE(70,*) c(1:3)
WRITE(70,*) '------------------------------'
WRITE(70,*) alpha(1:3)
WRITE(70,*) Mstar(1:3)
WRITE(70,*) Phistar(1:3)
WRITE(70,*) '------------------------------'
WRITE(70,*) rms(1:3)
 CLOSE (70)
 
 

  
  OPEN(70,file='plot_luminosity_function_'//TRIM(appendix)//'.txt')
WRITE(70,*) "f1(x)= ",a(1),"* (10.**(-0.4*x)) + ",b(1),"* x +",c(1)
WRITE(70,*) "f2(x)= ",a(2),"* (10.**(-0.4*x)) + ",b(2),"* x +",c(2)
WRITE(70,*) "f3(x)= ",a(3),"* (10.**(-0.4*x)) + ",b(3),"* x +",c(3)

 CLOSE (70)
 
 
 
 
  OPEN(70,file='log_luminosity_function_'//TRIM(appendix)//'.txt')
DO i=0,30
WRITE(70,*) binhelp2(i),logmag(i,1:3)
END DO
 CLOSE (70)
 
 
 
 
 
 
 

 
 
 
 

 
 
 

 DO i=0,30
 binhelp2power2(i)=binhelp2(i)**2
nbin(i,1)=bin(i,0)
nbin(i,2)=bin(i,17)
nbin(i,3)=bin(i,18)
END DO

 DO iii=1,3

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
DO i=0,30


Amatrix(1,1)=Amatrix(1,1)+binhelp2power2(i)*binhelp2power2(i)
Amatrix(2,2)=Amatrix(2,2)+binhelp2(i)*binhelp2(i)
Amatrix(3,3)=Amatrix(3,3)+1.D0
Amatrix(1,2)=Amatrix(1,2)+binhelp2power2(i)*binhelp2(i)
Amatrix(2,1)=Amatrix(2,1)+binhelp2power2(i)*binhelp2(i)
Amatrix(1,3)=Amatrix(1,3)+binhelp2power2(i)
Amatrix(3,1)=Amatrix(3,1)+binhelp2power2(i)
Amatrix(2,3)=Amatrix(2,3)+binhelp2(i)
Amatrix(3,2)=Amatrix(3,2)+binhelp2(i)
Yvector(1)=Yvector(1)+nbin(i,iii)*binhelp2power2(i)
Yvector(2)=Yvector(2)+nbin(i,iii)*binhelp2(i)
Yvector(3)=Yvector(3)+nbin(i,iii)

END DO

!get determinat of Amatrix

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


! get solutions

 a(iii)=detB1/detA
 b(iii)=detB2/detA
 c(iii)=detB3/detA


 
 
 END DO
 
 
 DO iii=1,3
 hint=0
 rms(iii)=0.D0
 DO i=0,30
 predicted_n(i,iii)=0.D0
 IF (binhelp2(i)>fit_max) THEN
 IF (binhelp2(i)<fit_min) THEN
 predicted_n(i,iii)=binhelp2power2(i)*a(iii)+binhelp2(i)*b(iii)+c(iii)
 hint=hint+1
 rms(iii)=rms(iii)+(predicted_n(i,iii)-(nbin(i,iii)))**2
 END IF
 END IF
 END DO
 rms(iii)=SQRT(rms(iii)/DBLE(hint))
 END DO
 
 
  OPEN(70,file='fitquad_luminosity_function_'//TRIM(appendix)//'.txt')
WRITE(70,*) a(1:3)
WRITE(70,*) b(1:3)
WRITE(70,*) c(1:3)
WRITE(70,*) '------------------------------'
WRITE(70,*) rms(1:3)
 CLOSE (70)
 
   OPEN(70,file='fitquad_plot_'//TRIM(appendix)//'.txt')
WRITE(70,*) "f1(x) = ",a(1),"* x*x + ",b(1),"* x + ",c(1)
WRITE(70,*) "f2(x) = ",a(2),"* x*x + ",b(2),"* x + ",c(2)
WRITE(70,*) "f3(x) = ",a(3),"* x*x + ",b(3),"* x + ",c(3)
 CLOSE (70)
 

 
deallocate(ra)
deallocate(dec)
deallocate(z)
deallocate(mag_g)
deallocate(mag_r)
deallocate(volume_weight)
deallocate(active)



OPEN(50,file='inputfilelist_2MRS.txt')
DO i=0,16
READ(50,*) directory(i),filename(i)
full_filename(i)=TRIM(adjustl(directory(i)))//'/'//filename(i)
END DO
READ(50,*)  mag_max,mag_min
READ(50,*) appendix
READ(50,*) mag_limit
READ(50,*) fit_max,fit_min
appendix=TRIM(adjustl(appendix))
 CLOSE(50)
 

DO ii=0,16
! WRITE(*,*) full_filename(ii)
! get length of file
OPEN(50,file=full_filename(ii))

io_err=0
n(ii)=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n(ii)=n(ii)+1
END DO
 CLOSE(50)
n(ii)=n(ii)-1

WRITE(*,*) n(ii),'galaxies will be used from '//directory(ii)

END DO



nmax=0
DO ii=0,16
IF (n(ii)>nmax) THEN
nmax=n(ii)
END IF
END DO


allocate(ra(1:nmax,0:16))
allocate(dec(1:nmax,0:16))
allocate(z(1:nmax,0:16))
allocate(mag_g(1:nmax,0:16))
allocate(mag_r(1:nmax,0:16))
allocate(volume_weight(1:nmax,0:16))
allocate(active(1:nmax,0:16))


DO ii=0,16
! read file
OPEN(50,file=full_filename(ii))
DO i=1,n(ii)
READ(50,*) ra(i,ii),dec(i,ii),z(i,ii),mag_g(i,ii),mag_r(i,ii)
END DO
CLOSE(50)

END DO


DO ii=0,16
V_sum=0.D0

DO i=1,n(ii)
active(i,ii)=.TRUE.
volume_weight(i,ii)=0.D0
IF (mag_r(i,ii)<mag_max) THEN
active(i,ii)=.FALSE.
END IF

IF (mag_r(i,ii)>mag_min) THEN
active(i,ii)=.FALSE.
END IF


IF (active(i,ii)) THEN
D_limit=-0.2D0*mag_r(i,ii)+((mag_limit+5.D0)/5.D0)
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




IF (ii==0) THEN

IF (TRIM(appendix)=='SDSS') THEN
V_limit=V_limit*sdss_cover
ELSE
V_limit=V_limit*twomass_cover
END IF

END IF

IF ((ii>=1).AND.(ii<=8)) THEN
V_limit=V_limit*mock_cover
END IF




! WRITE(*,*) V_limit,D_limit,z_limit
volume_weight(i,ii)=1.D0/V_limit
! V_sum=V_sum+1.D0/V_limit



END IF
END DO

IF (ii>8) THEN
DO i=1,n(ii)

 IF (ra(i,ii)<0.D0) THEN
 active(i,ii)=.FALSE. 
 END IF
 IF (ra(i,ii)>90.D0) THEN
 active(i,ii)=.FALSE. 
 END IF
 IF (dec(i,ii)<0.D0) THEN
 active(i,ii)=.FALSE. 
 END IF
 IF (z(i,ii)>0.11D0) THEN
 active(i,ii)=.FALSE. 
 END IF


IF (active(i,ii)) THEN

volume_weight(i,ii)=1.D0/Vmax

END IF
END DO
END IF

END DO





DO ii=0,16

DO i=0,110
bin(i,ii)=0
END DO

DO i=1,n(ii)


DO iii=0,30
binhelp=(DBLE(iii)/2.D0)+mag_max
IF ((mag_r(i,ii)>(binhelp-(1.D0/4.D0))).AND.(mag_r(i,ii)<(binhelp+(1.D0/4.D0)))) THEN
IF (active(i,ii)) THEN
bin(iii,ii)=bin(iii,ii)+volume_weight(i,ii)  
END IF
END IF
END DO


END DO
END DO

OPEN(70,file='corrected_luminosity_function_'//TRIM(appendix)//'.txt')
DO i=0,30
binhelp=(DBLE(i)/2.D0)+mag_max
WRITE(70,*) binhelp,bin(i,0:16)
END DO
 CLOSE (70)



 DO ii=17,18
DO i=0,30
bin(i,ii)=0
END DO
 END DO
 
 
  DO ii=1,8
  DO i=0,30
 bin(i,17)= bin(i,17)+ bin(i,ii)
 END DO
 END DO
 
 
   DO ii=9,16
  DO i=0,30
 bin(i,18)= bin(i,18)+ bin(i,ii)
 END DO
 END DO

 DO ii=17,18
DO i=0,30
bin(i,ii)=bin(i,ii)/8.D0
END DO
 END DO
 
 
 OPEN(70,file='average_luminosity_function_'//TRIM(appendix)//'.txt')
DO i=0,30
binhelp=(DBLE(i)/2.D0)+mag_max
binhelp2(i)=(DBLE(i)/2.D0)+mag_max
bin_help_power(i)=10.D0**(-0.4D0*binhelp2(i))
WRITE(70,*) binhelp,bin(i,0),bin(i,17),bin(i,18)
END DO
 CLOSE (70)
 
 
 DO i=0,30
 
 logmag(i,1)=-1001.D0
 logmag(i,2)=-1001.D0
 logmag(i,3)=-1001.D0

 IF (binhelp2(i)>fit_max) THEN
  IF (binhelp2(i)<fit_min) THEN
 IF (bin(i,0)>0.D0) THEN
 logmag(i,1)=LOG10(bin(i,0))
 END IF
 
 IF (bin(i,17)>0.D0) THEN
 logmag(i,2)=LOG10(bin(i,17))
 END IF
 
 IF (bin(i,18)>0.D0) THEN
 logmag(i,3)=LOG10(bin(i,18))
 END IF
 END IF
  END IF
 END DO
 
 
 DO iii=1,3

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
DO i=0,30

IF (logmag(i,iii)>-1000.D0) THEN
Amatrix(1,1)=Amatrix(1,1)+bin_help_power(i)*bin_help_power(i)
Amatrix(2,2)=Amatrix(2,2)+binhelp2(i)*binhelp2(i)
Amatrix(3,3)=Amatrix(3,3)+1.D0
Amatrix(1,2)=Amatrix(1,2)+bin_help_power(i)*binhelp2(i)
Amatrix(2,1)=Amatrix(2,1)+bin_help_power(i)*binhelp2(i)
Amatrix(1,3)=Amatrix(1,3)+bin_help_power(i)
Amatrix(3,1)=Amatrix(3,1)+bin_help_power(i)
Amatrix(2,3)=Amatrix(2,3)+binhelp2(i)
Amatrix(3,2)=Amatrix(3,2)+binhelp2(i)
Yvector(1)=Yvector(1)+logmag(i,iii)*bin_help_power(i)
Yvector(2)=Yvector(2)+logmag(i,iii)*binhelp2(i)
Yvector(3)=Yvector(3)+logmag(i,iii)
END IF

END DO

!get determinat of Amatrix

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


! get solutions

 a(iii)=detB1/detA
 b(iii)=detB2/detA
 c(iii)=detB3/detA

 
alpha(iii)=(b(iii)/(-0.4D0))-1.D0
Mstar(iii)=LOG10(-LOG(10.D0)/10.D0*a(iii))/0.4D0
Phistar(iii)=(1.D0/(0.4D0*LOG(10.D0)))*(10.D0**(c(iii)-(alpha(iii)+1.D0)*0.4D0*Mstar(iii)))
 
 
 
 
 END DO
 
 
 DO iii=1,3
 hint=0
 rms(iii)=0.D0
 DO i=0,30
 predicted_n(i,iii)=0.D0
 IF (binhelp2(i)>fit_max) THEN
 IF (binhelp2(i)<fit_min) THEN
 predicted_n(i,iii)=0.4D0*LOG(10.D0)*Phistar(iii)*((10.D0**(-0.4D0*(binhelp2(i)-Mstar(iii))))**(alpha(iii)+1.D0))
 predicted_n(i,iii)=predicted_n(i,iii)*EXP(-(10.D0**(-0.4D0*(binhelp2(i)-Mstar(iii))))) 
 hint=hint+1
 rms(iii)=rms(iii)+(predicted_n(i,iii)-(10.D0**logmag(i,iii)))**2
 END IF
 END IF
 END DO
 rms(iii)=SQRT(rms(iii)/DBLE(hint))
 END DO
 
 
  OPEN(70,file='fit_luminosity_function_'//TRIM(appendix)//'.txt')
WRITE(70,*) a(1:3)
WRITE(70,*) b(1:3)
WRITE(70,*) c(1:3)
WRITE(70,*) '------------------------------'
WRITE(70,*) alpha(1:3)
WRITE(70,*) Mstar(1:3)
WRITE(70,*) Phistar(1:3)
WRITE(70,*) '------------------------------'
WRITE(70,*) rms(1:3)
 CLOSE (70)
 
 

  
  OPEN(70,file='plot_luminosity_function_'//TRIM(appendix)//'.txt')
WRITE(70,*) "f1(x)= ",a(1),"* (10.**(-0.4*x)) + ",b(1),"* x +",c(1)
WRITE(70,*) "f2(x)= ",a(2),"* (10.**(-0.4*x)) + ",b(2),"* x +",c(2)
WRITE(70,*) "f3(x)= ",a(3),"* (10.**(-0.4*x)) + ",b(3),"* x +",c(3)

 CLOSE (70)
 
 
 
 
  OPEN(70,file='log_luminosity_function_'//TRIM(appendix)//'.txt')
DO i=0,30
WRITE(70,*) binhelp2(i),logmag(i,1:3)
END DO
 CLOSE (70)
 
 
 
 
 
 
 

 
 
 
 

 
 
 

 DO i=0,30
 binhelp2power2(i)=binhelp2(i)**2
nbin(i,1)=bin(i,0)
nbin(i,2)=bin(i,17)
nbin(i,3)=bin(i,18)
END DO

 DO iii=1,3

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
DO i=0,30


Amatrix(1,1)=Amatrix(1,1)+binhelp2power2(i)*binhelp2power2(i)
Amatrix(2,2)=Amatrix(2,2)+binhelp2(i)*binhelp2(i)
Amatrix(3,3)=Amatrix(3,3)+1.D0
Amatrix(1,2)=Amatrix(1,2)+binhelp2power2(i)*binhelp2(i)
Amatrix(2,1)=Amatrix(2,1)+binhelp2power2(i)*binhelp2(i)
Amatrix(1,3)=Amatrix(1,3)+binhelp2power2(i)
Amatrix(3,1)=Amatrix(3,1)+binhelp2power2(i)
Amatrix(2,3)=Amatrix(2,3)+binhelp2(i)
Amatrix(3,2)=Amatrix(3,2)+binhelp2(i)
Yvector(1)=Yvector(1)+nbin(i,iii)*binhelp2power2(i)
Yvector(2)=Yvector(2)+nbin(i,iii)*binhelp2(i)
Yvector(3)=Yvector(3)+nbin(i,iii)

END DO

!get determinat of Amatrix

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


! get solutions

 a(iii)=detB1/detA
 b(iii)=detB2/detA
 c(iii)=detB3/detA


 
 
 END DO
 
 
 DO iii=1,3
 hint=0
 rms(iii)=0.D0
 DO i=0,30
 predicted_n(i,iii)=0.D0
 IF (binhelp2(i)>fit_max) THEN
 IF (binhelp2(i)<fit_min) THEN
 predicted_n(i,iii)=binhelp2power2(i)*a(iii)+binhelp2(i)*b(iii)+c(iii)
 hint=hint+1
 rms(iii)=rms(iii)+(predicted_n(i,iii)-(nbin(i,iii)))**2
 END IF
 END IF
 END DO
 rms(iii)=SQRT(rms(iii)/DBLE(hint))
 END DO
 
 
  OPEN(70,file='fitquad_luminosity_function_'//TRIM(appendix)//'.txt')
WRITE(70,*) a(1:3)
WRITE(70,*) b(1:3)
WRITE(70,*) c(1:3)
WRITE(70,*) '------------------------------'
WRITE(70,*) rms(1:3)
 CLOSE (70)
 
   OPEN(70,file='fitquad_plot_'//TRIM(appendix)//'.txt')
WRITE(70,*) "f1(x) = ",a(1),"* x*x + ",b(1),"* x + ",c(1)
WRITE(70,*) "f2(x) = ",a(2),"* x*x + ",b(2),"* x + ",c(2)
WRITE(70,*) "f3(x) = ",a(3),"* x*x + ",b(3),"* x + ",c(3)
 CLOSE (70)
 

 
  
 
 
 
 
 
 
 

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


