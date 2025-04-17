PROGRAM recalibration
! declaration of variables
IMPLICIT NONE
! basic constants

double precision :: rms_sdss_s,rms_2mrs_s,rms_sdss_mh,rms_2mrs_mh,rms_sdss_ml,rms_2mrs_ml
double precision, dimension(1:7) :: solution_sdss_s,errorbar_sdss_s,solution_2mrs_s,errorbar_2mrs_s
double precision, dimension(1:9) :: solution_sdss_ml,errorbar_sdss_ml,solution_2mrs_ml,errorbar_2mrs_ml
double precision, dimension(1:10) :: solution_sdss_mh,errorbar_sdss_mh,solution_2mrs_mh,errorbar_2mrs_mh

double precision :: mass_fit_sum_sdss_s,mass_fit_sum_2mrs_s,mass_orig_sum_sdss_s,mass_orig_sum_2mrs_s
double precision :: mass_fit_sum_sdss_ml,mass_fit_sum_2mrs_ml,mass_orig_sum_sdss_ml,mass_orig_sum_2mrs_ml
double precision :: mass_fit_sum_sdss_mh,mass_fit_sum_2mrs_mh,mass_orig_sum_sdss_mh,mass_orig_sum_2mrs_mh

double precision :: mass_fit_sum_sdss_tot,mass_fit_sum_2mrs_tot
double precision :: mass_orig_sum_sdss_tot,mass_orig_sum_2mrs_tot
double precision :: mag_sol_r,mag_sol_g,mag_sol_Ks,mag_sol_J,dummy_var
double precision :: Omega_m,Omega_l,H0,PI,q0,light,tophat,h,G,rho_crit,z_max,divisor,d_com,area_cover,vol,m_expected
double precision :: sigma_lumtot_sdss_s,sigma_lumtot_sdss_ml,sigma_lumtot_sdss_mh
double precision :: sigma_lumtot_2mrs_s,sigma_lumtot_2mrs_ml,sigma_lumtot_2mrs_mh


integer :: io_err,i,ii,help_count,n_ml_sdss,n_mh_sdss,n_s_sdss,n_ml_2mrs,n_mh_2mrs,n_s_2mrs

double precision, allocatable ::  h_mass_mh_sdss(:),lumsum_cor_mh_sdss(:),lumdist_mh_sdss(:),m_dyn_mh_sdss(:)
double precision, allocatable :: sigma_mh_sdss(:)
double precision, allocatable :: radius_mh_sdss(:),n_visible_mh_sdss(:),lumsum_mh_sdss(:),lumsum_all_mh_sdss(:)

double precision, allocatable ::  h_mass_ml_sdss(:),lumsum_cor_ml_sdss(:),lumdist_ml_sdss(:),m_dyn_ml_sdss(:)
double precision, allocatable :: radius_ml_sdss(:),lumsum_ml_sdss(:),lumsum_all_ml_sdss(:),sigma_ml_sdss(:)

double precision, allocatable ::  h_mass_s_sdss(:),lumsum_cor_s_sdss(:),lumdist_s_sdss(:),lumsum_s_sdss(:)
double precision, allocatable :: lumsum_all_s_sdss(:)

double precision, allocatable ::  h_mass_mh_2mrs(:),lumsum_cor_mh_2mrs(:),lumdist_mh_2mrs(:),m_dyn_mh_2mrs(:)
double precision, allocatable :: sigma_mh_2mrs(:)
double precision, allocatable :: radius_mh_2mrs(:),n_visible_mh_2mrs(:),lumsum_mh_2mrs(:),lumsum_all_mh_2mrs(:)

double precision, allocatable ::  h_mass_ml_2mrs(:),lumsum_cor_ml_2mrs(:),lumdist_ml_2mrs(:),m_dyn_ml_2mrs(:)
double precision, allocatable :: radius_ml_2mrs(:),lumsum_ml_2mrs(:),lumsum_all_ml_2mrs(:),sigma_ml_2mrs(:)

double precision, allocatable ::  h_mass_s_2mrs(:),lumsum_cor_s_2mrs(:),lumdist_s_2mrs(:),lumsum_s_2mrs(:)
double precision, allocatable :: lumsum_all_s_2mrs(:)

double precision, allocatable :: otherside_sdss_s(:),otherside_sdss_ml(:),otherside_sdss_mh(:)
double precision, allocatable :: deviation_sdss_s(:),deviation_sdss_ml(:),deviation_sdss_mh(:)

double precision, allocatable :: otherside_2mrs_s(:),otherside_2mrs_ml(:),otherside_2mrs_mh(:)
double precision, allocatable :: deviation_2mrs_s(:),deviation_2mrs_ml(:),deviation_2mrs_mh(:)

logical, allocatable :: act_2mrs_s(:),act_2mrs_ml(:),act_2mrs_mh(:),act_sdss_s(:),act_sdss_ml(:),act_sdss_mh(:)

double precision, dimension(0:1000,0:1000) :: mapmap
double precision :: help_x2,help_y2
integer :: help_x,help_y

 Omega_m=0.25D0
 Omega_l=0.75D0
 H0=73.D0
 
 
! define constants
PI=ACOS(-1.D0)
 
q0=Omega_m/2.D0-Omega_l
light=3.D5
tophat=200.D0
h=H0/100.D0
G=4.302D-3 !in pc/Msol * (km/s)**2
rho_crit=3.D0*((1.D-6*H0)**2)/(8.D0*PI*G)

z_max=0.11D0
divisor=(SQRT(1.D0+2.D0*q0*z_max)+1.D0+q0*z_max)
d_com=light/H0*z_max*(1.D0+((z_max*(1.D0-q0))/divisor))
d_com=d_com/(1.D0+z_max)*1.D6

area_cover=(1.D0/8.D0)*8.D0
vol=4.D0*PI/3.D0*(d_com**3)/area_cover

m_expected=vol*rho_crit*Omega_m




! get length of file
OPEN(50,file='mass/log_dep_multi_high_2MRS_combi.txt')
io_err=0
n_mh_2mrs=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mh_2mrs=n_mh_2mrs+1
END DO
 CLOSE(50)
n_mh_2mrs=n_mh_2mrs-1

! get length of file
OPEN(50,file='mass/log_dep_multi_low_2MRS_combi.txt')
io_err=0
n_ml_2mrs=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_ml_2mrs=n_ml_2mrs+1
END DO
 CLOSE(50)
n_ml_2mrs=n_ml_2mrs-1

! get length of file
OPEN(50,file='mass/log_dep_single_2MRS_combi.txt')
io_err=0
n_s_2mrs=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_s_2mrs=n_s_2mrs+1
END DO
 CLOSE(50)
n_s_2mrs=n_s_2mrs-1



! get length of file
OPEN(50,file='mass/log_dep_multi_high_SDSS_combi.txt')
io_err=0
n_mh_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_mh_sdss=n_mh_sdss+1
END DO
 CLOSE(50)
n_mh_sdss=n_mh_sdss-1

! get length of file
OPEN(50,file='mass/log_dep_multi_low_SDSS_combi.txt')
io_err=0
n_ml_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_ml_sdss=n_ml_sdss+1
END DO
 CLOSE(50)
n_ml_sdss=n_ml_sdss-1

! get length of file
OPEN(50,file='mass/log_dep_single_SDSS_combi.txt')
io_err=0
n_s_sdss=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n_s_sdss=n_s_sdss+1
END DO
 CLOSE(50)
n_s_sdss=n_s_sdss-1




allocate(h_mass_mh_2mrs(1:n_mh_2mrs))
allocate(lumsum_cor_mh_2mrs(1:n_mh_2mrs))
allocate(lumdist_mh_2mrs(1:n_mh_2mrs))
allocate(m_dyn_mh_2mrs(1:n_mh_2mrs))
allocate(sigma_mh_2mrs(1:n_mh_2mrs))
allocate(radius_mh_2mrs(1:n_mh_2mrs))
allocate(n_visible_mh_2mrs(1:n_mh_2mrs))
allocate(lumsum_mh_2mrs(1:n_mh_2mrs))
allocate(lumsum_all_mh_2mrs(1:n_mh_2mrs))

allocate(h_mass_ml_2mrs(1:n_ml_2mrs))
allocate(lumsum_cor_ml_2mrs(1:n_ml_2mrs))
allocate(lumdist_ml_2mrs(1:n_ml_2mrs))
allocate(m_dyn_ml_2mrs(1:n_ml_2mrs))
allocate(radius_ml_2mrs(1:n_ml_2mrs))
allocate(lumsum_ml_2mrs(1:n_ml_2mrs))
allocate(lumsum_all_ml_2mrs(1:n_ml_2mrs))
allocate(sigma_ml_2mrs(1:n_ml_2mrs))

allocate(h_mass_s_2mrs(1:n_s_2mrs))
allocate(lumsum_cor_s_2mrs(1:n_s_2mrs))
allocate(lumdist_s_2mrs(1:n_s_2mrs))
allocate(lumsum_s_2mrs(1:n_s_2mrs))
allocate(lumsum_all_s_2mrs(1:n_s_2mrs))


allocate(h_mass_mh_sdss(1:n_mh_sdss))
allocate(lumsum_cor_mh_sdss(1:n_mh_sdss))
allocate(lumdist_mh_sdss(1:n_mh_sdss))
allocate(m_dyn_mh_sdss(1:n_mh_sdss))
allocate(sigma_mh_sdss(1:n_mh_sdss))
allocate(radius_mh_sdss(1:n_mh_sdss))
allocate(n_visible_mh_sdss(1:n_mh_sdss))
allocate(lumsum_mh_sdss(1:n_mh_sdss))
allocate(lumsum_all_mh_sdss(1:n_mh_sdss))

allocate(h_mass_ml_sdss(1:n_ml_sdss))
allocate(lumsum_cor_ml_sdss(1:n_ml_sdss))
allocate(lumdist_ml_sdss(1:n_ml_sdss))
allocate(m_dyn_ml_sdss(1:n_ml_sdss))
allocate(radius_ml_sdss(1:n_ml_sdss))
allocate(lumsum_ml_sdss(1:n_ml_sdss))
allocate(lumsum_all_ml_sdss(1:n_ml_sdss))
allocate(sigma_ml_sdss(1:n_ml_sdss))

allocate(h_mass_s_sdss(1:n_s_sdss))
allocate(lumsum_cor_s_sdss(1:n_s_sdss))
allocate(lumdist_s_sdss(1:n_s_sdss))
allocate(lumsum_s_sdss(1:n_s_sdss))
allocate(lumsum_all_s_sdss(1:n_s_sdss))


allocate(otherside_sdss_s(1:n_s_sdss))
allocate(otherside_sdss_ml(1:n_ml_sdss))
allocate(otherside_sdss_mh(1:n_mh_sdss))
allocate(deviation_sdss_s(1:n_s_sdss))
allocate(deviation_sdss_ml(1:n_ml_sdss))
allocate(deviation_sdss_mh(1:n_mh_sdss))

allocate(otherside_2mrs_s(1:n_s_2mrs))
allocate(otherside_2mrs_ml(1:n_ml_2mrs))
allocate(otherside_2mrs_mh(1:n_mh_2mrs))
allocate(deviation_2mrs_s(1:n_s_2mrs))
allocate(deviation_2mrs_ml(1:n_ml_2mrs))
allocate(deviation_2mrs_mh(1:n_mh_2mrs))


allocate(act_2mrs_s(1:n_s_2mrs))
allocate(act_2mrs_ml(1:n_ml_2mrs))
allocate(act_2mrs_mh(1:n_mh_2mrs))
allocate(act_sdss_s(1:n_s_sdss))
allocate(act_sdss_ml(1:n_ml_sdss))
allocate(act_sdss_mh(1:n_mh_sdss))



! get length of file
OPEN(50,file='mass/log_dep_multi_high_2MRS_combi.txt')
DO i=1,n_mh_2mrs
READ(50,*) h_mass_mh_2mrs(i),lumsum_cor_mh_2mrs(i),lumdist_mh_2mrs(i),m_dyn_mh_2mrs(i),sigma_mh_2mrs(i),&
radius_mh_2mrs(i),n_visible_mh_2mrs(i),lumsum_mh_2mrs(i),lumsum_all_mh_2mrs(i)
END DO
 CLOSE(50)
 
! get length of file
OPEN(50,file='mass/log_dep_multi_low_2MRS_combi.txt')
DO i=1,n_ml_2mrs
READ(50,*) h_mass_ml_2mrs(i),lumsum_cor_ml_2mrs(i),lumdist_ml_2mrs(i),m_dyn_ml_2mrs(i),sigma_ml_2mrs(i),&
radius_ml_2mrs(i),lumsum_ml_2mrs(i),lumsum_all_ml_2mrs(i)
END DO
 CLOSE(50)
 
! get length of file
OPEN(50,file='mass/log_dep_single_2MRS_combi.txt')
DO i=1,n_s_2mrs
READ(50,*) h_mass_s_2mrs(i),lumsum_cor_s_2mrs(i),lumdist_s_2mrs(i),lumsum_s_2mrs(i),lumsum_all_s_2mrs(i)
END DO
 CLOSE(50)


 
! get length of file
OPEN(50,file='mass/log_dep_multi_high_SDSS_combi.txt')
DO i=1,n_mh_sdss
READ(50,*) h_mass_mh_sdss(i),lumsum_cor_mh_sdss(i),lumdist_mh_sdss(i),m_dyn_mh_sdss(i),sigma_mh_sdss(i),&
radius_mh_sdss(i),n_visible_mh_sdss(i),lumsum_mh_sdss(i),lumsum_all_mh_sdss(i)
END DO
 CLOSE(50)

! get length of file
OPEN(50,file='mass/log_dep_multi_low_SDSS_combi.txt')
DO i=1,n_ml_sdss
READ(50,*) h_mass_ml_sdss(i),lumsum_cor_ml_sdss(i),lumdist_ml_sdss(i),m_dyn_ml_sdss(i),sigma_ml_sdss(i),&
radius_ml_sdss(i),lumsum_ml_sdss(i),lumsum_all_ml_sdss(i)
END DO
 CLOSE(50)

! get length of file
OPEN(50,file='mass/log_dep_single_SDSS_combi.txt')
DO i=1,n_s_sdss
READ(50,*) h_mass_s_sdss(i),lumsum_cor_s_sdss(i),lumdist_s_sdss(i),lumsum_s_sdss(i),lumsum_all_s_sdss(i)
END DO
 CLOSE(50)


 
 sigma_lumtot_sdss_s=0.D0
 sigma_lumtot_sdss_ml=0.D0
 sigma_lumtot_sdss_mh=0.D0
 sigma_lumtot_2mrs_s=0.D0
 sigma_lumtot_2mrs_ml=0.D0
 sigma_lumtot_2mrs_mh=0.D0
 
DO i=1,n_mh_2mrs
sigma_lumtot_2mrs_mh=sigma_lumtot_2mrs_mh+(lumsum_cor_mh_2mrs(i)-lumsum_all_mh_2mrs(i))**2
END DO

DO i=1,n_ml_2mrs
sigma_lumtot_2mrs_ml=sigma_lumtot_2mrs_ml+(lumsum_cor_ml_2mrs(i)-lumsum_all_ml_2mrs(i))**2
END DO

DO i=1,n_s_2mrs
sigma_lumtot_2mrs_s=sigma_lumtot_2mrs_s+(lumsum_cor_s_2mrs(i)-lumsum_all_s_2mrs(i))**2
END DO


DO i=1,n_mh_sdss
sigma_lumtot_sdss_mh=sigma_lumtot_sdss_mh+(lumsum_cor_mh_sdss(i)-lumsum_all_mh_sdss(i))**2
END DO

DO i=1,n_ml_sdss
sigma_lumtot_sdss_ml=sigma_lumtot_sdss_ml+(lumsum_cor_ml_sdss(i)-lumsum_all_ml_sdss(i))**2
END DO

DO i=1,n_s_sdss
sigma_lumtot_sdss_s=sigma_lumtot_sdss_s+(lumsum_cor_s_sdss(i)-lumsum_all_s_sdss(i))**2
END DO


sigma_lumtot_2mrs_mh=SQRT(sigma_lumtot_2mrs_mh/n_mh_2mrs)
sigma_lumtot_2mrs_ml=SQRT(sigma_lumtot_2mrs_ml/n_ml_2mrs)
sigma_lumtot_2mrs_s=SQRT(sigma_lumtot_2mrs_s/n_s_2mrs)

sigma_lumtot_sdss_mh=SQRT(sigma_lumtot_sdss_mh/n_mh_sdss)
sigma_lumtot_sdss_ml=SQRT(sigma_lumtot_sdss_ml/n_ml_sdss)
sigma_lumtot_sdss_s=SQRT(sigma_lumtot_sdss_s/n_s_sdss)


  
    
  OPEN(61,file='mass/fit_results.txt')
READ(61,*) rms_sdss_s
DO i=1,7
READ(61,*) solution_sdss_s(i),errorbar_sdss_s(i)
END DO
READ(61,*) 
READ(61,*) rms_2mrs_s
DO i=1,7
READ(61,*) solution_2mrs_s(i),errorbar_2mrs_s(i)
END DO
READ(61,*) 
READ(61,*) rms_sdss_ml
DO i=1,9
READ(61,*) solution_sdss_ml(i),errorbar_sdss_ml(i)
END DO
READ(61,*) 
READ(61,*) rms_2mrs_ml
DO i=1,9
READ(61,*) solution_2mrs_ml(i),errorbar_2mrs_ml(i)
END DO
READ(61,*) 
READ(61,*) rms_sdss_mh
DO i=1,10
READ(61,*) solution_sdss_mh(i),errorbar_sdss_mh(i)
END DO
READ(61,*) 
READ(61,*) rms_2mrs_mh
DO i=1,10
READ(61,*) solution_2mrs_mh(i),errorbar_2mrs_mh(i)
END DO
 CLOSE(61)
  
  
  
  
      WRITE(*,*) 'fit data loaded'

      
      
      
 OPEN(61,file='mass/lumdev.txt')
 WRITE(61,*) sigma_lumtot_2mrs_s,sigma_lumtot_2mrs_ml,sigma_lumtot_2mrs_mh 
 WRITE(61,*) sigma_lumtot_sdss_s,sigma_lumtot_sdss_ml,sigma_lumtot_sdss_mh
 CLOSE(61)
  
   OPEN(61,file='mass/lumdev_2MRS.txt')
 WRITE(61,*) sigma_lumtot_2mrs_s,sigma_lumtot_2mrs_ml,sigma_lumtot_2mrs_mh 
 CLOSE(61)
  
   OPEN(61,file='mass/lumdev_SDSS.txt')
 WRITE(61,*) sigma_lumtot_sdss_s,sigma_lumtot_sdss_ml,sigma_lumtot_sdss_mh
 CLOSE(61)
  
  
  
  

DO i=1,n_s_sdss
otherside_sdss_s(i)=solution_sdss_s(1)*lumsum_cor_s_sdss(i)+solution_sdss_s(2)*(lumsum_cor_s_sdss(i)**2)+&
solution_sdss_s(3)*(lumsum_cor_s_sdss(i)**3)+solution_sdss_s(4)*lumdist_s_sdss(i)+&
solution_sdss_s(5)*(lumdist_s_sdss(i)**2)+solution_sdss_s(6)*(lumdist_s_sdss(i)**3)+solution_sdss_s(7)
END DO

DO i=1,n_s_2mrs
otherside_2mrs_s(i)=solution_2mrs_s(1)*lumsum_cor_s_2mrs(i)+solution_2mrs_s(2)*(lumsum_cor_s_2mrs(i)**2)+&
solution_2mrs_s(3)*(lumsum_cor_s_2mrs(i)**3)+solution_2mrs_s(4)*lumdist_s_2mrs(i)+&
solution_2mrs_s(5)*(lumdist_s_2mrs(i)**2)+solution_2mrs_s(6)*(lumdist_s_2mrs(i)**3)+solution_2mrs_s(7)
END DO

DO i=1,n_ml_sdss
otherside_sdss_ml(i)=solution_sdss_ml(1)*lumsum_cor_ml_sdss(i)+solution_sdss_ml(2)*(lumsum_cor_ml_sdss(i)**2)+&
solution_sdss_ml(3)*(lumsum_cor_ml_sdss(i)**3)+solution_sdss_ml(4)*lumdist_ml_sdss(i)+&
solution_sdss_ml(5)*(lumdist_ml_sdss(i)**2)+solution_sdss_ml(6)*(lumdist_ml_sdss(i)**3)+&
solution_sdss_ml(7)*sigma_ml_sdss(i)+solution_sdss_ml(8)*radius_ml_sdss(i)+solution_sdss_ml(9)
END DO

DO i=1,n_ml_2mrs
otherside_2mrs_ml(i)=solution_2mrs_ml(1)*lumsum_cor_ml_2mrs(i)+solution_2mrs_ml(2)*(lumsum_cor_ml_2mrs(i)**2)+&
solution_2mrs_ml(3)*(lumsum_cor_ml_2mrs(i)**3)+solution_2mrs_ml(4)*lumdist_ml_2mrs(i)+&
solution_2mrs_ml(5)*(lumdist_ml_2mrs(i)**2)+solution_2mrs_ml(6)*(lumdist_ml_2mrs(i)**3)+&
solution_2mrs_ml(7)*sigma_ml_2mrs(i)+solution_2mrs_ml(8)*radius_ml_2mrs(i)+solution_2mrs_ml(9)
END DO

DO i=1,n_mh_sdss
otherside_sdss_mh(i)=solution_sdss_mh(1)*lumsum_cor_mh_sdss(i)+solution_sdss_mh(2)*(lumsum_cor_mh_sdss(i)**2)+&
solution_sdss_mh(3)*(lumsum_cor_mh_sdss(i)**3)+solution_sdss_mh(4)*lumdist_mh_sdss(i)+&
solution_sdss_mh(5)*(lumdist_mh_sdss(i)**2)+solution_sdss_mh(6)*(lumdist_mh_sdss(i)**3)+&
solution_sdss_mh(7)*sigma_mh_sdss(i)+solution_sdss_mh(8)*radius_mh_sdss(i)+solution_sdss_mh(9)*n_visible_mh_sdss(i)+&
solution_sdss_mh(10)
END DO

DO i=1,n_mh_2mrs
otherside_2mrs_mh(i)=solution_2mrs_mh(1)*lumsum_cor_mh_2mrs(i)+solution_2mrs_mh(2)*(lumsum_cor_mh_2mrs(i)**2)+&
solution_2mrs_mh(3)*(lumsum_cor_mh_2mrs(i)**3)+solution_2mrs_mh(4)*lumdist_mh_2mrs(i)+&
solution_2mrs_mh(5)*(lumdist_mh_2mrs(i)**2)+solution_2mrs_mh(6)*(lumdist_mh_2mrs(i)**3)+&
solution_2mrs_mh(7)*sigma_mh_2mrs(i)+solution_2mrs_mh(8)*radius_ml_2mrs(i)+solution_2mrs_mh(9)*n_visible_mh_2mrs(i)+&
solution_2mrs_mh(10)
END DO



DO i=1,n_s_sdss
act_sdss_s(i)=.TRUE.
deviation_sdss_s(i)=otherside_sdss_s(i)-h_mass_s_sdss(i)
IF (abs(deviation_sdss_s(i))>1.5D0) THEN
act_sdss_s(i)=.FALSE.
END IF
END DO
DO i=1,n_s_2mrs
act_2mrs_s(i)=.TRUE.
deviation_2mrs_s(i)=otherside_2mrs_s(i)-h_mass_s_2mrs(i)
IF (abs(deviation_2mrs_s(i))>1.5D0) THEN
act_2mrs_s(i)=.FALSE.
END IF
END DO
DO i=1,n_ml_sdss
act_sdss_ml(i)=.TRUE.
deviation_sdss_ml(i)=otherside_sdss_ml(i)-h_mass_ml_sdss(i)
IF (abs(deviation_sdss_ml(i))>1.5D0) THEN
act_sdss_ml(i)=.FALSE.
END IF
END DO
DO i=1,n_ml_2mrs
act_2mrs_ml(i)=.TRUE.
deviation_2mrs_ml(i)=otherside_2mrs_ml(i)-h_mass_ml_2mrs(i)
IF (abs(deviation_2mrs_ml(i))>1.5D0) THEN
act_2mrs_ml(i)=.FALSE.
END IF
END DO
DO i=1,n_mh_sdss
act_sdss_mh(i)=.TRUE.
deviation_sdss_mh(i)=otherside_sdss_mh(i)-h_mass_mh_sdss(i)
IF (abs(deviation_sdss_mh(i))>1.5D0) THEN
act_sdss_mh(i)=.FALSE.
END IF
END DO
DO i=1,n_mh_2mrs
act_2mrs_mh(i)=.TRUE.
deviation_2mrs_mh(i)=otherside_2mrs_mh(i)-h_mass_mh_2mrs(i)
IF (abs(deviation_2mrs_mh(i))>1.5D0) THEN
act_2mrs_mh(i)=.FALSE.
END IF
END DO





WRITE(*,*) 'residuals caluclated'

mass_fit_sum_sdss_s=0.D0
mass_fit_sum_2mrs_s=0.D0
mass_orig_sum_sdss_s=0.D0
mass_orig_sum_2mrs_s=0.D0

mass_fit_sum_sdss_ml=0.D0
mass_fit_sum_2mrs_ml=0.D0
mass_orig_sum_sdss_ml=0.D0
mass_orig_sum_2mrs_ml=0.D0

mass_fit_sum_sdss_mh=0.D0
mass_fit_sum_2mrs_mh=0.D0
mass_orig_sum_sdss_mh=0.D0
mass_orig_sum_2mrs_mh=0.D0




DO i=1,n_s_sdss
mass_fit_sum_sdss_s=mass_fit_sum_sdss_s+(10.D0**otherside_sdss_s(i))
mass_orig_sum_sdss_s=mass_orig_sum_sdss_s+(10.D0**h_mass_s_sdss(i))
END DO

DO i=1,n_ml_sdss
mass_fit_sum_sdss_ml=mass_fit_sum_sdss_ml+(10.D0**otherside_sdss_ml(i))
mass_orig_sum_sdss_ml=mass_orig_sum_sdss_ml+(10.D0**h_mass_ml_sdss(i))
END DO

DO i=1,n_mh_sdss
mass_fit_sum_sdss_mh=mass_fit_sum_sdss_mh+(10.D0**otherside_sdss_mh(i))
mass_orig_sum_sdss_mh=mass_orig_sum_sdss_mh+(10.D0**h_mass_mh_sdss(i))
END DO

DO i=1,n_s_2mrs
mass_fit_sum_2mrs_s=mass_fit_sum_2mrs_s+(10.D0**otherside_2mrs_s(i))
mass_orig_sum_2mrs_s=mass_orig_sum_2mrs_s+(10.D0**h_mass_s_2mrs(i))
END DO

DO i=1,n_ml_2mrs
mass_fit_sum_2mrs_ml=mass_fit_sum_2mrs_ml+(10.D0**otherside_2mrs_ml(i))
mass_orig_sum_2mrs_ml=mass_orig_sum_2mrs_ml+(10.D0**h_mass_ml_2mrs(i))
END DO

DO i=1,n_mh_2mrs
mass_fit_sum_2mrs_mh=mass_fit_sum_2mrs_mh+(10.D0**otherside_2mrs_mh(i))
mass_orig_sum_2mrs_mh=mass_orig_sum_2mrs_mh+(10.D0**h_mass_mh_2mrs(i))
END DO


mass_fit_sum_sdss_tot=mass_fit_sum_sdss_s+mass_fit_sum_sdss_ml+mass_fit_sum_sdss_mh
mass_fit_sum_2mrs_tot=mass_fit_sum_2mrs_s+mass_fit_sum_2mrs_ml+mass_fit_sum_2mrs_mh
mass_orig_sum_sdss_tot=mass_orig_sum_sdss_s+mass_orig_sum_sdss_ml+mass_orig_sum_sdss_mh
mass_orig_sum_2mrs_tot=mass_orig_sum_2mrs_s+mass_orig_sum_2mrs_ml+mass_orig_sum_2mrs_mh


OPEN(61,file='mass/sums.txt')
WRITE(61,*) mass_fit_sum_sdss_s,mass_orig_sum_sdss_s
WRITE(61,*) mass_fit_sum_sdss_ml,mass_orig_sum_sdss_ml
WRITE(61,*) mass_fit_sum_sdss_mh,mass_orig_sum_sdss_mh
WRITE(61,*) mass_fit_sum_2mrs_s,mass_orig_sum_2mrs_s
WRITE(61,*) mass_fit_sum_2mrs_ml,mass_orig_sum_2mrs_ml
WRITE(61,*) mass_fit_sum_2mrs_mh,mass_orig_sum_2mrs_mh
WRITE(61,*)
WRITE(61,*) mass_fit_sum_sdss_tot,mass_orig_sum_sdss_tot
WRITE(61,*) mass_fit_sum_2mrs_tot,mass_orig_sum_2mrs_tot
WRITE(61,*)
WRITE(61,*) (mass_fit_sum_sdss_tot/m_expected),(mass_orig_sum_sdss_tot/m_expected)
WRITE(61,*) (mass_fit_sum_2mrs_tot/m_expected),(mass_orig_sum_2mrs_tot/m_expected)
 CLOSE(61)





WRITE(*,*) 'masses summed up'




DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_ml_2mrs
help_x=NINT((h_mass_ml_2mrs(i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_ml.txt')
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

DO i=1,n_mh_2mrs
help_x=NINT((h_mass_mh_2mrs(i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/deviation_2mrs_mh.txt')
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

DO i=1,n_s_2mrs
help_x=NINT((h_mass_s_2mrs(i)-10.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
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

DO i=1,n_ml_sdss
help_x=NINT((h_mass_ml_sdss(i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_ml.txt')
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

DO i=1,n_mh_sdss
help_x=NINT((h_mass_mh_sdss(i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/deviation_sdss_mh.txt')
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

DO i=1,n_s_sdss
help_x=NINT((h_mass_s_sdss(i)-10.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
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

 
 
 WRITE(*,*) 'mass residuals written'
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_ml_2mrs
help_x=NINT((lumsum_cor_ml_2mrs(i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_ml.txt')
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

DO i=1,n_mh_2mrs
help_x=NINT((lumsum_cor_mh_2mrs(i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_2mrs_mh.txt')
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

DO i=1,n_s_2mrs
help_x=NINT((lumsum_cor_s_2mrs(i)-7.D0)*25.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
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

DO i=1,n_ml_sdss
help_x=NINT((lumsum_cor_ml_sdss(i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_ml.txt')
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

DO i=1,n_mh_sdss
help_x=NINT((lumsum_cor_mh_sdss(i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_lum_sdss_mh.txt')
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

DO i=1,n_s_sdss
help_x=NINT((lumsum_cor_s_sdss(i)-7.D0)*25.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<150)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
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

 

 
 
 
 
  WRITE(*,*) 'luminosity residuals written'
 
 
 
 
 
 
 
 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_ml_2mrs
help_x=NINT((lumdist_ml_2mrs(i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_ml.txt')
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

DO i=1,n_mh_2mrs
help_x=NINT((lumdist_mh_2mrs(i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_2mrs_mh.txt')
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

DO i=1,n_s_2mrs
help_x=NINT((lumdist_s_2mrs(i)-0.5D0)*50.D0)
help_y=NINT((deviation_2mrs_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
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

DO i=1,n_ml_sdss
help_x=NINT((lumdist_ml_sdss(i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_ml.txt')
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

DO i=1,n_mh_sdss
help_x=NINT((lumdist_mh_sdss(i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_dist_sdss_mh.txt')
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

DO i=1,n_s_sdss
help_x=NINT((lumdist_s_sdss(i)-0.5D0)*50.D0)
help_y=NINT((deviation_sdss_s(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<250)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
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

 

 
 
 
  WRITE(*,*) 'distance residuals written'
 
 
 

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_ml_2mrs
help_x=NINT((m_dyn_ml_2mrs(i)-6.D0)*25.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_2mrs_ml.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+6.D0
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

DO i=1,n_mh_2mrs
help_x=NINT((m_dyn_mh_2mrs(i)-6.D0)*25.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_2mrs_mh.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+6.D0
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

DO i=1,n_ml_sdss
help_x=NINT((m_dyn_ml_sdss(i)-6.D0)*25.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_sdss_ml.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+6.D0
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

DO i=1,n_mh_sdss
help_x=NINT((m_dyn_mh_sdss(i)-6.D0)*25.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_mdyn_sdss_mh.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/25.D0)+6.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 

 
 
 
  WRITE(*,*) 'dynamical mass residuals written'
 

 
  
  
  
  

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_ml_2mrs
help_x=NINT((sigma_ml_2mrs(i)-1.D0)*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_sigma_2mrs_ml.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.D0
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

DO i=1,n_mh_2mrs
help_x=NINT((sigma_mh_2mrs(i)-1.D0)*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_sigma_2mrs_mh.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.D0
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

DO i=1,n_ml_sdss
help_x=NINT((sigma_ml_sdss(i)-1.D0)*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_sigma_sdss_ml.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.D0
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

DO i=1,n_mh_sdss
help_x=NINT((sigma_mh_sdss(i)-1.D0)*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_sigma_sdss_mh.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)+1.D0
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 

 
 
 
  WRITE(*,*) 'velocity dispersion residuals written'
 

 

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_ml_2mrs
! WRITE(*,*) radius_ml_2mrs(i)
help_x=NINT((radius_ml_2mrs(i))*50.D0)
help_y=NINT((deviation_2mrs_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_rad_2mrs_ml.txt')
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

DO i=1,n_mh_2mrs
help_x=NINT((radius_mh_2mrs(i))*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_rad_2mrs_mh.txt')
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

DO i=1,n_ml_sdss
help_x=NINT((radius_ml_sdss(i))*50.D0)
help_y=NINT((deviation_sdss_ml(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_rad_sdss_ml.txt')
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

DO i=1,n_mh_sdss
help_x=NINT((radius_mh_sdss(i))*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
IF ((help_x>0).AND.(help_x<200)) THEN
IF ((help_y>0).AND.(help_y<150)) THEN
mapmap(help_x,help_y)=mapmap(help_x,help_y)+1.D0
END IF
END IF
END DO

OPEN(61,file='mass/res_rad_sdss_mh.txt')
DO i=0,200
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 

 
 
 
  WRITE(*,*) 'radius residuals written'
 

 
  
  
  
  
  
  
  
  
  

  
  
  
  
  

DO i=0,500
DO ii=0,500
mapmap(i,ii)=0
END DO
END DO

DO i=1,n_mh_2mrs
help_x=NINT((n_visible_mh_2mrs(i))*50.D0)
help_y=NINT((deviation_2mrs_mh(i)+1.5D0)*50.D0)
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

OPEN(61,file='mass/res_nvis_2mrs_mh.txt')
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

DO i=1,n_mh_sdss
help_x=NINT((n_visible_mh_sdss(i))*50.D0)
help_y=NINT((deviation_sdss_mh(i)+1.5D0)*50.D0)
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

OPEN(61,file='mass/res_nvis_sdss_mh.txt')
DO i=0,150
DO ii=0,150
help_x2=(DBLE(i)/50.D0)
help_y2=(DBLE(ii)/50.D0)-1.5D0
WRITE(61,*) help_x2,help_y2,mapmap(i,ii)
END DO
END DO
 CLOSE(61)


 
 
 
 
 
  WRITE(*,*) 'number residuals written'
 
 
 
 OPEN(61,file='mass/all_residuals_sdss_s.txt')
DO i=1,n_s_sdss
WRITE(61,*) deviation_sdss_s(i),h_mass_s_sdss(i),lumsum_cor_s_sdss(i),lumdist_s_sdss(i),&
lumsum_s_sdss(i),lumsum_all_s_sdss(i)
END DO
 CLOSE(61)

 
  OPEN(61,file='mass/all_residuals_sdss_ml.txt')
DO i=1,n_ml_sdss
WRITE(61,*) deviation_sdss_ml(i),h_mass_ml_sdss(i),lumsum_cor_ml_sdss(i),lumdist_ml_sdss(i),&
m_dyn_ml_sdss(i),sigma_ml_sdss(i),radius_ml_sdss(i),lumsum_ml_sdss(i),lumsum_all_ml_sdss(i)
END DO
 CLOSE(61)
 
 
   OPEN(61,file='mass/all_residuals_sdss_mh.txt')
DO i=1,n_mh_sdss
WRITE(61,*) deviation_sdss_mh(i),h_mass_mh_sdss(i),lumsum_cor_mh_sdss(i),lumdist_mh_sdss(i),m_dyn_mh_sdss(i),&
sigma_mh_sdss(i),radius_mh_sdss(i),n_visible_mh_sdss(i),lumsum_mh_sdss(i),lumsum_all_mh_sdss(i)
END DO
 CLOSE(61)
  
 
 OPEN(61,file='mass/all_residuals_2mrs_s.txt')
DO i=1,n_s_2mrs
WRITE(61,*) deviation_2mrs_s(i),h_mass_s_2mrs(i),lumsum_cor_s_2mrs(i),lumdist_s_2mrs(i),&
lumsum_s_2mrs(i),lumsum_all_s_2mrs(i)
END DO
 CLOSE(61)

 
  OPEN(61,file='mass/all_residuals_2mrs_ml.txt')
DO i=1,n_ml_2mrs
WRITE(61,*) deviation_2mrs_ml(i),h_mass_ml_2mrs(i),lumsum_cor_ml_2mrs(i),lumdist_ml_2mrs(i),&
m_dyn_ml_2mrs(i),sigma_ml_2mrs(i),radius_ml_2mrs(i),lumsum_ml_2mrs(i),lumsum_all_ml_2mrs(i)
END DO
 CLOSE(61)
 
 
   OPEN(61,file='mass/all_residuals_2mrs_mh.txt')
DO i=1,n_mh_2mrs
WRITE(61,*) deviation_2mrs_mh(i),h_mass_mh_2mrs(i),lumsum_cor_mh_2mrs(i),lumdist_mh_2mrs(i),m_dyn_mh_2mrs(i),&
sigma_mh_2mrs(i),radius_mh_2mrs(i),n_visible_mh_2mrs(i),lumsum_mh_2mrs(i),lumsum_all_mh_2mrs(i)
END DO
 CLOSE(61)
 
 
 
 
 
 
 
   OPEN(61,file='mass/mfit_SDSS_s.txt')
WRITE(61,*) rms_sdss_s
DO i=1,7
WRITE(61,*) solution_sdss_s(i),errorbar_sdss_s(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/mfit_SDSS_ml.txt')
WRITE(61,*) rms_sdss_ml
DO i=1,9
WRITE(61,*) solution_sdss_ml(i),errorbar_sdss_ml(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/mfit_SDSS_mh.txt')
WRITE(61,*) rms_sdss_mh
DO i=1,10
WRITE(61,*) solution_sdss_mh(i),errorbar_sdss_mh(i)
END DO
 CLOSE(61)
 

    OPEN(61,file='mass/mfit_2MRS_s.txt')
WRITE(61,*) rms_2mrs_s
DO i=1,7
WRITE(61,*) solution_2mrs_s(i),errorbar_2mrs_s(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/mfit_2MRS_ml.txt')
WRITE(61,*) rms_2mrs_ml
DO i=1,9
WRITE(61,*) solution_2mrs_ml(i),errorbar_2mrs_ml(i)
END DO
 CLOSE(61)
 
  OPEN(61,file='mass/mfit_2MRS_mh.txt')
WRITE(61,*) rms_2mrs_mh
DO i=1,10
WRITE(61,*) solution_2mrs_mh(i),errorbar_2mrs_mh(i)
END DO
 CLOSE(61)
 

 

 
 
 
 WRITE(*,*) 'all residuals written'



END PROGRAM 








