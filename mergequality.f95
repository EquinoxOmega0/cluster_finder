PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
character(20) :: nm,ecs
integer :: u,exec_counter,i,ii
double precision :: dummy,qmax
double precision, dimension(1:8) :: S_tot,b_factor,R_factor,lum_power
double precision, dimension(1:31) :: S_tot_av,b_factor0,R_factor0,lum_power0
 
 
 
 DO i=1,31
 DO ii=1,8
 
 
WRITE(ecs,*) i
ecs=adjustl(ecs)

WRITE(nm,*) ii
nm=adjustl(nm)
 
 

OPEN(50,file='quality_SDSS/group_cost_function_'//TRIM(nm)//'_'//TRIM(ecs)//'.txt')
READ(50,*) b_factor(ii),R_factor(ii),lum_power(ii)
READ(50,*) dummy,dummy,dummy
READ(50,*) dummy,dummy,dummy
READ(50,*) S_tot(ii)
 CLOSE(50)

 
 END DO
 
  
S_tot_av(i)=0.D0
DO ii=1,8
S_tot_av(i)=S_tot_av(i)+S_tot(ii)
END DO
S_tot_av(i)=S_tot_av(i)/8.D0

b_factor0(i)=b_factor(1)
R_factor0(i)=R_factor(1)
lum_power0(i)=lum_power(1)
 
 END DO
 qmax=0.D0
 
 
 OPEN(77,file='quality_SDSS/best_S_evolution.txt')
  DO i=1,31
  IF (S_tot_av(i)>qmax) THEN
WRITE(77,*) b_factor0(i),R_factor0(i),lum_power0(i),S_tot_av(i),i
qmax=S_tot_av(i)
END IF
 END DO
  CLOSE(77)

  
 OPEN(77,file='quality_SDSS/best_S_evolution_unfiltered.txt')
  DO i=1,31
WRITE(77,*) b_factor0(i),R_factor0(i),lum_power0(i),S_tot_av(i),i
 END DO
  CLOSE(77)
 
 
END PROGRAM


 
 
 
 
 
 
 
 
 


 

 