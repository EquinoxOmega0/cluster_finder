PROGRAM clusterfinder
! declaration of variables
IMPLICIT NONE
integer(kind=8) :: io_err,n,i,ii

integer(kind=8), allocatable :: members_clustercore(:),clusterindex(:),clustergalnmember(:)
integer(kind=8), allocatable :: id1(:),id2(:),numberofgal(:)
integer :: hwrite1,hwrite2

character(200), allocatable :: mrs_id(:)
character(40) :: str1,str2,pref1,pref2,signum

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme 2MRS_ID_merger started'
WRITE(*,*) '============================================================'



! get length of file
OPEN(50,file='catalogues/cluster_members_2MRS.txt')
io_err=0
n=0
DO WHILE (io_err==0) 
READ(50,*,iostat=io_err) 
n=n+1
END DO
 CLOSE(50)
n=n-1

WRITE(*,*) n,'galaxies will be used'

allocate(clusterindex(1:n))
allocate(id1(1:n))
allocate(id2(1:n))
allocate(numberofgal(1:n))
allocate(mrs_id(1:n))
 
OPEN(50,file='catalogues/cluster_members_2MRS.txt')
DO i=1,n
   READ(50,*)  numberofgal(i),clusterindex(i),id1(i),id2(i)
END DO
 CLOSE(50)
 
 
DO i=1,n

WRITE(str1,*) id1(i)
str1=adjustl(str1)

WRITE(str2,*) abs(id2(i))
str2=adjustl(str2)

pref1=''
pref2=''

IF (id1(i)<10000000) THEN
pref1='0'
END IF
IF (id1(i)<1000000) THEN
pref1='00'
END IF
IF (id1(i)<100000) THEN
pref1='000'
END IF
IF (id1(i)<10000) THEN
pref1='0000'
END IF
IF (id1(i)<1000) THEN
pref1='00000'
END IF
IF (id1(i)<100) THEN
pref1='000000'
END IF
IF (id1(i)<10) THEN
pref1='0000000'
END IF


IF (id2(i)<0) THEN
signum='-'
ELSE
signum='+'
END IF


IF (abs(id2(i))<1000000) THEN
pref2='0'
END IF
IF (abs(id2(i))<100000) THEN
pref2='00'
END IF
IF (abs(id2(i))<10000) THEN
pref2='000'
END IF
IF (abs(id2(i))<1000) THEN
pref2='0000'
END IF
IF (abs(id2(i))<100) THEN
pref2='00000'
END IF
IF (abs(id2(i))<10) THEN
pref2='000000'
END IF

mrs_id(i)=TRIM(adjustl(pref1))//TRIM(adjustl(str1))//TRIM(adjustl(signum))//TRIM(adjustl(pref2))//TRIM(adjustl(str2))

END DO

OPEN(50,file='catalogues/cluster_members_2MRSb.txt')
DO i=1,n
 hwrite1=numberofgal(i)
 hwrite2=clusterindex(i)
   WRITE(50,*)  hwrite1,hwrite2,mrs_id(i)
END DO
 CLOSE(50)

WRITE(*,*) '============================================================'
WRITE(*,*) '    programme complete'
WRITE(*,*) '============================================================'




END PROGRAM


 
 
 
 
 