! Must define these since decs.h is not included...
#define TRUE 1
#define FALSE 0

module table_module

  ! number of elements in AB table
  integer, parameter :: numel = 16

  ! table dimensions
  integer, parameter :: nt = 300
  integer, parameter :: nr = 300

  ! LS stellar collapse requires nt=300, nr=300, ny=50 (100MEV)
  integer, parameter :: ny = 50
  double precision :: r1  =  1.d0
  double precision :: r2  =  15.5d0
  double precision :: t1  = -2.0d0
  double precision :: t2  =  2.0d0
  double precision :: t12 = -2.0d0
  double precision :: t22 =  2.0d0
  double precision :: y1  =  0.035d0
  double precision :: y2  =  0.56d0

  ! LS stellar collapse requires nt=300, nr=300, ny=50 (HI TEMP)
!  integer, parameter :: ny = 50
!  double precision :: r1  =  1.d0  !  4.0d0
!  double precision :: r2  =  15.0d0
!  double precision :: t1  = -2.0d0
!  double precision :: t2  =  1.82d0 !0.8d0
!  double precision :: t12 = -2.0d0
!  double precision :: t22 =  1.82d0
!  double precision :: y1  =  0.035d0
!  double precision :: y2  =  0.56d0

  ! LS stellar collapse requires nt=300, nr=300, ny=50
!  integer, parameter :: ny = 50
!  double precision :: r1  =  1.d0  !  4.0d0
!  double precision :: r2  =  15.0d0
!  double precision :: t1  = -2.0d0
!  double precision :: t2  =  1.7d0 !0.8d0
!  double precision :: t12 = -2.0d0
!  double precision :: t22 =  1.7d0
!  double precision :: y1  =  0.035d0
!  double precision :: y2  =  0.56d0
 
  ! Shen 300 300 54, t lower 0.01 MeV

!  integer, parameter :: ny = 54
!  double precision :: r1 = 1.0d0
!  double precision :: r2  =  15.0d0
!  double precision :: t1  = -2.0d0  ! THIS IS FOR SHEN OCT 16
!  double precision :: t2 = 1.7d0
!  double precision :: t12 = -2.0d0  ! THIS IS FOR SHEN OCT 16
!  double precision :: t22 = 1.7d0
!  double precision :: y1 = 0.03d0
!  double precision :: y2 = 0.56d0

  double precision,allocatable,save :: table(:,:,:,:)  !table(numel,ny,nr,nt)
  double precision,allocatable :: rtable(:,:,:,:)      !rtable(numel,ny,nr,nt)
  double precision,allocatable :: ttable(:,:)          !ttable(numel,nt)  project for a specifc rho,ye

 contains

subroutine collapse_init(eos_name)

! use eos_module

  implicit none

  integer i, jy, jr, jt, irec, irecl, tndx, tndx0, ipass
  character(*) :: eos_name
  integer :: ncpu,mype,info,ierr

!  double precision, allocatable :: table1d(:)

#if(USE_MPI==TRUE)
  include 'mpif.h'
  ! Init MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,info)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,info)
#else
  ncpu = 1
  mype = 0
#endif

!--------------------------------------------
!  table(1,y,r,t)  :: energy per baryon      |
!  table(2,y,r,t)  :: pressure               |
!  table(3,y,r,t)  :: entropy per baryon     |
!  table(4,y,r,t)  :: cv                     |
!  table(5,y,r,t)  :: xn                     |
!  table(6,y,r,t)  :: xp                     |
!  table(7,y,r,t)  :: xa                     |
!  table(8,y,r,t)  :: xh                     |
!  table(9,y,r,t)  :: za                     |
!  table(10,y,r,t) :: aw                     |
!  table(11,y,r,t) :: muhat                  |
!  table(12,y,r,t) :: gamma!                 |
!  table(13,y,r,t) :: dhy                    |
!  table(14,y,r,t) :: zht                    |
!  table(15,y,r,t) :: mue                    |
!  table(16,y,r,t) :: dpde                   |
!--------------------------------------------

  irecl = 8 * ny * nr * nt
!  allocate(table1d(numel*nt*nr*ny))
  allocate(table(numel,ny,nr,nt), stat=ierr)
  if (ierr /= 0) write (*,*) '[collapse_init]:  proc',mype,'failed to allocate table!'
!  allocate(ttable(numel,nt))
!  allocate(rtable(numel,ny,nr,nt))

  if (mype==0) then
    open(10,file=eos_name, &
        form='unformatted', access='direct', recl=irecl, &
        err=1010)
    goto 1020
  1010 write(*,*) "failed to open file"
  1020 continue

    do irec = 1, numel
       read(10,rec=irec)(((table(irec,jy,jr,jt), jt = 1,nt),jr = 1, nr), jy = 1, ny)
    enddo


    do jt=1,nt
       do jr=1,nr
          do jy=1,ny
!             table(1,jy,jr,jt) = log(table(1,jy,jr,jt))
             table(2,jy,jr,jt) = log(table(2,jy,jr,jt))
          enddo
       enddo
    enddo
    close(unit=10)
!  write(5,*)
!  write(5,*)'finished reading eos table'
!  write(5,*)'table(numel,1,1,1) = ',table(numel,1,1,1)
!  write(5,*)'table(1,ny,1,1) = ',table(1,ny,1,1)
!  write(5,*)'table(1,1,nr,1) = ',table(1,1,nr,1)
!  write(5,*)'table(1,1,1,nt) = ',table(1,1,1,nt)
!  write(5,*)
!  write(5,*)
  endif

#if(USE_MPI==TRUE)
    ! Broadcast EOS table
    call MPI_BCAST(table,size(table),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
#endif


end subroutine collapse_init

subroutine get_numel(num_elements)

  integer, intent(out) :: num_elements
  num_elements = numel

end subroutine get_numel
 

!----------------------------------------------------------------------

subroutine finde(f,nelem,temp,rho,ye,jy_out,jq_out,jr_out,input,pt_index)

  implicit none

  integer, intent(in) :: nelem,input
  double precision, intent(inout) :: f(nelem)
  double precision, intent(in) :: rho, ye
  double precision, intent(inout) :: temp
  integer, intent(out) :: jy_out, jq_out, jr_out
  integer,optional,intent(in) :: pt_index(:)

  integer jt,jr,jy,nn,jtp,jrp,jyp
  double precision :: rl,tl,drho,dtemp,dye
  double precision :: dr,dt,dy
  double precision :: coeff(4)
  double precision :: elo,ehi,tlo,thi

  rl = dlog10(rho)
  tl = dlog10(temp)

  dr = (r2-r1)/(nr-1)
  dt = (t2-t1)/(nt-1)
  dy = (y2-y1)/(ny-1)

  jr = 1 + int((rl-r1)/dr)
  jt = 1 + int((tl-t1)/dt)
  jy = 1 + int((ye-y1)/dy)

  jyp = jy+1
  jrp = jr+1

  drho = (rl - (r1+(jr-1)*dr))/dr
  dye = (ye - (y1+(jy-1)*dy))/dy

  if (jy.eq.ny) then
     jyp = jy
     dye = 0.d0
  endif
  if (jr.eq.nr) then
     jrp = jr
     drho = 0.d0
  endif

  ! project into 1-D table
  coeff(1) = (1.d0 - drho)*(1.d0 - dye)
  coeff(2) = dye*(1.d0 - drho)
  coeff(3) = (1.d0 - dye)*drho
  coeff(4) = dye*drho

  !ttable(:,:) = coeff(1)*table(:,jy,jr,:)
  !ttable(:,:) = ttable(:,:) + coeff(2)*table(:,jyp,jr,:)
  !ttable(:,:) = ttable(:,:) + coeff(3)*table(:,jy,jrp,:)
  !ttable(:,:) = ttable(:,:) + coeff(4)*table(:,jyp,jrp,:)

  ! now find table entries that bound f(1), the energy we're after
  ehi = -1d99
  do
    elo = coeff(1)*table(1,jy,jr,jt)
    elo = elo + coeff(2)*table(1,jyp,jr,jt)
    elo = elo + coeff(3)*table(1,jy,jrp,jt)
    elo = elo + coeff(4)*table(1,jyp,jrp,jt)
    if (elo.lt.f(1) .or. jt.eq.1) then
      exit
    endif
    ehi = elo
    jt = jt - 1
  enddo

  if (elo.gt.f(1)) then
    write(*,*) "Off the bottom of the table in finde", rho, temp, ye, f(1), elo
  endif

  if (ehi.lt.-1d98) then
    do
      ehi = coeff(1)*table(1,jy,jr,jt+1)
      ehi = ehi + coeff(2)*table(1,jyp,jr,jt+1)
      ehi = ehi + coeff(3)*table(1,jy,jrp,jt+1)
      ehi = ehi + coeff(4)*table(1,jyp,jrp,jt+1)
      if (ehi.gt.f(1) .or. jt.eq.nt-1) then
        exit
      endif
      elo = ehi
      jt = jt + 1
    enddo
  endif

  if (ehi.lt.f(1)) then
    write(*,*) "Off the top of the table in finde", rho, temp, ye, f(1), ehi
    call exit(1)
  endif

  ! now use a linear model.  must be linear so it's invertible
  tlo = t1 + (jt-1)*dt
  thi = tlo + dt
  ! dtemp = de = (e-e0)/(e1-e0)
  dtemp = (f(1) - elo)/(ehi-elo)
  temp = tlo + (thi-tlo)*dtemp

  if (dtemp .lt. 0.d0 .or. dtemp .gt. 1.d0) then
    write(*,*) "dtemp out of bounds ", dtemp
    call exit(1)
  endif

  jtp = jt+1
  !f(:) = ttable(:,jt)*(1-dtemp) + ttable(:,jt+1)*dtemp
  f(:) = coeff(1)*(1-dtemp)*table(:,jy,jr,jt)
  f(:) = f(:) + coeff(2)*(1-dtemp)*table(:,jyp,jr,jt)
  f(:) = f(:) + coeff(3)*(1-dtemp)*table(:,jy,jrp,jt)
  f(:) = f(:) + coeff(4)*(1-dtemp)*table(:,jyp,jrp,jt)
  f(:) = f(:) + coeff(1)*dtemp*table(:,jy,jr,jtp)
  f(:) = f(:) + coeff(2)*dtemp*table(:,jyp,jr,jtp)
  f(:) = f(:) + coeff(3)*dtemp*table(:,jy,jrp,jtp)
  f(:) = f(:) + coeff(4)*dtemp*table(:,jyp,jrp,jtp)

  f(2) = exp(f(2))
  temp = 10**temp

  jy_out = jy
  jq_out = jt
  jr_out = jr

end subroutine finde


subroutine findthis(f,nelem,temp,rho,ye,jy_out,jq_out,jr_out,input,pt_index)

  implicit none

  integer, intent(in) :: nelem,input
  double precision, intent(inout) :: f(nelem)
  double precision, intent(in) :: temp, rho, ye
  integer, intent(out) :: jy_out, jq_out, jr_out
  integer,optional,intent(in) :: pt_index(:)

  integer jt,jr,jy,nn,jtp,jrp,jyp
  double precision :: rl,tl,drho,dtemp,dye
  double precision :: dr,dt,dy
  double precision :: coeff(8)

  rl = dlog10(rho)
  tl = dlog10(temp)

  dr = (r2-r1)/(nr-1)
  dt = (t2-t1)/(nt-1)
  dy = (y2-y1)/(ny-1)

  jr = 1 + int((rl-r1)/dr)
  jt = 1+ int((tl-t1)/dt)
  jy = 1 + int((ye-y1)/dy)

  jyp = jy+1
  jrp = jr+1
  jtp = jt+1

  drho = (rl - (r1+(jr-1)*dr))/dr
  dtemp = (tl - (t1+(jt-1)*dt))/dt
  dye = (ye - (y1+(jy-1)*dy))/dy

  if (jy.eq.ny) then
     jyp = jy
     dye = 0.d0
  endif
  if (jr.eq.nr) then
     jrp = jr
     drho = 0.d0
  endif
  if (jt.eq.nt) then
     jtp = jt
     dtemp = 0.d0
  endif

  coeff(1) = (1.d0 - drho)*(1.d0 - dtemp)*(1.d0 - dye)
  coeff(2) = (1.d0 - drho)*(1.d0 - dtemp)*dye
  coeff(3) = drho*(1.d0 - dtemp)*(1.d0 - dye)
  coeff(4) = drho*(1.d0 - dtemp)*dye
  coeff(5) = (1.d0 - drho)*dtemp*(1.d0 - dye)
  coeff(6) = (1.d0 - drho)*dtemp*dye
  coeff(7) = drho*dtemp*(1.d0 - dye)
  coeff(8) = drho*dtemp*dye

  f(:) = coeff(1)*table(:,jy,jr,jt)
  f(:) = f(:) + coeff(2)*table(:,jyp,jr,jt)
  f(:) = f(:) + coeff(3)*table(:,jy,jrp,jt)
  f(:) = f(:) + coeff(4)*table(:,jyp,jrp,jt)
  f(:) = f(:) + coeff(5)*table(:,jy,jr,jtp)
  f(:) = f(:) + coeff(6)*table(:,jyp,jr,jtp)
  f(:) = f(:) + coeff(7)*table(:,jy,jrp,jtp)
  f(:) = f(:) + coeff(8)*table(:,jyp,jrp,jtp)

!  f(1) = exp(f(1))
  f(2) = exp(f(2))

!  f(12) = (table(12,jy,jr,jt)    + table(12,jyp,jr,jt) &
!          + table(12,jy,jrp,jt)  + table(12,jrp,jrp,jt) &
!          + table(12,jy,jr,jtp)  + table(12,jyp,jr,jtp) &
!          + table(12,jy,jrp,jtp) + table(12,jyp,jrp,jtp))/8.d0

!  do nn=1,nelem
!
!     f(nn) =  coeff(1)*table(nn,jy,jr,jt) &
!            + coeff(2)*table(nn,jyp,jr,jt) &
!            + coeff(3)*table(nn,jy,jrp,jt) &
!            + coeff(4)*table(nn,jyp,jrp,jt) &
!            + coeff(5)*table(nn,jy,jr,jtp) &
!            + coeff(6)*table(nn,jyp,jr,jtp) &
!            + coeff(7)*table(nn,jy,jrp,jtp) &
!            + coeff(8)*table(nn,jyp,jrp,jtp)
!
!     if (nn.eq.5.or.nn.eq.6.or.nn.eq.7.or.nn.eq.8) then
!          f(nn)  = min(max(f(nn),0.d0),1.d0)
!     end if
!  enddo

  jy_out = jy
  jq_out = jt
  jr_out = jr

end subroutine findthis


subroutine findthis_quad(f,nelem,temp,rho,ye,jy_out,jq_out,jr_out,input,pt_index)

! use eos_module

  implicit none

  double precision, intent(inout) :: f(:)
  integer         , intent(in   ) :: nelem,input
  integer,optional, intent(in   ) :: pt_index(:)
  double precision, intent(in   ) :: temp,rho,ye
  integer         , intent(  out) :: jy_out,jq_out,jr_out

  double precision :: rl,tl,rfrac,yfrac
  double precision :: q,p,delty,yl0,yl1,yl2,told
  double precision :: ql,alpha,beta,delta,q0,q1,q2
  double precision :: rl0,rl1,rl2,tl0,tl1,tl2,t
  integer          :: jy,jq,jr,nn
  double precision :: coeff(10)
  double precision :: pq,pp,qq
  double precision :: dy10,dy20,dy21,yefac0,yefac1,yefac2,yefac3

  rl    = dlog10(rho)
  tl    = dlog10(temp)

! if ( (rl .lt. r1) .or. (rl .gt. r2) ) then
!    print *,'LOG(DENSITY) OUT OF BOUNDS ',rl
!    print *,'LIMITS ARE: R1 R2     ',r1,r2
!    print *,'FROM CALL WITH INPUT  ',input
!    if (present(pt_index)) &
!       print *,'AT POINT              ',pt_index(:)
!    stop
! end if

  if ( (tl .lt. t1) .or. (tl .gt. t2) ) then
     print *,'TEMP OUT OF BOUNDS   ',tl
     print *,'LIMITS ARE: T1 T2    ',t1,t2
     print *,'FROM CALL WITH INPUT ',input
     if (present(pt_index)) &
        print *,'AT POINT              ',pt_index(:)
     stop
  end if

! if ( (ye .lt. y1) .or. (ye .gt. y2) ) then
!    print *,'YE   OUT OF BOUNDS   ',ye
!    print *,'LIMITS ARE: Y1 Y2    ',y1,y2
!    print *,'FROM CALL WITH INPUT ',input
!    if (present(pt_index)) &
!       print *,'AT POINT              ',pt_index(:)
!    stop
! end if

  ! Checking limits on rho
  rfrac = (rl-r1)/(r2-r1)
  delta = dble(nr-1) * rfrac
  jr = 1 + int(delta)
  jr = max(2,min(jr,nr-1))
  p  = delta - dble(jr-1)

  if (p .lt. 0.d0) p = 0.d0
  if (p .gt. 1.d0) p = 1.d0

! if (p .lt. 0.d0 .or. p .gt. 1.d0) then
!    print *,'P OUT OF BOUNDS ',p
!    print *,'LOG(DEN) LOG(TEMP) ',rl,tl
!    print *,'FROM CALL WITH INPUT ',input
!    stop
! end if

  ! Checking limits on temp (given rho)
  alpha = t1    +            (t12-t1)  * rfrac
  beta  = t2-t1 + ((t22-t12)-( t2-t1)) * rfrac
  ql    = (tl - alpha)/beta
  jq    = 1 + idint(dble(nt-1)*ql)
  jq    = max(2,min(jq,nt-1))
  q     = dble(nt-1)*ql - dble(jq-1)

  if (q .lt. 0.d0) q = 0.d0
  if (q .gt. 1.d0) q = 1.d0

! if (q .lt. 0.d0 .or. q .gt. 1.d0) then
!    print *,'Q OUT OF BOUNDS ',q
!    print *,'LOG(DEN) LOG(TEMP) ',rl,tl
!    print *,'FROM CALL WITH INPUT ',input
!    print *,'    ********      '
!    print *,'ALPHA, BETA, QL   ',alpha, beta, ql
!    print *,'JQ BEFORE MIN/MAX ', 1 + idint(dble(nt-1)*ql)
!    print *,'    ********      '
!    stop
! end if
  
  ! Set the Y-related indices
  yfrac = (ye-y1)/(y2-y1)
  delty = dble(ny-1) * yfrac
  jy = 1 + int(delty)
  jy = max(2,min(jy,ny-1))
  yl0 = y1+(y2-y1)*dble(jy-1-1)/dble(ny-1)
  yl1 = y1+(y2-y1)*dble(jy-1  )/dble(ny-1)
  yl2 = y1+(y2-y1)*dble(jy+1-1)/dble(ny-1)


  !Calculate coefficients that will be common among the 
  pq = p*q
  pp = p*p
  qq = q*q
  
  coeff(1) = 1.D0-P-Q+PQ
  coeff(2) = P-PQ
  coeff(3) = Q-PQ
  coeff(4) = PQ
  coeff(5) = 0.5D0*(QQ-Q)
  coeff(6) = 0.5D0*(PP-P)
  coeff(7) = 1.D0+PQ-PP-QQ
  coeff(8) = 0.5D0*(PP-2.D0*PQ+P)
  coeff(9) = 0.5D0*(QQ-2.D0*PQ+Q)
  coeff(10) = coeff(4)

  !coeff(1) = (1.d0-p)*(1.d0-q)
  !coeff(2) = p*(1.d0-q)
  !coeff(3) = (1.d0-p)*q
  !coeff(4) = p*q
  !coeff(5) = 0.5d0*q*(q-1.d0)
  !coeff(6) = 0.5d0*p*(p-1.d0)
  !coeff(7) = (1.d0+p*q-p*p-q*q)
  !coeff(8) = 0.5d0*p*(p-2.d0*q+1.d0)
  !coeff(9) = 0.5d0*q*(q-2.d0*p+1.d0)
  !coeff(10) = coeff(4)
  !factors for quadratic interpolation
  dy10 = 1.d0/(yl1-yl0)
  dy20 = 1.d0/(yl2-yl0)
  dy21 = 1.d0/(yl2-yl1)
  !yefac0 = (ye-yl1)*dy21
  yefac1 = (ye-yl1)*(ye-yl2)*dy10*dy20
  yefac2 = -(ye-yl0)*(ye-yl2)*dy10*dy21
  yefac3 = (ye-yl0)*(ye-yl1)*dy20*dy21


  do nn= 1,nelem

   if (nn.eq.20) then

     q0 = coeff(1)*table(nn,jy-1,jr  ,jq  ) &
         + coeff(2)*table(nn,jy-1,jr+1,jq  ) &
         + coeff(3)*table(nn,jy-1,jr  ,jq+1) &
         + coeff(4)*table(nn,jy-1,jr+1,jq+1)

     q1 = coeff(1)*table(nn,jy,jr  ,jq  ) &
         +coeff(2)*table(nn,jy,jr+1,jq  ) &
         +coeff(3)*table(nn,jy,jr  ,jq+1) &
         +coeff(4)*table(nn,jy,jr+1,jq+1)

     q2 = coeff(1)*table(nn,jy+1,jr  ,jq  ) &
         +coeff(2)*table(nn,jy+1,jr+1,jq  ) &
         +coeff(3)*table(nn,jy+1,jr  ,jq+1) &
         +coeff(4)*table(nn,jy+1,jr+1,jq+1)

!    f(nn)  = linter(ye,yl1,q1,yl2,q2)
     f(nn) = q0*yefac1+ &
             q1*yefac2+ &
             q2*yefac3

     !f(nn)  = qinter(ye,yl0,q0,yl1,q1,yl2,q2)

   else if (nn .ne. 20) then

     q0 = coeff(5)*table(nn,jy-1,jr  ,jq-1) &
          + coeff(6)*table(nn,jy-1,jr-1,jq  ) &
          + coeff(7)*table(nn,jy-1,jr  ,jq  ) &
          + coeff(8)*table(nn,jy-1,jr+1,jq  ) &
          + coeff(9)*table(nn,jy-1,jr  ,jq+1) &
          + coeff(10)*table(nn,jy-1,jr+1,jq+1)
     q1 = coeff(5)*table(nn,jy,jr,jq-1) &
          + coeff(6)*table(nn,jy,jr-1,jq) &
          + coeff(7)*table(nn,jy,jr,jq) &
          + coeff(8)*table(nn,jy,jr+1,jq) &
          + coeff(9)*table(nn,jy,jr,jq+1) &
          + coeff(10)*table(nn,jy,jr+1,jq+1)
     q2 = coeff(5)*table(nn,jy+1,jr,jq-1) &
          + coeff(6)*table(nn,jy+1,jr-1,jq) &
          + coeff(7)*table(nn,jy+1,jr,jq) &
          + coeff(8)*table(nn,jy+1,jr+1,jq) &
          + coeff(9)*table(nn,jy+1,jr,jq+1) &
          + coeff(10)*table(nn,jy+1,jr+1,jq+1)

     f(nn) = q0*yefac1+ &
             q1*yefac2+ &
             q2*yefac3
     !f(nn)  = qinter(ye,yl0,q0,yl1,q1,yl2,q2)
     if (nn.eq.5.or.nn.eq.6.or.nn.eq.7.or.nn.eq.8) &
          f(nn)  = min(max(f(nn),0.d0),1.d0)

   end if

  end do

  f(2) = exp(f(2))

  jy_out = jy
  jq_out = jq
  jr_out = jr

contains

  double precision function qinter(x,x1,y1,x2,y2,x3,y3)
    implicit none
    double precision, intent(in) ::  x,x1,y1,x2,y2,x3,y3
    qinter = y1*(x-x2)*(x-x3)/(x1-x2)/(x1-x3)+ &
         y2*(x-x1)*(x-x3)/(x2-x1)/(x2-x3)+ &
         y3*(x-x1)*(x-x2)/(x3-x1)/(x3-x2)
  end function qinter
  
  double precision function linter(x,x0,y0,x1,y1)
    implicit none
    double precision, intent(in) :: x0,y0,x1,y1,x
    linter = y0+(y1-y0)/(x1-x0)*(x-x0)
  end function linter

end subroutine findthis_quad

end module table_module
