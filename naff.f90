module naffModule

  implicit none
  real(8), parameter :: pi    = 3.141592653589793d0
  real(8), parameter :: twopi = 6.283185307179586d0
  complex(8), parameter :: i1 = (0d0,1d0)
  interface fft
    module procedure fft_r, fft_c
  end interface
  private :: fft_recursive, quadratic_interpolation
contains

subroutine naff(tune,amplitude,y,n_mode,x,n,winID,tol)
!=========================================================
! Laskars NAFF algorithm with minor improvement
! input
!   x(n) : length n signal
!   winID(optional)
!   : 0 = rectangular window (=no window)
!     1 = hanning window (1+cos(x))
!     2 = (1+cos(x))^2
!   tol(optional) tolderance 1.0e-4/n by default
! output
!   tune : measured tune
!   amplitude : measured amplitude
!   y : subtracted data
!       y = x - Sum[amplitude[n]*exp(i*2pi*tune[n]*T),{n,1,n_mode}]
!=========================================================
  implicit none
  integer,   intent(in) :: n_mode,n
  real(8),   intent(out):: tune(n_mode)
  complex(8),intent(out):: amplitude(n_mode)
  complex(8),intent(in) :: x(n)
  complex(8),intent(out):: y(n)
  integer,optional,intent(in) :: winID
  real(8),optional,intent(in) :: tol
  complex(8) :: ffty(n)
  integer :: mode, j, iMax, closeMode
  real(8) :: closeTune
  y(:) = x(:)
  if(present(winID)) then
    do mode=1,n_mode
      ffty = fft(y)
      iMax = maxloc(abs(ffty),dim=1)-1
      tune(mode) = iMax/dble(N)
      call quadratic_interpolation(tune(mode),amplitude(mode),y,winID,tol)
      y = y - amplitude(mode)*exp(i1*twopi*tune(mode)*[(j,j=1,N)])
    enddo
  else
    ffty = fft(y)
    iMax = maxloc(abs(ffty),dim=1)-1
    tune(1) = iMax/dble(N)
    call quadratic_interpolation(tune(1),amplitude(1),y,1,tol)
    y = y - amplitude(1)*exp(i1*twopi*tune(1)*[(j,j=1,N)])
    do mode=2,n_mode
      ffty = fft(y)
      iMax = maxloc(abs(ffty),dim=mode)-1
      tune(mode) = iMax/dble(N)
      call quadratic_interpolation(tune(mode),amplitude(mode),y,1,tol)
      closeMode = minloc(abs(tune(mode)-tune(:mode-1)),dim=1)
      closeTune = abs(tune(mode)-tune(closeMode))
      if(closeTune >= 6.5/N) then
        y = y + amplitude(closeMode)*exp(i1*twopi*tune(closeMode)*[(j,j=1,N)])
        call quadratic_interpolation(tune(closeMode),amplitude(closeMode),y,4,tol)
        y = y - amplitude(closeMode)*exp(i1*twopi*tune(closeMode)*[(j,j=1,N)])
        call quadratic_interpolation(tune(mode),amplitude(mode),y,4,tol)
      elseif(closeTune >= 3.2/N) then
        y = y + amplitude(closeMode)*exp(i1*twopi*tune(closeMode)*[(j,j=1,N)])
        call quadratic_interpolation(tune(closeMode),amplitude(closeMode),y,-1,tol)
        y = y - amplitude(closeMode)*exp(i1*twopi*tune(closeMode)*[(j,j=1,N)])
        call quadratic_interpolation(tune(mode),amplitude(mode),y,-1,tol)        
      elseif(closeTune < 1.95/N) then
        y = y + amplitude(closeMode)*exp(i1*twopi*tune(closeMode)*[(j,j=1,N)])
        call quadratic_interpolation(tune(closeMode),amplitude(closeMode),y,0,tol)
        y = y - amplitude(closeMode)*exp(i1*twopi*tune(closeMode)*[(j,j=1,N)])
        call quadratic_interpolation(tune(mode),amplitude(mode),y,0,tol)
      endif
      y = y - amplitude(mode)*exp(i1*twopi*tune(mode)*[(j,j=1,N)])
    enddo
  endif
end subroutine naff

subroutine quadratic_interpolation(tune,amplitude,x,winID,tol)
  implicit none
  complex(8),intent(in) :: x(:)
  complex(8),intent(out):: amplitude
  real(8), intent(inout):: tune
  integer,optional,intent(in) :: winID
  real(8),optional,intent(in) :: tol
  
  integer :: i,j,N,winID1
  integer,parameter :: maxiter = 50
  real(8) :: window(size(x)), tuneStep, tuneRef, amp, tol1
  real(8) :: deltaSqure23,deltaSqure31,deltaSqure12,delta23,delta31,delta12
  real(8) :: TuneScan(3), AmplitudeScan(3)
  complex(8) :: dummy(size(x))
  
  N = size(x)
  if(present(winID)) then
    winID1=winID
  else
    winID1=1
  endif
  if(present(tol)) then
    tol1=tol
  else
    tol1=1d-4/dble(N)
  endif
  ! selection of windows
  select case(winID1)
    case (1:4)    ! (1+cos)^p window 
      window = (1d0+cos(pi*(-1d0+2d0/(N+1d0)*[(i,i=1,N)])))**winID1
    case (-1)     ! blackman-nuttall window
      window = twopi/dble(N-1)*[(i,i=0,N-1)]
      window = 0.3635819-0.4891775*cos(    window)&
                        +0.1365995*cos(2d0*window)&
                        -0.0106411*cos(3d0*window)
    case default  ! rectangular window : best for very close tune separation
      window = 1d0
  end select
  window = window/sum(window)
  
  ! ====  if FFT fails... ====
  TuneScan(1) = tune-1.0/dble(N)
  TuneScan(2) = tune
  TuneScan(3) = tune+1.0/dble(N)
  do i=1,3
    dummy = exp(-i1*twopi*(TuneScan(i))*[(j,j=1,N)])
    AmplitudeScan(i) = -abs(sum(dummy*x*window))
  enddo
  i = minloc(AmplitudeScan,1)
  tune = TuneScan(i)
  ! ==========================
  
  tuneStep = 1.6/dble(N)
  TuneScan(1) = tune-tuneStep
  TuneScan(2) = tune
  TuneScan(3) = tune+tuneStep
  do i=1,3
    dummy = exp(-i1*twopi*TuneScan(i)*[(j,j=1,N)])
    AmplitudeScan(i) = -abs(sum(dummy*x*window))
  enddo
  tune = TuneScan(2)+ ((AmplitudeScan(1)-AmplitudeScan(3))*tuneStep)&
                     /(2*(AmplitudeScan(1)-2*AmplitudeScan(2)+AmplitudeScan(3)))
  dummy = exp(-i1*twopi*tune*[(j,j=1,N)])
  amplitude = sum(dummy*x*window)
  amp = -abs(amplitude)
  if (TuneScan(1) < tune .and. tune < TuneScan(2)) then
      if (amp <= AmplitudeScan(2)) then
          TuneScan(3)=TuneScan(2)
          AmplitudeScan(3)=AmplitudeScan(2)
          TuneScan(2)=tune
          AmplitudeScan(2)=amp
      else
          TuneScan(1)=tune
          AmplitudeScan(1)=amp
      endif
  elseif (TuneScan(2) < tune .and. tune < TuneScan(3)) then
      if (amp <= AmplitudeScan(2)) then
          TuneScan(1)=TuneScan(2)
          AmplitudeScan(1)=AmplitudeScan(2)
          TuneScan(2)=tune
          AmplitudeScan(2)=amp
      else
          TuneScan(3)=tune
          AmplitudeScan(3)=amp
      endif
  endif
  tuneRef = 1d0
  do i=1,maxiter
    if(abs(tune-tuneRef) <= tol1) exit
    tuneRef=tune;
    deltaSqure23 = TuneScan(2)**2-TuneScan(3)**2
    deltaSqure31 = TuneScan(3)**2-TuneScan(1)**2
    deltaSqure12 = TuneScan(1)**2-TuneScan(2)**2
    delta23 = TuneScan(2)-TuneScan(3)
    delta31 = TuneScan(3)-TuneScan(1)
    delta12 = TuneScan(1)-TuneScan(2)
    tune = (deltaSqure23*AmplitudeScan(1)+deltaSqure31*AmplitudeScan(2)+deltaSqure12*AmplitudeScan(3)) &
          /(2d0*(delta23*AmplitudeScan(1)+delta31*AmplitudeScan(2)+delta12*AmplitudeScan(3)))
    dummy = exp(-i1*twopi*tune*[(j,j=1,N)])
    amplitude = sum(dummy*x*Window)
    amp = -abs(amplitude)
    if (TuneScan(1) < tune .and. tune < TuneScan(2)) then
        if (amp <= AmplitudeScan(2)) then
            TuneScan(3)=TuneScan(2)
            AmplitudeScan(3)=AmplitudeScan(2)
            TuneScan(2)=tune
            AmplitudeScan(2)=amp
        else
            TuneScan(1)=tune;
            AmplitudeScan(1)=amp;
        endif
    elseif (TuneScan(2) < tune .and. tune < TuneScan(3)) then
        if (amp <= AmplitudeScan(2)) then
            TuneScan(1)=TuneScan(2)
            AmplitudeScan(1)=AmplitudeScan(2)
            TuneScan(2)=tune
            AmplitudeScan(2)=amp
        else
            TuneScan(3)=tune
            AmplitudeScan(3)=amp
        endif
    else
        tune=TuneScan(2)
    endif
  enddo
end subroutine quadratic_interpolation

recursive subroutine fft_recursive(x)
  complex(8), intent(inout)  :: x(:)
  integer    :: N
  integer    :: i
  complex(8) :: t
  complex(8), dimension(:), allocatable    :: even, odd

  N=size(x)

  if(N .le. 1) return

  allocate(odd((N+1)/2))
  allocate(even(N/2))

  ! divide
  odd =x(1:N:2)
  even=x(2:N:2)

  ! conquer
  call fft_recursive(odd)
  call fft_recursive(even)

  ! combine
  do i=1,N/2
     t=exp(-i1*2d0*pi*(i-1)/dble(N))*even(i)
     x(i)     = odd(i) + t
     x(i+N/2) = odd(i) - t
  end do

  deallocate(odd)
  deallocate(even)
end subroutine fft_recursive

function fft_c(x)
  complex(8), intent(in) :: x(:)
  complex(8) :: fft_c(size(x))
  
  fft_c(:) = x(:)
  call fft_recursive(fft_c)
  return
end function fft_c

function fft_r(x)
  real(8), intent(in) :: x(:)
  complex(8) :: fft_r(size(x))
  
  fft_r(:) = x(:)
  call fft_recursive(fft_r)
  return
end function fft_r


end module 
