program runNaff
  use naffModule
  implicit none
  
  integer :: n
  integer, parameter :: n_mode = 2, T(1024) = (/ (n, n = 1, 1024) /)
  real(8), parameter :: pi    = 3.141592653589793d0
  real(8), parameter :: twopi = 6.283185307179586d0
  complex(8), parameter :: i1 = (0d0,1d0)
  
  real(8)    :: tune(n_mode)
  complex(8) :: data_in(1024), data_out(1024), amplitude(n_mode)
  
  
  data_in = i1*1.0*exp(i1*twopi*0.333333333*T) &
          &  +  0.5*exp(i1*twopi*0.888888888*T) &
          &  -  0.1*exp(i1*twopi*0.111111111*T)
  
  call naff(tune,amplitude,data_out,n_mode,data_in,n=1024)
  
  print*, tune
  print*, amplitude
  
  call naff(tune(:1),amplitude(:1),data_out,1,data_out,n=1024)
  print*, tune(:1)
  print*, amplitude(:1)
  
end program