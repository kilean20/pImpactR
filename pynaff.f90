subroutine pynaff(tune,amplitude,y,n_mode,x,n,window_id,tol)
  use naffmodule
  implicit none
  integer,   intent(in) :: n_mode, n
  real(8),   intent(out):: tune(n_mode)
  complex(8),intent(out):: amplitude(n_mode)
  complex(8),intent(out):: y(n)
  complex(8),intent(in) :: x(n)
  integer,optional,intent(in) :: window_id
  real(8),optional,intent(in) :: tol
  integer :: p1,maxiter1
  real(8) :: tol1
  !========== !!!do not remove!!! == f2py compiler directive ============
  !f2py integer :: window_id = 1
  !f2py real*8  :: tol = 1.0d-3/n
  !======================================================================
  p1 = 1
  maxiter1 = 50
  tol1 = 1.0d-4/n
  if(present(window_id)) p1=window_id
  if(present(tol)) tol1=tol

  if(p1==1 .and. tol1 == 1.0d-3/n) then
      call naff(tune,amplitude,y,n_mode,x,n)
  else
    call naff(tune,amplitude,y,n_mode,x,n,winID=p1,tol=tol1)
  endif
  
end subroutine pynaff
