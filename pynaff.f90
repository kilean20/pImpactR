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
  !f2py integer :: window_id = -10
  !f2py real*8  :: tol = -10.0
  !======================================================================
  p1 = 1
  maxiter1 = 50
  tol1 = 1.0d-4/size(x)
  if(present(window_id)) p1=window_id
  if(present(tol)) tol1=tol

  if(p1==-10) then
    if (tol1==-10.0) then
      call naff(tune,amplitude,y,n_mode,x,n)
    else
      call naff(tune,amplitude,y,n_mode,x,n,tol=tol1)
    endif
  elseif(tol1==-10.0) then
    call naff(tune,amplitude,y,n_mode,x,n,winID=p1)
  else
    call naff(tune,amplitude,y,n_mode,x,n,winID=p1,tol=tol1)
  endif
  
end subroutine pynaff