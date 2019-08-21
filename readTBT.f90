subroutine get_TBTsize(fID,nturn,npt)
  integer, intent(in) :: fID
  integer, intent(out):: nturn,npt
  
  logical :: file_open 
  integer :: iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 346192 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  nturn = 0
  loop1 : DO
    READ(iUnit,iostat=eastat) npt
    IF (eastat < 0) THEN
      EXIT loop1
    ELSE IF (eastat > 0) THEN
      STOP 'IO-error'
    ENDIF
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    nturn = nturn+1
  END DO loop1
  close(iUnit)
end subroutine get_TBTsize

subroutine get_TBTsize_npt(fID,nturn,npt)
  integer, intent(in) :: fID, nturn
  integer, intent(out):: npt
  
  logical :: file_open 
  integer :: i,iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 469242 
  file_open = .true. 
  loop0 : do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do loop0
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  loop1 : DO i=1,nturn-1
    READ(iUnit,iostat=eastat) npt
    IF (eastat > 0) STOP 'IO-error'    
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
  END DO loop1
  READ(iUnit,iostat=eastat) npt
  close(iUnit)
end subroutine get_TBTsize_npt


subroutine get_TBTdata(fID,nturn,npt,ke,mass,freq,pIndex,pData)
  integer, intent(in) :: fID,nturn,npt
  integer*8, intent(out) :: pIndex(npt)
  double precision, intent(in)  :: ke, mass, freq
  double precision, intent(out) :: pData(nturn,6,npt)
  logical :: file_open 
  integer :: iUnit,eastat,mpt
  double precision :: x_norm, px_norm, gamma, beta
  integer*8, allocatable :: pIndexTmp(:)
  double precision, allocatable :: pDataTmp(:,:)
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  iUnit = 465192 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  DO i=1,nturn-1
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
  END DO
  
  READ(iUnit,iostat=eastat)
  READ(iUnit,iostat=eastat) pIndex
  READ(iUnit,iostat=eastat) pData(nturn,:,:)
  
  rewind(iUnit)
  READ(iUnit,iostat=eastat) mpt
  allocate(pIndexTmp(mpt))
  READ(iUnit,iostat=eastat) pIndexTmp(1:mpt)
  allocate(pDataTmp(6,mpt))
  READ(iUnit,iostat=eastat) pDataTmp(1:6,1:mpt)
  k=1
  do j=1,mpt
    if(pIndex(k)==pIndexTmp(j)) then
      pData(1,:,k) = pDataTmp(:,j)
      k=k+1
    endif
  enddo
  
  do i=2,nturn-1
    READ(iUnit,iostat=eastat) mpt
    READ(iUnit,iostat=eastat) pIndexTmp(1:mpt)
    READ(iUnit,iostat=eastat) pDataTmp(1:6,1:mpt)
    k=1
    do j=1,mpt
      if(pIndex(k)==pIndexTmp(j)) then
        pData(i,:,k) = pDataTmp(:,j)
        k=k+1
      endif
    enddo
  enddo
  
  close(iUnit)

  gamma = ke/mass+1.0
  beta = sqrt(1.0-1.0/(gamma*gamma))
  x_norm = 2*freq*3.141592653589793/299792458
  px_norm = gamma*beta
  
  pData(:,1,:) = pData(:,1,:)/x_norm
  pData(:,2,:) = pData(:,2,:)/px_norm
  pData(:,3,:) = pData(:,3,:)/x_norm
  pData(:,4,:) = pData(:,4,:)/px_norm
  pData(:,5,:) = pData(:,5,:)*(180.0/3.14159265359)
  pData(:,6,:) = pData(:,6,:)*mass
  
end subroutine get_TBTdata


subroutine get_TBTsize_integral(fID,nturn,npt)
  integer, intent(in) :: fID
  integer, intent(out):: nturn,npt
  
  logical :: file_open 
  integer :: iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 541692 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.integral.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  nturn = 0
  loop1 : DO
    READ(iUnit,iostat=eastat) npt
    IF (eastat < 0) THEN
      EXIT loop1
    ELSE IF (eastat > 0) THEN
      STOP 'IO-error'
    ENDIF
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    nturn = nturn+1
  END DO loop1
  close(iUnit)
end subroutine get_TBTsize_integral

subroutine get_TBTsize_npt_integral(fID,nturn,npt)
  integer, intent(in) :: fID,nturn
  integer, intent(out):: npt
  
  logical :: file_open 
  integer :: i,iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 834692 
  file_open = .true. 
  loop0 : do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do loop0
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.integral.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  loop1 : DO i=1,nturn-1
    READ(iUnit,iostat=eastat)
    IF (eastat > 0) STOP 'IO-error'
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
  END DO loop1
  READ(iUnit,iostat=eastat) npt
  close(iUnit)
end subroutine get_TBTsize_npt_integral



subroutine get_TBTdata_integral(fID,nturn,npt,pIndex,Integral)
  integer, intent(in) :: fID,nturn,npt
  integer*8, intent(out) :: pIndex(npt)
  double precision, intent(out) :: Integral(nturn,2,npt)
  logical :: file_open 
  integer :: iUnit,eastat,mpt
  double precision :: x_norm, px_norm, gamma, beta
  integer*8, allocatable :: pIndexTmp(:)
  double precision, allocatable :: intTmp(:,:)
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  iUnit = 964692 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.integral.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  DO i=1,nturn-1
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
  END DO
  
  READ(iUnit,iostat=eastat)
  READ(iUnit,iostat=eastat) pIndex
  READ(iUnit,iostat=eastat) Integral(nturn,:,:)
  
  rewind(iUnit)
  READ(iUnit,iostat=eastat) mpt
  allocate(pIndexTmp(mpt))
  READ(iUnit,iostat=eastat) pIndexTmp(1:mpt)
  allocate(intTmp(2,mpt))
  READ(iUnit,iostat=eastat) intTmp(1:2,1:mpt)
  k=1
  do j=1,mpt
    if(pIndex(k)==pIndexTmp(j)) then
      Integral(1,:,k) = intTmp(:,j)
      k=k+1
    endif
  enddo
  
  do i=2,nturn-1
    READ(iUnit,iostat=eastat) mpt
    READ(iUnit,iostat=eastat) pIndexTmp(1:mpt)
    READ(iUnit,iostat=eastat) intTmp(1:2,1:mpt)
    k=1
    do j=1,mpt
      if(pIndex(k)==pIndexTmp(j)) then
        Integral(i,:,k) = intTmp(:,j)
        k=k+1
      endif
    enddo
  enddo
  close(iUnit)
end subroutine get_TBTdata_integral


subroutine get_TBTsize_integral_onMomentum(fID,nturn,npt)
  integer, intent(in) :: fID
  integer, intent(out):: nturn,npt
  
  logical :: file_open 
  integer :: iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 734692 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.integral.onMomentum.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  nturn = 0
  loop1 : DO
    READ(iUnit,iostat=eastat) npt
    IF (eastat < 0) THEN
      EXIT loop1
    ELSE IF (eastat > 0) THEN
      STOP 'IO-error'
    ENDIF
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    nturn = nturn+1
  END DO loop1
  close(iUnit)
end subroutine get_TBTsize_integral_onMomentum

subroutine get_TBTsize_npt_integral_onMomentum(fID,nturn,npt)
  integer, intent(in) :: fID,nturn
  integer, intent(out):: npt
  
  logical :: file_open 
  integer :: i,iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 744692 
  file_open = .true. 
  loop0 : do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do loop0
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.integral.onMomentum.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  loop1 : DO i=1,nturn-1
    READ(iUnit,iostat=eastat)
    IF (eastat > 0) STOP 'IO-error'
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
  END DO loop1
  READ(iUnit,iostat=eastat) npt
  close(iUnit)
end subroutine get_TBTsize_npt_integral_onMomentum



subroutine get_TBTdata_integral_onMomentum(fID,nturn,npt,pIndex,Integral)
  integer, intent(in) :: fID,nturn,npt
  integer*8, intent(out) :: pIndex(npt)
  double precision, intent(out) :: Integral(nturn,2,npt)
  logical :: file_open 
  integer :: iUnit,eastat,mpt
  double precision :: x_norm, px_norm, gamma, beta
  integer*8, allocatable :: pIndexTmp(:)
  double precision, allocatable :: intTmp(:,:)
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  iUnit = 974692 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.integral.onMomentum.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  DO i=1,nturn-1
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
  END DO
  
  READ(iUnit,iostat=eastat)
  READ(iUnit,iostat=eastat) pIndex
  READ(iUnit,iostat=eastat) Integral(nturn,:,:)
  
  rewind(iUnit)
  READ(iUnit,iostat=eastat) mpt
  allocate(pIndexTmp(mpt))
  READ(iUnit,iostat=eastat) pIndexTmp(1:mpt)
  allocate(intTmp(2,mpt))
  READ(iUnit,iostat=eastat) intTmp(1:2,1:mpt)
  k=1
  do j=1,mpt
    if(pIndex(k)==pIndexTmp(j)) then
      Integral(1,:,k) = intTmp(:,j)
      k=k+1
    endif
  enddo
  
  do i=2,nturn-1
    READ(iUnit,iostat=eastat) mpt
    READ(iUnit,iostat=eastat) pIndexTmp(1:mpt)
    READ(iUnit,iostat=eastat) intTmp(1:2,1:mpt)
    k=1
    do j=1,mpt
      if(pIndex(k)==pIndexTmp(j)) then
        Integral(i,:,k) = intTmp(:,j)
        k=k+1
      endif
    enddo
  enddo
  close(iUnit)
end subroutine get_TBTdata_integral_onMomentum



subroutine get_rawTBTsize(fID,nturn,npt)
  integer, intent(in) :: fID
  integer, intent(out):: nturn,npt
  
  logical :: file_open 
  integer :: iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 564692 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  READ(iUnit,iostat=eastat) npt
  READ(iUnit,iostat=eastat)
  READ(iUnit,iostat=eastat)
  nturn = 1
  loop1 : DO
    READ(iUnit,iostat=eastat)
    IF (eastat < 0) THEN
      EXIT loop1
    ELSE IF (eastat > 0) THEN
      STOP 'IO-error'
    ENDIF
    READ(iUnit,iostat=eastat)
    READ(iUnit,iostat=eastat)
    nturn = nturn+1
  END DO loop1
  close(iUnit)
end subroutine get_rawTBTsize


subroutine get_rawTBTsize_npt(fID,npt)
  integer, intent(in) :: fID
  integer, intent(out):: npt
  
  logical :: file_open 
  integer :: i,iUnit,eastat
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 754692 
  file_open = .true. 
  loop0 : do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do loop0
  
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')

  READ(iUnit,iostat=eastat) npt
  close(iUnit)
end subroutine get_rawTBTsize_npt


subroutine get_rawTBTdata(fID,nturn,npt,ke,mass,freq,pIndex,pData)
!================================================================
! test purpose (check if working correctly)
! get all raw data including pIndex every turn
! not only alive (at the end) particle data
!================================================================
  integer, intent(in) :: fID,nturn,npt
  integer*8, intent(out) :: pIndex(nturn,npt)
  double precision, intent(in)  :: ke, mass, freq
  double precision, intent(out) :: pData(nturn,6,npt)
  logical :: file_open 
  integer :: iUnit,eastat,mpt
  double precision :: x_norm, px_norm, gamma, beta
  double precision :: pDataTmp(6,npt)
  character(len=6) :: num2str
  character(len=8), parameter :: fmt_ = "(I0)"
  
  pData = 0d0
  iUnit = 624692 
  file_open = .true. 
  do while ( file_open ) 
    iUnit = iUnit + 1 
    inquire(iUnit, opened = file_open ) 
  end do 
  write(num2str,fmt_) fID
  OPEN(iUnit, file='TBT.'//trim(num2str), status='old', action='read', form='unformatted', position='rewind')
  DO i=1,nturn
    READ(iUnit,iostat=eastat) mpt
    READ(iUnit,iostat=eastat) pIndex(i,1:mpt)
    READ(iUnit,iostat=eastat) pDataTmp(1:6,1:mpt)
    pData(i,1:6,1:mpt) = pDataTmp(1:6,1:mpt)
  END DO
  close(iUnit)

  gamma = ke/mass+1.0
  beta = sqrt(1.0-1.0/(gamma*gamma))
  x_norm = 2*freq*3.141592653589793/299792458
  px_norm = gamma*beta
  
  pData(:,1,:) = pData(:,1,:)/x_norm
  pData(:,2,:) = pData(:,2,:)/px_norm
  pData(:,3,:) = pData(:,3,:)/x_norm
  pData(:,4,:) = pData(:,4,:)/px_norm
  pData(:,5,:) = pData(:,5,:)*(180.0/3.14159265359)
  pData(:,6,:) = pData(:,6,:)*mass
  
end subroutine get_rawTBTdata