subroutine get_TBTsize(fID,nturn,npt)
  integer, intent(in) :: fID
  integer, intent(out):: nturn,npt
  
  logical :: file_open 
  integer :: iUnit,eastat
  character(len=4) :: num2str
  character(len=6), parameter :: fmt_ = "(I0)"
  
  
  iUnit = 4692 
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



subroutine get_TBTdata(fID,nturn,npt,pData)
  integer, intent(in) :: fID,nturn,npt
  double precision, intent(out) :: pData(nturn,6,npt)
  logical :: file_open 
  integer :: iUnit,eastat,mpt,pIndex(npt)
  integer, allocatable :: pIndexTmp(:)
  double precision, allocatable :: pDataTmp(:,:)
  character(len=4) :: num2str
  character(len=6), parameter :: fmt_ = "(I0)"
  
  iUnit = 4692 
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

end subroutine get_TBTdata
