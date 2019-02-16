subroutine read_phasespace_size(fID,npt)
  integer, intent(in) :: fID
  integer, intent(out):: npt
  integer :: eastat
  
  if(fID>0) error stop 'Error <- get_raw_phasespace_size :: fID must be negative for binary phase-space output'
  OPEN(-fID, status='old', action='read', form='unformatted')
  READ(-fID,iostat=eastat) npt
  !
  !npt = 0
  !loop1 : DO
  !  READ(-fID,iostat=eastat) 
  !  IF (eastat < 0) THEN
  !    EXIT loop1
  !  ELSE IF (eastat > 0) THEN
  !    error STOP 'IO-error'
  !  ENDIF
  !  npt = npt+1
  !END DO loop1
  close(-fID)
end subroutine read_phasespace_size

subroutine read_phasespace(fID,npt,pData)
  integer, intent(in) :: fID,npt
  double precision, intent(out) :: pData(9,npt)
  integer :: eastat
  
  if(fID>0) error stop 'Error <- get_raw_phasespace_size :: fID must be negative for binary phase-space output'
  OPEN(-fID, status='old', action='read', form='unformatted')
  READ(-fID,iostat=eastat) 
  DO i=1,npt
    READ(-fID,iostat=eastat) pData(1:9,i)
  END DO
  close(-fID)
  
end subroutine read_phasespace
