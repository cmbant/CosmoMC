    module MpiUtils
    implicit none
    contains

    function GetMpiRank()
    integer GetMpiRank
#ifdef MPI 
    integer ierror
    call mpi_comm_rank(mpi_comm_world,GetMPIrank,ierror)
#else
    GetMpiRank=0
#endif    

    end function GetMpiRank

    function IsMainMPI()
    logical IsMainMPI

    IsMainMPI =  GetMpiRank() == 0

    end function IsMainMPI

    subroutine MpiStop(Msg)
    character(LEN=*), intent(in), optional :: Msg
    integer i
#ifdef MPI 
    integer ierror, MpiRank
#endif

    if (present(Msg)) write(*,*) trim(Msg)

#ifdef MPI
    call mpi_comm_rank(mpi_comm_world,MPIrank,ierror)
    write (*,*) 'MpiStop: ', MpiRank
    call MPI_ABORT(MPI_COMM_WORLD,i)
#endif
    i=1     !put breakpoint on this line to debug
    stop

    end subroutine MpiStop

    subroutine MpiStat(MpiID, MpiSize)
    implicit none
    integer MpiID,MpiSize  
#ifdef MPI  
    integer ierror
    call mpi_comm_rank(mpi_comm_world,MpiID,ierror)
    if (ierror/=MPI_SUCCESS) stop 'MpiStat: MPI rank'
    call mpi_comm_size(mpi_comm_world,MpiSize,ierror)
#else
    MpiID=0
    MpiSize=1   
#endif
    end subroutine MpiStat

    subroutine MpiQuietWait
    !Set MPI thread to sleep, e.g. so can run openmp on cpu instead
#ifdef MPI  
    integer flag, ierr, STATUS(MPI_STATUS_SIZE)
    integer i, MpiId, MpiSize

    call MpiStat(MpiID, MpiSize)
    if (MpiID/=0) then  
        do
            call MPI_IPROBE(0,0,MPI_COMM_WORLD,flag, MPI_STATUS_IGNORE,ierr)
            if (flag/=0) then
                call MPI_RECV(i,1,MPI_INTEGER, 0,0,MPI_COMM_WORLD,status,ierr)
                exit
            end if
            call sleep(1)
        end do
    end if 
#endif
    end subroutine

    subroutine MpiWakeQuietWait
#ifdef MPI  
    integer j, MpiId, MpiSize, ierr,r

    call MpiStat(MpiID, MpiSize)
    if (MpiID==0) then
        do j=1, MpiSize-1              
            call MPI_ISSEND(MpiId,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,r,ierr)
        end do  
    end if
#endif
    end subroutine MpiWakeQuietWait

    end module MpiUtils