program amuse_worker_program
  
  use Mikkola
  use StoppingConditions
  use iso_c_binding
  implicit none
  
  integer HEADER_FLAGS, HEADER_CALL_ID, HEADER_FUNCTION_ID, HEADER_CALL_COUNT, & 
        HEADER_INTEGER_COUNT, HEADER_LONG_COUNT, HEADER_FLOAT_COUNT, & 
        HEADER_DOUBLE_COUNT, HEADER_BOOLEAN_COUNT, HEADER_STRING_COUNT, & 
        HEADER_SIZE, MAX_COMMUNICATORS

  parameter (HEADER_FLAGS=1, HEADER_CALL_ID=2, HEADER_FUNCTION_ID=3, & 
        HEADER_CALL_COUNT=4, HEADER_INTEGER_COUNT=5, HEADER_LONG_COUNT=6, & 
        HEADER_FLOAT_COUNT=7, HEADER_DOUBLE_COUNT=8, & 
        HEADER_BOOLEAN_COUNT=9, HEADER_STRING_COUNT=10, & 
        HEADER_SIZE=11, MAX_COMMUNICATORS = 2048)

  logical NEEDS_MPI
  parameter (NEEDS_MPI=.true.)
  
  
  integer MAX_INTEGERS_IN, MAX_INTEGERS_OUT, MAX_LONGS_IN, MAX_LONGS_OUT, &
  MAX_FLOATS_IN, MAX_FLOATS_OUT, MAX_DOUBLES_IN,MAX_DOUBLES_OUT, &
  MAX_BOOLEANS_IN,MAX_BOOLEANS_OUT, MAX_STRINGS_IN, MAX_STRINGS_OUT
  
  parameter (MAX_INTEGERS_IN=2)
  parameter (MAX_INTEGERS_OUT=3)
  parameter (MAX_LONGS_IN=0)
  parameter (MAX_LONGS_OUT=0)
  parameter (MAX_FLOATS_IN=0)
  parameter (MAX_FLOATS_OUT=0)
  parameter (MAX_DOUBLES_IN=8)
  parameter (MAX_DOUBLES_OUT=8)
  parameter (MAX_BOOLEANS_IN=1)
  parameter (MAX_BOOLEANS_OUT=1)
  parameter (MAX_STRINGS_IN=2)
  parameter (MAX_STRINGS_OUT=1)
  
  
  integer, save :: polling_interval = 0
  integer, save :: last_communicator_id = 0
  integer, save  :: communicators(MAX_COMMUNICATORS)
  integer, save  :: id_to_activate = -1
  integer, save  :: active_communicator_id = -1

  
  integer (c_int32_t), target :: header_in(HEADER_SIZE)
  integer (c_int32_t), target :: header_out(HEADER_SIZE)
  
  integer (c_int32_t), allocatable, target :: integers_in(:)
  integer (c_int32_t), allocatable, target :: integers_out(:)
  
  integer (c_int64_t), allocatable, target :: longs_in(:)
  integer (c_int64_t), allocatable, target :: longs_out(:)
  
  real (c_float), allocatable, target :: floats_in(:)
  real (c_float), allocatable, target :: floats_out(:)
  
  real (c_double), allocatable, target :: doubles_in(:)
  real (c_double), allocatable, target :: doubles_out(:)
  
  logical (c_bool), allocatable, target :: c_booleans_in(:)
  logical (c_bool), allocatable, target :: c_booleans_out(:)

  logical, allocatable, target :: booleans_in(:)
  logical, allocatable, target :: booleans_out(:)
  
  integer (c_int32_t), allocatable, target :: string_sizes_in(:)
  integer (c_int32_t), allocatable, target :: string_sizes_out(:)

  character (c_char), allocatable, target :: strings_in(:) * 256
  character (c_char), allocatable, target :: strings_out(:) * 256

  character (len=1000000) :: characters_in
  character (len=1000000) :: characters_out
  
  character (kind=c_char), target :: c_characters_in(1000000)
  character (kind=c_char), target :: c_characters_out(1000000)

  
  
  integer :: count
  logical :: use_mpi
  character(len=32) :: use_mpi_string
 
  count = command_argument_count()
  
  use_mpi = NEEDS_MPI

  if (count .eq. 0) then
    call run_loop_mpi()
  else if (count .eq. 3) then
    call get_command_argument(3, use_mpi_string)
      
    if (use_mpi_string .eq. 'true') then
      use_mpi = .true.
    else if (use_mpi_string .eq. 'false') then
      use_mpi = .false.
    else
      print*, 'fortran worker: need either true or false as mpi enable arguments, not', use_mpi_string
      stop
    end if
  
    if (use_mpi) then
      call run_loop_sockets_mpi()
    else 
      call run_loop_sockets()
    end if
  else
    print*, 'fortran worker: need either 0 or 3 arguments, not', count
    stop
  end if

  
  CONTAINS
    FUNCTION internal__get_message_polling_interval(outval)
        INTEGER,intent(out) :: outval
        INTEGER :: internal__get_message_polling_interval
        outval = polling_interval
        internal__get_message_polling_interval = 0
    END FUNCTION
    FUNCTION internal__set_message_polling_interval(inval)
        INTEGER,intent(in) :: inval
        INTEGER :: internal__set_message_polling_interval
        polling_interval = inval
        internal__set_message_polling_interval = 0
    END FUNCTION

FUNCTION internal__open_port(outval)
    USE mpi
    IMPLICIT NONE
    character(len=MPI_MAX_PORT_NAME+1), intent(out) :: outval
    INTEGER :: internal__open_port
    INTEGER :: ierror
    call MPI_Open_port(MPI_INFO_NULL, outval, ierror);
    internal__open_port = 0
END FUNCTION

FUNCTION internal__accept_on_port(port_identifier, comm_identifier)
    USE mpi
    IMPLICIT NONE
    character(len=*), intent(in) :: port_identifier
    INTEGER, intent(out) :: comm_identifier
    INTEGER :: internal__accept_on_port
    INTEGER :: ierror, rank
    INTEGER :: mcommunicator, communicator
    last_communicator_id = last_communicator_id + 1
    IF (last_communicator_id .GE. MAX_COMMUNICATORS) THEN
        last_communicator_id = last_communicator_id - 1
        comm_identifier = -1
        internal__accept_on_port = -1
        return;
    END IF
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror);
    IF (rank .EQ. 0) THEN
        call MPI_Comm_accept(port_identifier, MPI_INFO_NULL, 0,  MPI_COMM_SELF, communicator, ierror)
        call MPI_Intercomm_merge(communicator, .FALSE., mcommunicator, ierror)
        call MPI_Intercomm_create(MPI_COMM_WORLD, 0, mcommunicator, 1, 65, communicators(last_communicator_id), ierror)
        call MPI_Comm_free(mcommunicator, ierror)
        call MPI_Comm_free(communicator, ierror)
    ELSE
        call MPI_Intercomm_create(MPI_COMM_WORLD,0, MPI_COMM_NULL, 1, 65, communicators(last_communicator_id), ierror)
    END IF
    comm_identifier = last_communicator_id;
    
    internal__accept_on_port = 0
END FUNCTION

FUNCTION internal__connect_to_port(port_identifier, comm_identifier)
    USE MPI
    IMPLICIT NONE
    character(len=*), intent(in) :: port_identifier
    INTEGER, intent(out) :: comm_identifier
    INTEGER :: internal__connect_to_port
    INTEGER :: ierror, rank
    INTEGER :: mcommunicator, communicator
    last_communicator_id = last_communicator_id + 1
    IF (last_communicator_id .GE. MAX_COMMUNICATORS) THEN
        last_communicator_id = last_communicator_id - 1
        comm_identifier = -1
        internal__connect_to_port = -1
        return;
    END IF
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror);
    
    IF (rank .EQ. 0) THEN
        call MPI_Comm_connect(port_identifier, MPI_INFO_NULL, 0,  MPI_COMM_SELF, communicator, ierror)
        call MPI_Intercomm_merge(communicator, .TRUE., mcommunicator, ierror)
        call MPI_Intercomm_create(MPI_COMM_WORLD, 0, mcommunicator, 0, 65, communicators(last_communicator_id), ierror)
        call MPI_Comm_free(mcommunicator, ierror)
        call MPI_Comm_free(communicator, ierror)
    ELSE
        call MPI_Intercomm_create(MPI_COMM_WORLD,0, MPI_COMM_NULL, 1, 65, communicators(last_communicator_id), ierror)
    END IF
    comm_identifier = last_communicator_id;
    
    internal__connect_to_port = 0
END FUNCTION


FUNCTION  internal__activate_communicator(comm_identifier)
    USE mpi
    IMPLICIT NONE
    INTEGER, intent(in) :: comm_identifier
    INTEGER :: internal__activate_communicator
    
    if ((comm_identifier .LT. 0) .OR. (comm_identifier .GT. last_communicator_id)) then
        internal__activate_communicator = -1
        return 
    end if
    internal__activate_communicator = 0
    id_to_activate = comm_identifier
END FUNCTION

FUNCTION internal__become_code(number_of_workers, modulename, classname)
    IMPLICIT NONE
    character(len=*), intent(in) :: modulename, classname
    integer, intent(in) :: number_of_workers
    INTEGER :: internal__become_code
    
    internal__become_code = 0
END FUNCTION


    SUBROUTINE mpi_recv_header(parent, ioerror)
        use iso_c_binding
        use mpi
        implicit none
              
        integer,intent(in) :: parent
        integer,intent(inout) :: ioerror
        integer :: request_status(MPI_STATUS_SIZE),header_request
        logical is_finished
        
        INTERFACE
          INTEGER (C_INT) FUNCTION usleep(useconds) bind(C)
          !SUBROUTINE usleep(useconds) bind(C)
            use iso_c_binding
            implicit none
            INTEGER(c_int32_t), value  :: useconds
          END
        END INTERFACE
        
        call MPI_Irecv(header_in, HEADER_SIZE, MPI_INTEGER, 0, 989, parent, header_request, ioerror)
        if(polling_interval.GT.0) then
            is_finished = .false.
            call MPI_Test(header_request, is_finished, request_status, ioerror)
            DO WHILE(.NOT. is_finished)
                ioerror =  usleep(int(polling_interval, c_int32_t))
                call MPI_Test(header_request, is_finished, request_status, ioerror)
            END DO
            call MPI_Wait(header_request, request_status, ioerror)
        else
            call MPI_Wait(header_request, request_status, ioerror)
        endif
    END SUBROUTINE

    SUBROUTINE run_loop_mpi
      use mpi
      implicit none
            
      integer :: provided
      integer :: rank, parent, ioerror, max_call_count = 255
      integer :: must_run_loop, maximum_size, total_string_length
      integer i, offset, call_count
      
      call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, ioerror)
      
      ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
      ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
      ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
      ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
      ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
      ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
      ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
      ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
      ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
      ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
      
      ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
      !ensure there is at least one string to return an error code in
      ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      
      call MPI_COMM_GET_PARENT(parent, ioerror)
      call MPI_COMM_RANK(parent, rank, ioerror)
      last_communicator_id = last_communicator_id + 1
      communicators(1) = parent
      active_communicator_id = 1
      
      must_run_loop = 1
      
      do while (must_run_loop .eq. 1)
        if ((id_to_activate .GE. 0) .AND. (id_to_activate .NE. active_communicator_id)) then
            active_communicator_id = id_to_activate
            id_to_activate = -1
            parent = communicators(active_communicator_id)
            call MPI_COMM_RANK(parent, rank, ioerror)
        end if
      
        call mpi_recv_header(parent, ioerror)
        
        !print*, 'fortran: got header ', header_in
        
        call_count = header_in(HEADER_CALL_COUNT)
        
        IF (call_count .gt. max_call_count) THEN
          max_call_count = call_count + 255;
          DEALLOCATE(integers_in)
          DEALLOCATE(integers_out)
          DEALLOCATE(longs_in)
          DEALLOCATE(longs_out)
          DEALLOCATE(floats_in)
          DEALLOCATE(floats_out)
          DEALLOCATE(doubles_in)
          DEALLOCATE(doubles_out)
          DEALLOCATE(c_booleans_in)
          DEALLOCATE(c_booleans_out)
          DEALLOCATE(booleans_in)
          DEALLOCATE(booleans_out)
          DEALLOCATE(string_sizes_in)
          DEALLOCATE(string_sizes_out)
          DEALLOCATE(strings_in)
          DEALLOCATE(strings_out)
          ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
          ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
          ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
          ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
          ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
          ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
          ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
          ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
          ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
          ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
        END IF
    
        if (header_in(HEADER_INTEGER_COUNT) .gt. 0) then
          call MPI_BCast(integers_in, header_in(HEADER_INTEGER_COUNT), MPI_INTEGER, 0, parent, ioError);
        end if
        if (header_in(HEADER_LONG_COUNT) .gt. 0) then
          call MPI_BCast(longs_in, header_in(HEADER_LONG_COUNT), MPI_INTEGER8, 0, parent, ioError);
        end if
        if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
          call MPI_BCast(floats_in, header_in(HEADER_FLOAT_COUNT), MPI_REAL, 0, parent, ioError);
        end if
        if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
          call MPI_BCast(doubles_in, header_in(HEADER_DOUBLE_COUNT), MPI_REAL8, 0, parent, ioError);
        end if
        if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
          ! some older MPI do not define MPI_C_BOOL; this seems to work ok
          ! maybe booleans_in in this call should be replaced by char (more portable) or logical*1  
          call MPI_BCast(c_booleans_in, header_in(HEADER_BOOLEAN_COUNT), MPI_BYTE, 0, parent, ioError);
          do i=1,header_in(HEADER_BOOLEAN_COUNT)
              booleans_in(i)=logical(c_booleans_in(i))
          enddo
        end if
        if (header_in(HEADER_STRING_COUNT) .gt. 0) then
          strings_in = ' '
          call MPI_BCast(string_sizes_in, header_in(HEADER_STRING_COUNT), MPI_INTEGER, 0, parent, ioError);

          maximum_size = 0
          total_string_length = 0
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              total_string_length = total_string_length + string_sizes_in(i) + 1
              if (string_sizes_in(i) .gt. maximum_size) then
                maximum_size = string_sizes_in(i)
              end if
          end do

          if(maximum_size.GT.256) then
            print*, "fortran_worker reports too large string"
            stop          
          endif

          if(total_string_length.GT.1000000) then
            print*, "fortran_worker reports too large string message"
            stop
          endif
          
          call MPI_BCast(characters_in, total_string_length, MPI_CHARACTER, 0, parent, ioError);
          
          offset = 1
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              strings_in(i) = ' '
              strings_in(i)  = characters_in(offset : (offset + string_sizes_in(i)))
              strings_in(i)((string_sizes_in(i) + 1):(string_sizes_in(i) + 1)) = ' ' 
              offset = offset + string_sizes_in(i) + 1
              !print*, 'fortran: strings_in(i) ', i, strings_in(i) , ' of length ', string_sizes_in(i), &
              !' actually of size ', len_trim(strings_in(i))
          end do
          
        end if
        
        header_out = 0
        header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
        header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
        header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
        
        strings_out = ' '
        
        must_run_loop = handle_call()
        
        !print*, 'fortran: sending header ', header_out
    
        if (rank .eq. 0 ) then
    
          call MPI_SEND(header_out, HEADER_SIZE, MPI_INTEGER, 0, 999, parent, ioerror);
    
          if (header_out(HEADER_INTEGER_COUNT) .gt. 0) then
            call MPI_SEND(integers_out,  header_out(HEADER_INTEGER_COUNT), MPI_INTEGER, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_LONG_COUNT) .gt. 0) then
            call MPI_SEND(longs_out,  header_out(HEADER_LONG_COUNT), MPI_INTEGER8, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_FLOAT_COUNT) .gt. 0) then
            call MPI_SEND(floats_out,  header_out(HEADER_FLOAT_COUNT), MPI_REAL, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_DOUBLE_COUNT) .gt. 0) then
            call MPI_SEND(doubles_out, header_out(HEADER_DOUBLE_COUNT), MPI_REAL8, 0, 999, parent, ioerror)
          end if
          if (header_out(HEADER_BOOLEAN_COUNT) .gt. 0) then
            do i=1,header_out(HEADER_BOOLEAN_COUNT)
              c_booleans_out(i)=booleans_out(i)
            enddo
            call MPI_SEND(c_booleans_out, header_out(HEADER_BOOLEAN_COUNT), MPI_BYTE, 0, 999, parent, ioerror)
          end if
       
          if (header_out(HEADER_STRING_COUNT) .gt. 0) then
          
            offset = 1
            do i = 1, header_out(HEADER_STRING_COUNT),1
              
              string_sizes_out(i) = len_trim(strings_out(i))
              
              !print*, 'fortran: sending strings, strings_out(i) ', i, strings_out(i) , ' of length ', string_sizes_out(i), &
              !' actually of size ', len_trim(strings_out(i))
              
              characters_out(offset:offset+string_sizes_out(i)) = strings_out(i)
              offset = offset + string_sizes_out(i) + 1
              characters_out(offset-1:offset-1) = char(0)
            end do

          total_string_length=offset-1
          if(total_string_length.GT.1000000) then
            print*, "fortran_worker reports too large string message"
            stop
          endif
            
            call MPI_SEND(string_sizes_out, header_out(HEADER_STRING_COUNT), MPI_INTEGER, 0, 999, parent, ioerror)
            call MPI_SEND(characters_out, offset -1, MPI_CHARACTER, 0, 999, parent, ioerror)
          end if
        end if
      end do
    
      DEALLOCATE(integers_in)
      DEALLOCATE(integers_out)
      DEALLOCATE(longs_in)
      DEALLOCATE(longs_out)
      DEALLOCATE(floats_in)
      DEALLOCATE(floats_out)
      DEALLOCATE(doubles_in)
      DEALLOCATE(doubles_out)
      DEALLOCATE(booleans_in)
      DEALLOCATE(booleans_out)
      DEALLOCATE(string_sizes_in)
      DEALLOCATE(string_sizes_out)
      DEALLOCATE(strings_in)
      DEALLOCATE(strings_out)
      
      do i = 1, last_communicator_id, 1
            call MPI_COMM_DISCONNECT(communicators(i), ioerror);
      end do
      call MPI_FINALIZE(ioerror)
      return
    end subroutine

  
    SUBROUTINE run_loop_sockets
      use iso_c_binding
      use FortranSocketsInterface
      
      implicit none
    
      integer :: max_call_count = 255
      integer :: must_run_loop, maximum_size, total_string_length
      integer :: i, offset, call_count, port
      character(len=32) :: port_string
      character(kind=c_char, len=64) :: host
      logical (c_bool), allocatable, target :: c_booleans_in(:)
      logical (c_bool), allocatable, target :: c_booleans_out(:)
      
      ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
      ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
      ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
      ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
      ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
      ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
      ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
      ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
      ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      
      ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      
      ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
      ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
      
      !ensure there is at least one string to return an error code in
      ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      ALLOCATE(string_sizes_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      
      call get_command_argument(1, port_string)
      call get_command_argument(2, host)

      read (port_string,*) port
      !add a null character to the end of the string so c knows when the string ends
      host = trim(host) // c_null_char


      call forsockets_init(host, port)
      
      must_run_loop = 1
      
      do while (must_run_loop .eq. 1)
        call receive_integers(c_loc(header_in), HEADER_SIZE)
        
        !print*, 'fortran sockets: got header ', header_in
        
        call_count = header_in(HEADER_CALL_COUNT)
        
        IF (call_count .gt. max_call_count) THEN
          max_call_count = call_count + 255;
          DEALLOCATE(integers_in)
          DEALLOCATE(integers_out)
          DEALLOCATE(longs_in)
          DEALLOCATE(longs_out)
          DEALLOCATE(floats_in)
          DEALLOCATE(floats_out)
          DEALLOCATE(doubles_in)
          DEALLOCATE(doubles_out)
          DEALLOCATE(booleans_in)
          DEALLOCATE(booleans_out)
          DEALLOCATE(c_booleans_in)
          DEALLOCATE(c_booleans_out)
          DEALLOCATE(string_sizes_in)
          DEALLOCATE(string_sizes_out)
          DEALLOCATE(strings_in)
          DEALLOCATE(strings_out)
          ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
          ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
          ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
          ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
          ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
          ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
          ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
          ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
          ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
          ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
        END IF
    
        if (header_in(HEADER_INTEGER_COUNT) .gt. 0) then
          call receive_integers(c_loc(integers_in), header_in(HEADER_INTEGER_COUNT))
        end if
        if (header_in(HEADER_LONG_COUNT) .gt. 0) then
          call receive_longs(c_loc(longs_in), header_in(HEADER_LONG_COUNT))
        end if
        if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
          call receive_floats(c_loc(floats_in), header_in(HEADER_FLOAT_COUNT))
        end if
        if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
          call receive_doubles(c_loc(doubles_in), header_in(HEADER_DOUBLE_COUNT))
        end if
        if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
          call receive_booleans(c_loc(c_booleans_in), header_in(HEADER_BOOLEAN_COUNT))
          do i = 1, header_in(HEADER_BOOLEAN_COUNT), 1
              booleans_in(i) = logical(c_booleans_in(i))
          end do
        end if
        if (header_in(HEADER_STRING_COUNT) .gt. 0) then
          strings_in = ' '
          call receive_integers(c_loc(string_sizes_in), header_in(HEADER_STRING_COUNT))

          maximum_size = 0
          total_string_length = 0
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              total_string_length = total_string_length + string_sizes_in(i) + 1
              if (string_sizes_in(i) .gt. maximum_size) then
                maximum_size = string_sizes_in(i)
              end if
          end do

          if(maximum_size.GT.256) then
            print*, "fortran_worker reports too large string"
            stop          
          endif

          if(total_string_length.GT.1000000) then
            print*, "fortran_worker reports too large string message"
            stop
          endif
          
          call receive_string(c_loc(c_characters_in), total_string_length)
          
          ! this trick is necessary on older gfortran compilers (~<4.9)
          ! as c_loc needs character(len=1)
          do i=1, total_string_length
            characters_in(i:i)=c_characters_in(i)
          enddo
          
          offset = 1
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              strings_in(i) = ' '
              strings_in(i)  = characters_in(offset : (offset + string_sizes_in(i)))
              strings_in(i)((string_sizes_in(i) + 1):(string_sizes_in(i) + 1)) = ' ' 
              offset = offset + string_sizes_in(i) + 1
              !print*, 'fortran: strings_in(i) ', i, strings_in(i) , ' of length ', string_sizes_in(i), &
              !' actually of size ', len_trim(strings_in(i))
          end do

        end if
        
        header_out = 0
        header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
        header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
        header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
        
        strings_out = ' '
        
        must_run_loop = handle_call()
        
        !print*, 'fortran: sending header ', header_out
    
        call send_integers(c_loc(header_out), HEADER_SIZE)

        if (header_out(HEADER_INTEGER_COUNT) .gt. 0) then
          call send_integers(c_loc(integers_out), header_out(HEADER_INTEGER_COUNT))
        end if
        if (header_out(HEADER_LONG_COUNT) .gt. 0) then
          call send_longs(c_loc(longs_out), header_out(HEADER_LONG_COUNT))
        end if
        if (header_out(HEADER_FLOAT_COUNT) .gt. 0) then
          call send_floats(c_loc(floats_out), header_out(HEADER_FLOAT_COUNT))
        end if
        if (header_out(HEADER_DOUBLE_COUNT) .gt. 0) then
          call send_doubles(c_loc(doubles_out), header_out(HEADER_DOUBLE_COUNT))
        end if
        if (header_out(HEADER_BOOLEAN_COUNT) .gt. 0) then
          do i = 1, header_out(HEADER_BOOLEAN_COUNT), 1
              c_booleans_out(i) = logical(booleans_out(i), c_bool) 
          end do

          call send_booleans(c_loc(c_booleans_out), header_out(HEADER_BOOLEAN_COUNT))
        end if
   
        if (header_out(HEADER_STRING_COUNT) .gt. 0) then
          offset = 1
          do i = 1, header_out(HEADER_STRING_COUNT),1
              
            string_sizes_out(i) = len_trim(strings_out(i))
              
              !print*, 'fortran: sending strings, strings_out(i) ', i, strings_out(i) , ' of length ', string_sizes_out(i), &
              !' actually of size ', len_trim(strings_out(i))
              
            characters_out(offset:offset+string_sizes_out(i)) = strings_out(i)
            offset = offset + string_sizes_out(i) + 1
            characters_out(offset-1:offset-1) = char(0)
          end do

          total_string_length=offset-1

          if(total_string_length.GT.1000000) then
            print*, "fortran_worker reports too large string message"
            stop
          endif
          
          do i=1, total_string_length
            c_characters_out(i)=characters_out(i:i)
          enddo

          call send_integers(c_loc(string_sizes_out), header_out(HEADER_STRING_COUNT))
          call send_string(c_loc(c_characters_out), offset-1 )
        end if
      end do
    
      DEALLOCATE(integers_in)
      DEALLOCATE(integers_out)
      DEALLOCATE(longs_in)
      DEALLOCATE(longs_out)
      DEALLOCATE(floats_in)
      DEALLOCATE(floats_out)
      DEALLOCATE(doubles_in)
      DEALLOCATE(doubles_out)
      DEALLOCATE(booleans_in)
      DEALLOCATE(booleans_out)
      DEALLOCATE(c_booleans_in)
      DEALLOCATE(c_booleans_out)
      DEALLOCATE(string_sizes_in)
      DEALLOCATE(string_sizes_out)
      DEALLOCATE(strings_in)
      DEALLOCATE(strings_out)

      call forsockets_close()
      return
    end subroutine

  
    SUBROUTINE run_loop_sockets_mpi
      use iso_c_binding
      use FortranSocketsInterface
      use mpi

      implicit none
            
      integer :: provided
      integer :: max_call_count = 255
      integer :: must_run_loop, maximum_size, total_string_length
      integer :: i, offset, call_count, port, rank, ioerror
      character(len=32) :: port_string
      character(kind=c_char, len=64) :: host
      logical (c_bool), allocatable, target :: c_booleans_in(:)
      logical (c_bool), allocatable, target :: c_booleans_out(:)

      
      ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
      ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
      ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
      ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
      ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
      ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
      ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
      ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
      ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
      ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
      
      ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
      ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
      
      !ensure there is at least one string to return an error code in
      ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      ALLOCATE(string_sizes_out(max(1, max_call_count * MAX_STRINGS_OUT)))
      
      call mpi_init_thread(mpi_thread_multiple, provided, ioerror)
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ioerror)

      if (rank .eq. 0) then
        call get_command_argument(1, port_string)
        call get_command_argument(2, host)

        read (port_string,*) port
        !add a null character to the end of the string so c knows when the string ends
        host = trim(host) // c_null_char

        call forsockets_init(host, port)
      end if
      
      must_run_loop = 1
      
      do while (must_run_loop .eq. 1)
        if (rank .eq. 0) then
          call receive_integers(c_loc(header_in), HEADER_SIZE)
        end if
        call MPI_BCast(header_in, HEADER_SIZE , MPI_INTEGER, 0, MPI_COMM_WORLD, ioerror)
        
        !print*, 'fortran sockets mpi: got header ', header_in
        
        call_count = header_in(HEADER_CALL_COUNT)
        
        IF (call_count .gt. max_call_count) THEN
          max_call_count = call_count + 255;
          DEALLOCATE(integers_in)
          DEALLOCATE(integers_out)
          DEALLOCATE(longs_in)
          DEALLOCATE(longs_out)
          DEALLOCATE(floats_in)
          DEALLOCATE(floats_out)
          DEALLOCATE(doubles_in)
          DEALLOCATE(doubles_out)
          DEALLOCATE(booleans_in)
          DEALLOCATE(booleans_out)
          DEALLOCATE(c_booleans_in)
          DEALLOCATE(c_booleans_out)
          DEALLOCATE(string_sizes_in)
          DEALLOCATE(string_sizes_out)
          DEALLOCATE(strings_in)
          DEALLOCATE(strings_out)
          ALLOCATE(integers_in(max_call_count * MAX_INTEGERS_IN))
          ALLOCATE(integers_out(max_call_count * MAX_INTEGERS_OUT))
          ALLOCATE(longs_in(max_call_count * MAX_LONGS_IN))
          ALLOCATE(longs_out(max_call_count * MAX_LONGS_OUT))
          ALLOCATE(floats_in(max_call_count * MAX_FLOATS_IN))
          ALLOCATE(floats_out(max_call_count * MAX_FLOATS_OUT))
          ALLOCATE(doubles_in(max_call_count * MAX_DOUBLES_IN))
          ALLOCATE(doubles_out(max_call_count * MAX_DOUBLES_OUT))
          ALLOCATE(booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(c_booleans_in(max_call_count * MAX_BOOLEANS_IN))
          ALLOCATE(c_booleans_out(max_call_count * MAX_BOOLEANS_OUT))
          ALLOCATE(string_sizes_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(string_sizes_out(max_call_count * MAX_STRINGS_OUT))
          ALLOCATE(strings_in(max_call_count * MAX_STRINGS_IN))
          ALLOCATE(strings_out(max(1, max_call_count * MAX_STRINGS_OUT)))
        END IF
    
        if (header_in(HEADER_INTEGER_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_integers(c_loc(integers_in), header_in(HEADER_INTEGER_COUNT))
          end if
          call MPI_BCast(integers_in, header_in(HEADER_INTEGER_COUNT), MPI_INTEGER, 0, MPI_COMM_WORLD, ioError);
        end if
        
        if (header_in(HEADER_LONG_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_longs(c_loc(longs_in), header_in(HEADER_LONG_COUNT))
          end if
            call MPI_BCast(longs_in, header_in(HEADER_LONG_COUNT), MPI_INTEGER8, 0, MPI_COMM_WORLD, ioError);
        end if
        
        if (header_in(HEADER_FLOAT_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_floats(c_loc(floats_in), header_in(HEADER_FLOAT_COUNT))
          end if
          call MPI_BCast(floats_in,  header_in(HEADER_FLOAT_COUNT), MPI_REAL, 0, MPI_COMM_WORLD, ioerror)
        end if
        
        if (header_in(HEADER_DOUBLE_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_doubles(c_loc(doubles_in), header_in(HEADER_DOUBLE_COUNT))
          end if
          call MPI_BCast(doubles_in, header_in(HEADER_DOUBLE_COUNT), MPI_REAL8, 0, MPI_COMM_WORLD, ioerror)
        end if
        
        if (header_in(HEADER_BOOLEAN_COUNT) .gt. 0) then
          if (rank .eq. 0) then
            call receive_booleans(c_loc(c_booleans_in), header_in(HEADER_BOOLEAN_COUNT))
            do i = 1, header_in(HEADER_BOOLEAN_COUNT), 1
              booleans_in(i) = logical(c_booleans_in(i))
            end do
          end if
          call MPI_BCast(booleans_in, header_in(HEADER_BOOLEAN_COUNT), MPI_LOGICAL, 0, MPI_COMM_WORLD, ioerror)
        end if
        
        if (header_in(HEADER_STRING_COUNT) .gt. 0) then
          strings_in = ' '
                    
          if (rank .eq. 0) then
            call receive_integers(c_loc(string_sizes_in), header_in(HEADER_STRING_COUNT))
          end if
          call MPI_BCast(string_sizes_in, header_in(HEADER_STRING_COUNT), MPI_INTEGER, 0, MPI_COMM_WORLD, ioError);

          maximum_size = 0
          total_string_length = 0
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              total_string_length = total_string_length + string_sizes_in(i) + 1
              if (string_sizes_in(i) .gt. maximum_size) then
                maximum_size = string_sizes_in(i)
              end if
          end do

          if(maximum_size.GT.256) then
            print*, "fortran_worker reports too large string"
            stop          
          endif

          if(total_string_length.GT.1000000) then
            print*, "fortran_worker reports too large string message"
            stop
          endif

          if (rank .eq. 0) then
            call receive_string(c_loc(c_characters_in), total_string_length)
          endif
          
          do i=1, total_string_length
            characters_in(i:i)=c_characters_in(i)
          enddo
          
          call MPI_BCast(characters_in, total_string_length, MPI_CHARACTER, 0, MPI_COMM_WORLD, ioError);
          
          offset = 1
          do i = 1, header_in(HEADER_STRING_COUNT), 1
              strings_in(i) = ' '
              strings_in(i)  = characters_in(offset : (offset + string_sizes_in(i)))
              strings_in(i)((string_sizes_in(i) + 1):(string_sizes_in(i) + 1)) = ' ' 
              offset = offset + string_sizes_in(i) + 1
              !print*, 'fortran: strings_in(i) ', i, strings_in(i) , ' of length ', string_sizes_in(i), &
              !' actually of size ', len_trim(strings_in(i))
          end do

        end if
        
        header_out = 0
        header_out(HEADER_CALL_ID) = header_in(HEADER_CALL_ID)
        header_out(HEADER_FUNCTION_ID) = header_in(HEADER_FUNCTION_ID)
        header_out(HEADER_CALL_COUNT) = header_in(HEADER_CALL_COUNT)
        
        strings_out = ' '
        
        must_run_loop = handle_call()
        
        call MPI_Barrier(MPI_COMM_WORLD, ioerror)
        
        if (rank .eq. 0) then
        
          !print*, 'fortran: sending header ', header_out
    
          call send_integers(c_loc(header_out), HEADER_SIZE)

          if (header_out(HEADER_INTEGER_COUNT) .gt. 0) then
            call send_integers(c_loc(integers_out), header_out(HEADER_INTEGER_COUNT))
          end if
          if (header_out(HEADER_LONG_COUNT) .gt. 0) then
            call send_longs(c_loc(longs_out), header_out(HEADER_LONG_COUNT))
          end if
          if (header_out(HEADER_FLOAT_COUNT) .gt. 0) then
            call send_floats(c_loc(floats_out), header_out(HEADER_FLOAT_COUNT))
          end if
          if (header_out(HEADER_DOUBLE_COUNT) .gt. 0) then
            call send_doubles(c_loc(doubles_out), header_out(HEADER_DOUBLE_COUNT))
          end if
          if (header_out(HEADER_BOOLEAN_COUNT) .gt. 0) then
            do i = 1, header_out(HEADER_BOOLEAN_COUNT), 1
              c_booleans_out(i) = logical(booleans_out(i), c_bool)
              !print*, 'fortran sockets mpi: sending boolean', booleans_out(i) , i, ' send as ', c_booleans_out(i) 
            end do
        
            call send_booleans(c_loc(c_booleans_out), header_out(HEADER_BOOLEAN_COUNT))
          end if
   
          if (header_out(HEADER_STRING_COUNT) .gt. 0) then
            offset = 1
            do i = 1, header_out(HEADER_STRING_COUNT),1
                
              string_sizes_out(i) = len_trim(strings_out(i))
                
                !print*, 'fortran: sending strings, strings_out(i) ', i, strings_out(i) , ' of length ', string_sizes_out(i), &
                !' actually of size ', len_trim(strings_out(i))
                
              characters_out(offset:offset+string_sizes_out(i)) = strings_out(i)
              offset = offset + string_sizes_out(i) + 1
              characters_out(offset-1:offset-1) = char(0)
            end do
  
          total_string_length=offset-1

          if(total_string_length.GT.1000000) then
            print*, "fortran_Worker reports too large string message"
            stop
          endif

          do i=1, total_string_length
            c_characters_out(i)=characters_out(i:i)
          enddo  
  
            call send_integers(c_loc(string_sizes_out), header_out(HEADER_STRING_COUNT))
            call send_string(c_loc(c_characters_out), offset-1 )
         end if
        end if
      end do
    
      DEALLOCATE(integers_in)
      DEALLOCATE(integers_out)
      DEALLOCATE(longs_in)
      DEALLOCATE(longs_out)
      DEALLOCATE(floats_in)
      DEALLOCATE(floats_out)
      DEALLOCATE(doubles_in)
      DEALLOCATE(doubles_out)
      DEALLOCATE(booleans_in)
      DEALLOCATE(booleans_out)
      DEALLOCATE(string_sizes_in)
      DEALLOCATE(string_sizes_out)
      DEALLOCATE(strings_in)
      DEALLOCATE(strings_out)

      if (rank .eq. 0) then
        call forsockets_close()
      end if
      
      call MPI_FINALIZE(ioerror)
      return
    end subroutine

  integer function handle_call()
    
    implicit none
    integer i, call_count
    call_count = header_in(HEADER_CALL_COUNT)
    handle_call = 1
    SELECT CASE (header_in(HEADER_FUNCTION_ID))
      
      CASE(0)
        handle_call = 0
      CASE(20284990)
        do i = 1, call_count, 1
          integers_out(i) = get_mass( &
            integers_in(i) ,&
            doubles_out(i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(20920053)
        integers_out(1) = commit_particles( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(37921492)
        integers_out(1) = set_stopping_condition_timeout_parameter( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(44188957)
        integers_out(1) = get_time( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(102282623)
        do i = 1, call_count, 1
          integers_out(i) = get_children_of_particle( &
            integers_in(i) ,&
            integers_out(( 1 * call_count) + i) ,&
            integers_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 3 * call_count
        
      
      CASE(104547857)
        do i = 1, call_count, 1
          integers_out(i) = set_mass( &
            integers_in(i) ,&
            doubles_in(i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(107212404)
        integers_out(1) = get_stopping_condition_maximum_density_parameter( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(120233917)
        integers_out(1) = set_stopping_condition_out_of_box_use_center_of_mass_parameter( &
          booleans_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(128926247)
        integers_out(1) = get_index_of_first_particle( &
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(146223729)
        integers_out(1) = get_stopping_condition_out_of_box_use_center_of_mass_parameter( &
          booleans_out(1) &
        )
        header_out(HEADER_BOOLEAN_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(159095171)
        do i = 1, call_count, 1
          integers_out(i) = enable_stopping_condition( &
            integers_in(i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(205426934)
        integers_out(1) = get_total_radius( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(205548127)
        integers_out(1) = get_stopping_condition_maximum_internal_energy_parameter( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(219539042)
        do i = 1, call_count, 1
          integers_out(i) = get_number_of_stopping_conditions_set( &
            integers_out(( 1 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(271440754)
        do i = 1, call_count, 1
          integers_out(i) = is_stopping_condition_set( &
            integers_in(i) ,&
            integers_out(( 1 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(290264013)
        do i = 1, call_count, 1
          integers_out(i) = new_particle( &
            integers_out(( 1 * call_count) + i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * call_count) + i) ,&
            doubles_in(( 2 * call_count) + i) ,&
            doubles_in(( 3 * call_count) + i) ,&
            doubles_in(( 4 * call_count) + i) ,&
            doubles_in(( 5 * call_count) + i) ,&
            doubles_in(( 6 * call_count) + i) ,&
            doubles_in(( 7 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(353555706)
        integers_out(1) = get_radiated_gravitational_energy( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(383112453)
        integers_out(1) = internal__connect_to_port( &
          strings_in(1) ,&
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(384567015)
        integers_out(1) = get_total_mass( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(408512489)
        integers_out(1) = get_lightspeed( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(416721977)
        integers_out(1) = evolve_model( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(508605261)
        integers_out(1) = set_stopping_condition_out_of_box_parameter( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(542058817)
        integers_out(1) = set_eps2( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(632979349)
        integers_out(1) = set_stopping_condition_number_of_steps_parameter( &
          integers_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(639951605)
        integers_out(1) = get_stopping_condition_timeout_parameter( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(653682513)
        integers_out(1) = get_begin_time( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(658631024)
        integers_out(1) = get_eps2( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(662348285)
        integers_out(1) = get_stopping_condition_minimum_internal_energy_parameter( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(678380482)
        integers_out(1) = get_index_of_next_particle( &
          integers_in(1) ,&
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(727361823)
        integers_out(1) = internal__set_message_polling_interval( &
          integers_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(728786188)
        do i = 1, call_count, 1
          integers_out(i) = delete_particle( &
            integers_in(i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(733749514)
        do i = 1, call_count, 1
          integers_out(i) = is_stopping_condition_enabled( &
            integers_in(i) ,&
            integers_out(( 1 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(835969050)
        do i = 1, call_count, 1
          integers_out(i) = get_potential( &
            integers_in(i) ,&
            doubles_out(i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(859923342)
        do i = 1, call_count, 1
          integers_out(i) = get_id_of_added_particle( &
            integers_in(i) ,&
            integers_out(( 1 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(887125873)
        integers_out(1) = synchronize_model( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(919768251)
        integers_out(1) = internal__get_message_polling_interval( &
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(929181341)
        do i = 1, call_count, 1
          integers_out(i) = set_state( &
            integers_in(i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * call_count) + i) ,&
            doubles_in(( 2 * call_count) + i) ,&
            doubles_in(( 3 * call_count) + i) ,&
            doubles_in(( 4 * call_count) + i) ,&
            doubles_in(( 5 * call_count) + i) ,&
            doubles_in(( 6 * call_count) + i) ,&
            doubles_in(( 7 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(953558780)
        integers_out(1) = get_stopping_condition_minimum_density_parameter( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(967950880)
        do i = 1, call_count, 1
          integers_out(i) = get_state( &
            integers_in(i) ,&
            doubles_out(i) ,&
            doubles_out(( 1 * call_count) + i) ,&
            doubles_out(( 2 * call_count) + i) ,&
            doubles_out(( 3 * call_count) + i) ,&
            doubles_out(( 4 * call_count) + i) ,&
            doubles_out(( 5 * call_count) + i) ,&
            doubles_out(( 6 * call_count) + i) ,&
            doubles_out(( 7 * call_count) + i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 8 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1024680297)
        integers_out(1) = get_time_step( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1050085724)
        integers_out(1) = recommit_particles( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1060979074)
        integers_out(1) = set_tolerance( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1071152125)
        integers_out(1) = get_kinetic_energy( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1082457792)
        integers_out(1) = get_number_of_particles( &
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(1098838617)
        integers_out(1) = get_stopping_condition_number_of_steps_parameter( &
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(1137702459)
        do i = 1, call_count, 1
          integers_out(i) = disable_stopping_condition( &
            integers_in(i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1221190952)
        do i = 1, call_count, 1
          integers_out(i) = set_acceleration( &
            integers_in(i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * call_count) + i) ,&
            doubles_in(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1231060790)
        do i = 1, call_count, 1
          integers_out(i) = get_center_of_mass_position( &
            doubles_out(i) ,&
            doubles_out(( 1 * call_count) + i) ,&
            doubles_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 3 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1236118821)
        integers_out(1) = set_time_step( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1245044653)
        integers_out(1) = get_number_of_particles_added( &
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(1266912718)
        integers_out(1) = get_tolerance( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1315680918)
        do i = 1, call_count, 1
          integers_out(i) = get_center_of_mass_velocity( &
            doubles_out(i) ,&
            doubles_out(( 1 * call_count) + i) ,&
            doubles_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 3 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1317242279)
        do i = 1, call_count, 1
          integers_out(i) = get_radius( &
            integers_in(i) ,&
            doubles_out(i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1393481104)
        integers_out(1) = set_stopping_condition_minimum_internal_energy_parameter( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1436937242)
        integers_out(1) = set_maximum_number_of_particles( &
          integers_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1508430886)
        integers_out(1) = set_begin_time( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1513820140)
        integers_out(1) = internal__become_code( &
          integers_in(1) ,&
          strings_in(1) ,&
          strings_in(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1544727346)
        integers_out(1) = set_stopping_condition_minimum_density_parameter( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1605138627)
        integers_out(1) = get_total_energy( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1623630901)
        do i = 1, call_count, 1
          integers_out(i) = set_radius( &
            integers_in(i) ,&
            doubles_in(i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1625942852)
        do i = 1, call_count, 1
          integers_out(i) = has_stopping_condition( &
            integers_in(i) ,&
            integers_out(( 1 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(1644113439)
        integers_out(1) = cleanup_code( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1655137210)
        integers_out(1) = set_stopping_condition_maximum_density_parameter( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1706735752)
        integers_out(1) = internal__activate_communicator( &
          integers_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1732760736)
        integers_out(1) = set_lightspeed( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1744145122)
        integers_out(1) = recommit_parameters( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1768994498)
        integers_out(1) = initialize_code( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1802157412)
        integers_out(1) = set_evolve_to_exact_time( &
          booleans_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1829880540)
        integers_out(1) = internal__open_port( &
          strings_out(1) &
        )
        header_out(HEADER_STRING_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1852958273)
        integers_out(1) = get_potential_energy( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1891926464)
        integers_out(1) = internal__accept_on_port( &
          strings_in(1) ,&
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(1892689129)
        do i = 1, call_count, 1
          integers_out(i) = get_velocity( &
            integers_in(i) ,&
            doubles_out(i) ,&
            doubles_out(( 1 * call_count) + i) ,&
            doubles_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 3 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(1937183958)
        integers_out(1) = get_stopping_condition_out_of_box_parameter( &
          doubles_out(1) &
        )
        header_out(HEADER_DOUBLE_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2010900811)
        do i = 1, call_count, 1
          integers_out(i) = get_position( &
            integers_in(i) ,&
            doubles_out(i) ,&
            doubles_out(( 1 * call_count) + i) ,&
            doubles_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 3 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2016681522)
        integers_out(1) = set_stopping_condition_maximum_internal_energy_parameter( &
          doubles_in(1) &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2026192840)
        do i = 1, call_count, 1
          integers_out(i) = set_position( &
            integers_in(i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * call_count) + i) ,&
            doubles_in(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2046699524)
        do i = 1, call_count, 1
          integers_out(i) = get_stopping_condition_info( &
            integers_in(i) ,&
            integers_out(( 1 * call_count) + i) ,&
            integers_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 3 * call_count
        
      
      CASE(2061473599)
        do i = 1, call_count, 1
          integers_out(i) = get_acceleration( &
            integers_in(i) ,&
            doubles_out(i) ,&
            doubles_out(( 1 * call_count) + i) ,&
            doubles_out(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_DOUBLE_COUNT) = 3 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2069478464)
        integers_out(1) = commit_parameters( &
        )
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2097634310)
        integers_out(1) = get_evolve_to_exact_time( &
          booleans_out(1) &
        )
        header_out(HEADER_BOOLEAN_COUNT) = 1 * call_count
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE(2129795713)
        do i = 1, call_count, 1
          integers_out(i) = get_stopping_condition_particle_index( &
            integers_in(i) ,&
            integers_in(( 1 * call_count) + i) ,&
            integers_out(( 1 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(2142977458)
        integers_out(1) = get_maximum_number_of_particles( &
          integers_out(2) &
        )
        header_out(HEADER_INTEGER_COUNT) = 2 * call_count
        
      
      CASE(2144268908)
        do i = 1, call_count, 1
          integers_out(i) = set_velocity( &
            integers_in(i) ,&
            doubles_in(i) ,&
            doubles_in(( 1 * call_count) + i) ,&
            doubles_in(( 2 * call_count) + i) &
          )
        end do
        header_out(HEADER_INTEGER_COUNT) = 1 * call_count
        
      
      CASE DEFAULT
        header_out(HEADER_STRING_COUNT) = 1
        header_out(HEADER_FLAGS) = IOR(header_out(HEADER_FLAGS), 256) 
        strings_out(1) = 'error, illegal function id'
    END SELECT
    return
  end function
end program amuse_worker_program