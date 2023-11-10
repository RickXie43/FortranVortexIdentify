!============================================================
!Version:
!  Date            Author              Description
!  ========        ========            ===================
!  09/02/23        Rick                Initialize
!  01/03/23        Rick                receive inputdir and output dir from command line
!  21/04/23        Rick                add calculating flitervalue in the first flitervalueframe
!  12/08/23        Rick                enlarge outputdir, realoutputdir to length 150
!
!============================================================
!
! Usage:        #! sh ./make_MPI_CMD.sh
!               #! mpirun -np processer_number ./vortex_mpi_cmd inputdir outputdir

Program mainfindvortex
        Use find_vortex
        Use mpi_f08

        Implicit None
        !MPI Parameters
        Integer::ierror, rank, num_proc, i, j

        !Compute Paremeters
        Character(len=150)::inputdir
        Character(len=150)::outputdir
        Character(len=150)::realoutputdir
        Integer::fileamount = 0
        Character(len=10)::dataname
        Real::fliterthreshold
        Integer, Parameter::flitervalueframe = 30
        Real, Dimension(flitervalueframe)::fliterthresholdlist

        !count the number fileamount of the original data and create directory "Vortexdata_Name" in the outputdir
        !set the input and output directory

        !Read Commend Argument through String 1 to Argu Argu_Name
        Call GET_COMMAND_ARGUMENT(1, inputdir)
        Call GET_COMMAND_ARGUMENT(2, outputdir)
        !inputdir = '/home/rick/Desktop/170509_5V_3.23K_A_3.33s'            !The directory contains original data (No / in the end)
        !outputdir = '/home/rick/Desktop'           !The directory to output result (No / in the end)
        realoutputdir = trim(outputdir)//'/VortexData_'//inputdir(index(inputdir, '/', .true.) + 1:len_trim(inputdir))

        !Use MPI Parallel to compute data
        Call MPI_Init(ierror)
        Call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierror)
        Call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierror)
        rank0print: If (rank == 0) Then
                Write (*, *) 'MPI Initialize Successful.'
                Write (*, *) 'The number of MPI processes is ', num_proc

                call system('mkdir '//realoutputdir)
                call system('(find '//trim(inputdir)//' -type f | wc -l) >'//trim(realoutputdir)//'/fileamount')
                open (unit=16, file=trim(realoutputdir)//'/fileamount', status='old', iostat=ierror)
                read (16, *) fileamount
                close (unit=16)
                call system('rm -f '//trim(realoutputdir)//'/fileamount')

                !send fileamount to other processes
                sendfileamount: Do i = 1, num_proc - 1
                        Call MPI_Send(fileamount, 1, MPI_INT, i, 1, MPI_COMM_WORLD, ierror)
                End Do sendfileamount
        Else
                !receive fileamount from process rank 0
                Call MPI_Recv(fileamount, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_Status_Ignore, ierror)
        End If rank0print

        !calculate flitervalue====================================================================
        calculateonrank0: If (rank == 0) Then
                Write (*, *) "Step 1: Calculating flitervalue ..."
                averagethreshold: Do i = 1, flitervalueframe
                        Write (*, '(a,I2,a,I2,a)') '(', i, '/', flitervalueframe, ')'
                        Call find_flitervalue(trim(inputdir)//'/'//dataname(i), fliterthresholdlist(i))
                End Do averagethreshold
                fliterthreshold = 0
                totalthres: Do i = 1, flitervalueframe
                        fliterthreshold = fliterthreshold + fliterthresholdlist(i)
                End Do totalthres
                fliterthreshold = fliterthreshold/flitervalueframe
                Write (*, *) "the flitervalue calculated by first", flitervalueframe, "frame is ", fliterthreshold
        End If calculateonrank0

        sendrank0flitervalue: If (rank == 0) Then
                !send fliterthreshold to other processes
                sendfliterthreshold: Do i = 1, num_proc - 1
                        Call MPI_Send(fliterthreshold, 1, MPI_REAL, i, 1, MPI_COMM_WORLD, ierror)
                End Do sendfliterthreshold
                Write (*, *) "Step 2: Identifying vortex ..."
        Else
                !receive fliterthreshold from process rank 0
                Call MPI_Recv(fliterthreshold, 1, MPI_REAL, 0, 1, MPI_COMM_WORLD, MPI_Status_Ignore, ierror)
        End If sendrank0flitervalue

        !Paralleling process =====================================================================
        rankfindvortex: Do i = rank + 1, fileamount, num_proc
                Call findvortex(trim(inputdir)//'/'//dataname(i), &
                        &outputfilename=trim(realoutputdir)//'/Vortex_'//dataname(i),&
                        &flitervalue=fliterthreshold)
                Write (*, *) dataname(i), ' Complete!'
        End Do rankfindvortex

        !MPI Finalize
        Call MPI_Finalize(ierror)

End Program mainfindvortex

!convert interger to dataname
Character(len=*) Function dataname(time)
        Implicit None
        Integer, Intent(In)::time
        Character(len=5)::middlestring
        Write (middlestring, 100) time
100     Format(I5.5)
        dataname = 'B'//middlestring//'.txt'
End Function dataname
