!============================================================
!Version:
!  Date            Author              Description
!  ========        ========            ===================
!  09/02/23        Rick                Initialize
!
!============================================================
Program mainfindvortex
        Use find_vortex

        Implicit None

        !Compute Paremeters
        Character(len=100)::inputdir
        Character(len=100)::outputdir
        Character(len=100)::realoutputdir
        Integer::fileamount = 0
        Integer::i, ierror
        Character(len=10)::dataname

        !count the number fileamount of the original data and create directory "Vortexdata_Name" in the outputdir
        !set the input and output directory
        inputdir = '/home/server/Desktop/2023_datavortex/PIV_Data_Excel/2.13K_A_G3.8/180720_5V_2.13K_A_4s'            !The directory contains original data (No / in the end)
        outputdir = '/home/server/Desktop'           !The directory to output result (No / in the end)
        realoutputdir = trim(outputdir)//'/VortexData_'//inputdir(index(inputdir, '/', .true.) + 1:len_trim(inputdir))

        call system('mkdir '//realoutputdir)
        call system('(find '//trim(inputdir)//' -type f | wc -l) >'//trim(realoutputdir)//'/fileamount')
        open (unit=16, file=trim(realoutputdir)//'/fileamount', status='old', iostat=ierror)
        read (16, *) fileamount
        close (unit=16)
        call system('rm -f '//trim(realoutputdir)//'/fileamount')

        !process =====================================================================
        rankfindvortex: Do i = 1, fileamount
                Call findvortex(trim(inputdir)//'/'//dataname(i), outputfilename=trim(realoutputdir)//'/Vortex_'//dataname(i))
                Write (*, *) 'Task Completed ', dataname(i), ' ', i, '/', fileamount
        End Do rankfindvortex

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
