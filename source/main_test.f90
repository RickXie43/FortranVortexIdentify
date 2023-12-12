
!============================================================
!Version:
!  Date            Author              Description
!  ========        ========            ===================
!  04/04/23        Rick                Initialize
!
!============================================================
!Test Program Using Subroutine findvortex

Program testfindvortex
        Use find_vortex
        Implicit None
        Real::value1 = 0
        Call find_flitervalue('./data/B00001.txt', value1)
        Write (*, *) value1
        Call findvortex('./data/B00001.txt',&
                &outputfilename='vortexdata.txt',&
                &flitervalue=value1)
End Program testfindvortex
