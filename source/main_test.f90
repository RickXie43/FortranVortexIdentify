
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
        Character(len=50)::testdataname
        testdataname='./data/B00001.txt'
        Call find_flitervalue(trim(testdataname), value1)
        Write (*, *) value1
        Call findvortex(trim(testdataname),&
                &outputfilename='vortexdata.txt',&
                &flitervalue=value1)
End Program testfindvortex
