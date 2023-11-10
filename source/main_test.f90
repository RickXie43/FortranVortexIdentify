
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
        Call find_flitervalue('~/Data/PIV_Data_Excel/2.13K_A_G3.8/180720_5V_2.13K_A_4s/B00534.txt', value1)
        Write (*, *) value1
        Call findvortex('~/Data/PIV_Data_Excel/2.13K_A_G3.8/180720_5V_2.13K_A_4s/B00534.txt',&
                &outputfilename='vortexdata.txt',&
                &flitervalue=value1)
End Program testfindvortex
