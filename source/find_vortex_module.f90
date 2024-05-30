!============================================================
!Version:
!  Date            Author              Description
!  ========        ========            ===================
!  09/02/23        Rick                Initialize
!  11/02/23        Rick                fix bug when input dataname is absolute dataname
!  12/02/23        Rick                debug with initialize variable
!  20/02/23        Rick                set dimension m, n automatically
!  24/02/23        Rick                debug the auto new line at the beginning of vortexdata
!  31/03/23        Rick                modulize writetofile part with subroutine writedata
!  31/03/23        Rick                add nearest vortex distance to type VORTEX
!  04/04/23        Rick                add interpolation function using bicubic spine line
!  07/04/23        Rick                add subroutine q-criterion option instead of liutex-omega
!  08/04/23        Rick                fixbug of insane frame with big big Q and Qc determination
!  18/04/23        Rick                record the max vorticity into type vortex
!  19/04/23        Rick                use subroutine vortexcenter and vortexinfo1 to record to vortexarray
!  19/04/23        Rick                add subroutine in vortexcenter to interpolate accurate_center
!  19/04/23        Rick                add subroutine vortexfliter to drop plume vortex
!  20/04/23        Rick                comments and eff_radius setting in first several lines
!  21/04/23        Rick                add subroutine find_flitervalue to determine flitervalue in main
!  12/08/23        Rick                add parameter fliter_beta to recuce fake vprtex, default should be 1.0
!  09/11/23        Rick                add include macrodefinitions in parameters.F90 to set constants
!  12/12/23        Rick                add subroutine vortexarea to calculate and record the area of a vortex
!  30/05/24        Rick                fix bug, the order of x_array and y_array should be strictly increasing in
!  bicubicinterpolation
!
!============================================================
!Usage:
!       The subroutine findvortex(dataname, outputvortexarray_outside, outputfilename) is used to find vortex
!               The input is the txt dataname, and the output contains two options:
!                       1.Call an array of derive type VORTEX (The definition is below), the length of array should large enough.
!                       2.Write into a file(ASCII), the filename is outputfilename.
!       a. In order to get more accurate data, bicubic interpolation are used on velocity field. Thus, x_min, x_max, y_min, y_max, and
!               interpolationpoints should be set depend on your computing resource before calculation.
!       b. nx, ny in SUBROUTINE bicubicinterpolation should also be set, or there will be a warning of m, n which is the nx, ny.
!
!============================================================

!Use macro-definition to set parameters
#include "parameters.F90"

Module find_vortex
        Implicit None
        Private
        Integer::m, n                                                           !set the dimension of data
        !the data quality of interpolation data
        Real, Parameter::x_min = X_MIN_DEF
        Real, Parameter::x_max = X_MAX_DEF
        Real, Parameter::y_min = Y_MIN_DEF
        Real, Parameter::y_max = Y_MAX_DEF
        Integer, Parameter::interpolationpoints = INTERPOLATIONPOINTS_DEF       !interpolation points in two dimension using bicubic
        Integer, Parameter::nx = NX_DEF                                         !the original point number in x, in interpolation
        Integer, Parameter::ny = NY_DEF                                         !the original point number in y
        Real, Parameter::eff_radius = EFF_RADIUS_DEF                            !the effective radius(Unit:M), data out of range will be drop
        Integer, Parameter::fliter_switch = 1                                   !0 to turn of fliter, 1 to turn on fliter with flitervalue
        Real, Parameter::fliter_beta = FLITER_BETA_DEF                           !multiple a parameter on fliter value to reduce fake vortex.
        !the definition of the element(derive type) of outputvortexarray
        Type::vortex
                Real::xc = 0.
                Real::yc = 0.
                Real::vorticity = 0.
                Real::neardistance = 0.
                Integer::direction = 0
                Integer::time = 0
                Real::maximumvorticity = 0.
                Real::area = 0.
                Type(vortex), Pointer::previous_vortex => Null()
                Type(vortex), Pointer::next_vortex => Null()
                Integer::beginmarker = 1
        End Type vortex
        !Set Public Variable and Routine that canbe used in Program
        Public::vortex, findvortex, find_flitervalue

Contains
        Subroutine findvortex(dataname, outputvortexarray_outside, outputfilename, flitervalue)
                Implicit None
                Character(len=*), Intent(In)::dataname                                          !the original data name "b00001.txt"
                Type(vortex), Intent(Out), Dimension(:), Optional::outputvortexarray_outside    !an array of type vortex
                Type(vortex), Allocatable, Dimension(:)::outputvortexarray                      !an array of type vortex
                Integer::vortex_amount                                                          !the amount of vortexs
                Character(len=*), Intent(In), Optional::outputfilename                          !write result in outputfile if input
                Real, Intent(In), Optional::flitervalue                                         !fliterthreshold in vortexfliter

                !basic physic quantity
                Real, Dimension(:, :), Allocatable::x, y, u, v
                Real, Dimension(:, :), Allocatable::ux, uy, vx, vy                              !ux is partial derive of u with respect to x
                Real, Dimension(:, :), Allocatable::vorticity

                !interpolation quantity
                Real, Dimension(:, :), Allocatable::x_interpolation, y_interpolation, u_interpolation, v_interpolation

                !the groupjudgevalue matrix is the value matrix used to group vortex
                Real, Dimension(:, :), Allocatable::groupjudgevalue
                !the groupvalue matrix tag same value for the point in the same vortex using groupjudgevalue
                Integer, Dimension(:, :), Allocatable::groupvalue

                !set dimension m, n from dataname
                Call setdatamn(dataname, m, n)

                !Read data from dataname to matrix x,y,u,v
                !And Convert the Unit of x, y from mm to m, and delete the data out of experiment range
                Allocate (x(m, n))
                Allocate (y(m, n))
                Allocate (u(m, n))
                Allocate (v(m, n))
                x = 0
                y = 0
                u = 0
                v = 0
                Call readconvertdata(dataname, x, y, u, v)

                !set bicubic interpolation data data_interpolation, and reset m, n with interpolation data
                Allocate (x_interpolation(interpolationpoints, interpolationpoints))
                Allocate (y_interpolation(interpolationpoints, interpolationpoints))
                Allocate (u_interpolation(interpolationpoints, interpolationpoints))
                Allocate (v_interpolation(interpolationpoints, interpolationpoints))
                x_interpolation = 0
                y_interpolation = 0
                u_interpolation = 0
                v_interpolation = 0
                Call bicubicinterpolation(x, y, u, v, m, n,&
                        &x_interpolation, y_interpolation, u_interpolation, v_interpolation)
                Deallocate (x, y, u, v)
                Allocate (x(m, n))
                Allocate (y(m, n))
                Allocate (u(m, n))
                Allocate (v(m, n))
                x = x_interpolation
                y = y_interpolation
                u = u_interpolation
                v = v_interpolation
                Deallocate (x_interpolation, y_interpolation, u_interpolation, v_interpolation)

                !Calculate Basic matrix ux, uy, vx, vy, vorticity
                Allocate (ux(m, n))
                Allocate (uy(m, n))
                Allocate (vx(m, n))
                Allocate (vy(m, n))
                Allocate (vorticity(m, n))
                ux = 0
                uy = 0
                vx = 0
                vy = 0
                vorticity = 0
                Call calculatebasicmatrix(x, y, u, v, ux, uy, vx, vy, vorticity)

                !Calculate groupjudgevalue
                !Using Liutex-Omega Identify Method to findvortex
                Allocate (groupjudgevalue(m, n))
                groupjudgevalue = 0

                Call LiutexOmega(groupjudgevalue, ux, uy, vx, vy, vorticity)
                !                Call qcriterion(groupjudgevalue, ux, uy, vx, vy, vorticity)

                !Group Vortex using groupjudgevalue matrix given by specific vortex identify method
                Allocate (groupvalue(m, n))
                groupvalue = 0
                Call groupallvortex(groupvalue, groupjudgevalue)

                !Allocate outputvortexarray for all the vortexs
                vortex_amount = maxval(groupvalue)
                Allocate (outputvortexarray(vortex_amount))

                !record center of vortex, total vorticity and maximumvorticity
                Call vortexcenter(outputvortexarray, groupvalue, vorticity, x, y)

                !record the area of vortex
                Call vortexarea(outputvortexarray, groupvalue)

                !vortex fliter: drop turbulence plume with maximumvorticity less than most probable maximumvorticity
                ifpresent_flitervalue: If (present(flitervalue)) Then
                        Call vortexfliter_liutexomega(outputvortexarray, fliter_beta*flitervalue)
                End If ifpresent_flitervalue

                !record time, neardistance and direction of vortex
                Call vortexinfo1(outputvortexarray, dataname)

                !count and write down all vortexs which is identify and going to be recorded
                vortex_amount = size(outputvortexarray)
                Write (*, *) 'There are ', vortex_amount, ' vortexs.'

                !if outputvortexarray_outside present, set the vortex data (call subroutine setdatatoarray)
                ifoutputarrayoutside_present: If (Present(outputvortexarray_outside) .And. fliter_switch == 1) Then
                        Call setdatatoarray(outputvortexarray_outside, outputvortexarray)
                End If ifoutputarrayoutside_present

                !if outputfilename present, write the vortex data to file (call subroutine writedata)
                ifoutputfilename_present: If (Present(outputfilename)) Then
                        Call writedata(dataname, outputfilename, outputvortexarray)
                End If ifoutputfilename_present

        End Subroutine findvortex

        Subroutine find_flitervalue(dataname, flitervalue)
                Implicit None
                Character(len=*), Intent(In)::dataname                                          !the original data name "b00001.txt"
                Real, Intent(Out)::flitervalue
                Type(vortex), Allocatable, Dimension(:)::outputvortexarray                      !an array of type vortex
                Integer::vortex_amount                                                          !the amount of vortexs

                !basic physic quantity
                Real, Dimension(:, :), Allocatable::x, y, u, v
                Real, Dimension(:, :), Allocatable::ux, uy, vx, vy                              !ux is partial derive of u with respect to x
                Real, Dimension(:, :), Allocatable::vorticity

                !interpolation quantity
                Real, Dimension(:, :), Allocatable::x_interpolation, y_interpolation, u_interpolation, v_interpolation

                !the groupjudgevalue matrix is the value matrix used to group vortex
                Real, Dimension(:, :), Allocatable::groupjudgevalue
                !the groupvalue matrix tag same value for the point in the same vortex using groupjudgevalue
                Integer, Dimension(:, :), Allocatable::groupvalue

                !set dimension m, n from dataname
                Call setdatamn(dataname, m, n)

                !Read data from dataname to matrix x,y,u,v
                !And Convert the Unit of x, y from mm to m, and delete the data out of experiment range
                Allocate (x(m, n))
                Allocate (y(m, n))
                Allocate (u(m, n))
                Allocate (v(m, n))
                x = 0
                y = 0
                u = 0
                v = 0
                Call readconvertdata(dataname, x, y, u, v)

                !set bicubic interpolation data data_interpolation, and reset m, n with interpolation data
                Allocate (x_interpolation(interpolationpoints, interpolationpoints))
                Allocate (y_interpolation(interpolationpoints, interpolationpoints))
                Allocate (u_interpolation(interpolationpoints, interpolationpoints))
                Allocate (v_interpolation(interpolationpoints, interpolationpoints))
                x_interpolation = 0
                y_interpolation = 0
                u_interpolation = 0
                v_interpolation = 0
                Call bicubicinterpolation(x, y, u, v, m, n,&
                        &x_interpolation, y_interpolation, u_interpolation, v_interpolation)
                Deallocate (x, y, u, v)
                Allocate (x(m, n))
                Allocate (y(m, n))
                Allocate (u(m, n))
                Allocate (v(m, n))
                x = x_interpolation
                y = y_interpolation
                u = u_interpolation
                v = v_interpolation
                Deallocate (x_interpolation, y_interpolation, u_interpolation, v_interpolation)

                !Calculate Basic matrix ux, uy, vx, vy, vorticity
                Allocate (ux(m, n))
                Allocate (uy(m, n))
                Allocate (vx(m, n))
                Allocate (vy(m, n))
                Allocate (vorticity(m, n))
                ux = 0
                uy = 0
                vx = 0
                vy = 0
                vorticity = 0
                Call calculatebasicmatrix(x, y, u, v, ux, uy, vx, vy, vorticity)

                !Calculate groupjudgevalue
                !Using Liutex-Omega Identify Method to findvortex
                Allocate (groupjudgevalue(m, n))
                groupjudgevalue = 0

                Call LiutexOmega(groupjudgevalue, ux, uy, vx, vy, vorticity)

                !Group Vortex using groupjudgevalue matrix given by specific vortex identify method
                Allocate (groupvalue(m, n))
                groupvalue = 0
                Call groupallvortex(groupvalue, groupjudgevalue)

                !Allocate outputvortexarray for all the vortexs
                vortex_amount = maxval(groupvalue)
                Allocate (outputvortexarray(vortex_amount))

                !record center of vortex, total vorticity and maximumvorticity
                Call vortexcenter(outputvortexarray, groupvalue, vorticity, x, y)

                !vortex fliter: drop turbulence plume with maximumvorticity less than most probable maximumvorticity
                Call flitervalue_liutexomega(outputvortexarray, flitervalue)
        End Subroutine find_flitervalue
        !==========================================================================================================

        !the subroutine to set the dimension m, n of dataname
        Subroutine setdatamn(dataname, m, n)
                Implicit None
                Character(len=*), Intent(In)::dataname
                Integer, Intent(Out)::m
                Integer, Intent(Out)::n
                Real::throw, check, checkbefore
                Integer::countm, countn, ierror
                throw = 0.
                check = 0.
                checkbefore = 0.
                countm = 0
                countn = 0
                Open (Unit=99, File=dataname, Status='Old', Iostat=ierror)
                Read (99, *)
                countmnumber: Do
                Read (99, *) throw, check, throw, throw
                endcheckm: If (Abs(checkbefore - check) > 0.01 .And. countm /= 0) Then
                        Exit countmnumber
                End If endcheckm
                countm = countm + 1
                checkbefore = check
                End Do countmnumber
                m = countm
                Do While (.Not. eof(99))
                Read (99, *)
                countn = countn + 1
                End do
                n = (countn + countm + 1)/m
                Close (Unit=99)
        End Subroutine setdatamn

        !the subroutine to read and convert data
        Subroutine readconvertdata(dataname, x, y, u, v)
                Implicit None
                Character(len=*), Intent(In)::dataname
                Real, Dimension(m, n), Intent(Out)::x, y, u, v
                Integer::ierror
                Integer::i, j
                Open (Unit=9, File=dataname, Status='Old', Iostat=ierror)
                Read (9, *)
                readdata_n: Do j = 1, n
                readdata_m: Do i = 1, m
                Read (9, *) x(i, j), y(i, j), u(i, j), v(i, j)
                End Do readdata_m
                End Do readdata_n
                Close (Unit=9)
                x = x/1000
                y = y/1000
                Forall (i=1:m, j=1:n, x(i, j)**2 + y(i, j)**2 > eff_radius**2)
                        u(i, j) = 0.
                        v(i, j) = 0.
                End Forall
        End Subroutine readconvertdata

        !The subroutine to calculate interpolation data of x, y, u ,v(Using bicubicinterpolation method), and reset m, n
        Subroutine bicubicinterpolation(x, y, u, v, m, n,&
                        &x_interpolation, y_interpolation, u_interpolation, v_interpolation)
                Use bspline_kinds_module, only: wp, ip
                Use bspline_module
                Implicit None
                Real, Intent(In), Dimension(:, :)::x, y, u, v
                Integer, Intent(InOut)::m, n
                Real, Intent(Out), Dimension(:, :)::x_interpolation, y_interpolation, u_interpolation, v_interpolation
                Integer::i, j
                Integer::xorder, yorder
                !                Integer(ip), parameter :: nx = 103     !! number of points in x
                !                Integer(ip), parameter :: ny = 86     !! number of points in y
                Integer(ip), parameter :: kx = 4     !! order in x
                Integer(ip), parameter :: ky = 4     !! order in y
                Integer(ip), parameter :: iknot = 0  !! automatically select the knots
                Integer(ip) :: inbvx = 1  !! automatically select the knots
                Integer(ip) :: inbvy = 1  !! automatically select the knots
                Integer(ip) :: iloy = 1  !! automatically select the knots
                Integer(ip)::iflag
                Real(wp) :: x_array(nx), y_array(ny)
                Real(wp) :: u_2d(nx, ny), v_2d(nx, ny)
                Real(wp) :: coefficient_u(nx, ny), coefficient_v(nx, ny)
                Real(wp) :: tx(nx + kx), ty(ny + ky)
                Real(wp)::valuee
                Real(wp), dimension(ky)                       :: w1_2du, w1_2dv
                Real(wp), dimension(3*max(kx, ky))             :: w2_2du, w2_2dv
                Real(wp)            :: xx, yy, uu, vv
                x_interpolation = 0
                y_interpolation = 0
                u_interpolation = 0
                v_interpolation = 0

                !Judge the nx ny setting conform the setting of originaldata
                conformmn: If (nx /= m .Or. ny /= n) Then
                        Write (*, *) 'WARNING: The nx, ny in Subroutine bicubicinterpolation is not equal to m, n!'
                        Write (*, *) 'm = ', m
                        Write (*, *) 'n = ', n
                        Write (*, *) 'PLEASE reset the nx, ny!'
                End If conformmn

                !the x,y array should be strictly increasing, due to interpolation function, test order
                xorderjudge: If(x(2,1)-x(1,1)>0) Then
                        xorder = 1
                Else
                        xorder = -1
                End If xorderjudge
                yorderjudge: If(y(1,2)-y(1,1)>0) Then
                        yorder = 1
                Else
                        yorder = -1
                End If yorderjudge

                !Set x,y of originaldata
                xoriginaldata: If(xorder > 0) Then
                        setx1: Do i = 1, nx
                        x_array(i) = x(i, 1)
                        End Do setx1
                Else
                        setx2: Do i = 1, nx
                        x_array(i) = x(nx - i + 1, 1)
                        End Do setx2
                End If xoriginaldata
                yoriginaldata: If(yorder > 0) Then
                        sety1: Do j = 1, ny
                        y_array(j) = y(1, j)
                        End Do sety1
                Else
                        sety2: Do j = 1, ny
                        y_array(j) = y(1, ny - j + 1)
                        End Do sety2
                End If yoriginaldata

                !Set u,v of orginaldata, same as x,y order
                iforder1: If(xorder>0) Then
                        iforder2: If(yorder>0) Then
                                setuv1: Do i = 1, nx
                                setuvj1: Do j = 1, ny
                                u_2d(i, j) = u(i, j)
                                v_2d(i, j) = v(i, j)
                                End Do setuvj1
                                End Do setuv1
                        Else
                                setuv2: Do i = 1, nx
                                setuvj2: Do j = 1, ny
                                u_2d(i, j) = u(i, ny - j + 1)
                                v_2d(i, j) = v(i, ny - j + 1)
                                End Do setuvj2
                                End Do setuv2
                        End If iforder2
                Else
                        iforder3: If(yorder>0) Then
                                setuv3: Do i = 1, nx
                                setuvj3: Do j = 1, ny
                                u_2d(i, j) = u(nx-i+1, j)
                                v_2d(i, j) = v(nx-i+1, j)
                                End Do setuvj3
                                End Do setuv3
                        Else
                                setuv4: Do i = 1, nx
                                setuvj4: Do j = 1, ny
                                u_2d(i, j) = u(nx-i+1, ny - j + 1)
                                v_2d(i, j) = v(nx-i+1, ny - j + 1)
                                End Do setuvj4
                                End Do setuv4
                        End If iforder3
                End If iforder1


                !calculate interpolation coefficient
                Call db2ink(x_array, nx, y_array, ny, u_2d, kx, ky, iknot, tx, ty, coefficient_u, iflag)
                Call db2ink(x_array, nx, y_array, ny, v_2d, kx, ky, iknot, tx, ty, coefficient_v, iflag)

                !set x, y after interpolation
                setxinter: Do i = 1, interpolationpoints
                setyinter: Do j = 1, interpolationpoints
                x_interpolation(i, j) = x_min + i*(x_max - x_min)/interpolationpoints
                y_interpolation(i, j) = y_min + j*(y_max - y_min)/interpolationpoints
                End Do setyinter
                End Do setxinter

                !set u, v after interpolation
                setvelocityi: Do i = 1, interpolationpoints
                setvelocityj: Do j = 1, interpolationpoints
                xx = x_interpolation(i, j)
                yy = y_interpolation(i, j)
                Call db2val(xx, yy, 0, 0, tx, ty, nx, ny, kx, ky,&
                        &coefficient_u, uu, iflag, inbvx, inbvy, iloy, w1_2du, w2_2du)
                u_interpolation(i, j) = uu
                Call db2val(xx, yy, 0, 0, tx, ty, nx, ny, kx, ky,&
                        &coefficient_v, vv, iflag, inbvx, inbvy, iloy, w1_2dv, w2_2dv)
                v_interpolation(i, j) = vv
                End Do setvelocityj
                End Do setvelocityi

                !reset m, n
                m = interpolationpoints
                n = interpolationpoints
        End Subroutine bicubicinterpolation

        !calculate Basic matrix ux, uy, vx, vy, vorticity
        Subroutine calculatebasicmatrix(x, y, u, v, ux, uy, vx, vy, vorticity)
                Implicit None
                Real, Intent(In), Dimension(m, n)::x, y, u, v
                Real, Intent(Out), Dimension(m, n)::ux, uy, vx, vy, vorticity
                Integer::i, j
                calculate_bm: Do i = 2, m - 1
                calculate_bn: Do j = 2, n - 1
                ux(i, j) = (u(i + 1, j) - u(i - 1, j))/(x(i + 1, j) - x(i - 1, j))
                uy(i, j) = (u(i, j + 1) - u(i, j - 1))/(y(i, j + 1) - y(i, j - 1))
                vx(i, j) = (v(i + 1, j) - v(i - 1, j))/(x(i + 1, j) - x(i - 1, j))
                vy(i, j) = (v(i, j + 1) - v(i, j - 1))/(y(i, j + 1) - y(i, j - 1))
                vorticity(i, j) = vx(i, j) - uy(i, j)
                End Do calculate_bn
                End Do calculate_bm
        End Subroutine calculatebasicmatrix

        !Subroutine of Liutex-Omega identify method, the output is groupjudgevalue
        Subroutine LiutexOmega(groupjudgevalue, ux, uy, vx, vy, vorticity)
                Implicit None
                Real, Intent(Out), Dimension(m, n)::groupjudgevalue
                Real, Intent(In), Dimension(m, n)::ux, uy, vx, vy, vorticity
                Integer::i, j
                !the parameter of the Vortex Identify Method
                Real, Dimension(:, :), Allocatable::alpha, beta

                Real, Parameter::omegarc = 0.7                        !omegarc is the critical point of omegar for Liutex theory
                Real::b = 0.001                                         !the parameter in L-O Method
                Real::epsil = 0
                Real, Dimension(:, :), Allocatable::ba                  !ba is the matrix of beta**2-alpha**2
                Real::maxba = 0                                         !maxba is the maximum of beta**2-alpha**2
                Real, Dimension(:, :), Allocatable::omegar
                Real, Dimension(:, :), Allocatable::omegar_limit        !omegar_limit is the matrix of omegar higher then omegarc

                !Initialize
                Allocate (omegar_limit(m, n))
                Allocate (omegar(m, n))
                Allocate (alpha(m, n))
                Allocate (beta(m, n))
                Allocate (ba(m, n))
                alpha = 0
                beta = 0
                b = 0.001
                epsil = 0
                ba = 0
                maxba = 0
                omegar_limit = 0
                omegar = 0

                !1.calculate matrix: alpha, beta
                calculate_m: Do i = 2, m - 1
                calculate_n: Do j = 2, n - 1
                alpha(i, j) = 0.5*sqrt((vy(i, j) - ux(i, j))**2 + (vx(i, j) + uy(i, j))**2)
                beta(i, j) = 0.5*vorticity(i, j)
                End Do calculate_n
                End Do calculate_m

                !2.calculate parameter epsilon
                ba_m: do i = 2, m - 1
                ba_n: do j = 2, n - 1
                ba(i, j) = beta(i, j)**2 - alpha(i, j)**2
                end do ba_n
                end do ba_m
                maxba = maxval(ba)
                epsil = b*maxba

                !3.calculate omegar and omegar_limit(the matrix used to group vortex)
                omegar_m: Do i = 2, m - 1
                omegar_n: Do j = 2, n - 1
                omegar(i, j) = (beta(i, j))**2/(alpha(i, j)**2 + beta(i, j)**2 + epsil)
                End Do omegar_n
                End Do omegar_m

                omegar_limit = 0
                omegarlimit_m: Do i = 2, m - 1
                omegarlimit_n: Do j = 2, n - 1
                omegacondition: If (omegar(i, j) > omegarc) Then
                        omegar_limit(i, j) = omegar(i, j)
                End If omegacondition
                End Do omegarlimit_n
                End Do omegarlimit_m

                !4.set the matrix used to groupvalue (which is omegar_limit if using Liutex-Omega Identify Method)
                groupjudgevalue = omegar_limit
        End Subroutine LiutexOmega

        !Subroutine of Q-Criterion identify method, the output is groupjudgevalue
        Subroutine qcriterion(groupjudgevalue, ux, uy, vx, vy, vorticity)
                Implicit None
                Real, Intent(Out), Dimension(m, n)::groupjudgevalue
                Real, Intent(In), Dimension(m, n)::ux, uy, vx, vy, vorticity
                Integer::i, j
                !the parameter of the Vortex Identify Method
                Real, Dimension(:, :), Allocatable:: q_matrix
                Real, Dimension(:, :), Allocatable:: q_limit
                Real::qc, maxq, minq, averq, varianceq
                Integer, Parameter::binnumber = 1800
                Integer::bincount
                Integer, Dimension(binnumber)::q_hist

                !Initialize
                Allocate (q_matrix(m, n))
                Allocate (q_limit(m, n))
                q_matrix = 0
                qc = 0
                maxq = 0
                minq = 0
                q_hist = 0

                !1. calculate matrix: q_matrix
                qcalculate_m: Do i = 2, m - 1
                qcalculate_n: Do j = 2, n - 1
                q_matrix(i, j) = (ux(i, j) + vy(i, j))**2 - 4*(ux(i, j)*vy(i, j) - vx(i, j)*uy(i, j))
                End Do qcalculate_n
                End Do qcalculate_m

                !2. calculate averq and varianceq
                averq = 0
                qsscalculate_m: Do i = 2, m - 1
                qsscalculate_n: Do j = 2, n - 1
                averq = averq + q_matrix(i, j)
                End Do qsscalculate_n
                End Do qsscalculate_m
                averq = averq/((m - 2)*(n - 2))
                varianceq = 0
                qscalculate_m: Do i = 2, m - 1
                qscalculate_n: Do j = 2, n - 1
                varianceq = varianceq + (q_matrix(i, j) - averq)**2
                End Do qscalculate_n
                End Do qscalculate_m
                varianceq = sqrt(varianceq/((m - 2)*(n - 2)))

                !2.1 calculate the max bin q value - maxq
                maxqcalculate_m: Do i = 2, m - 1
                maxqcalculate_n: Do j = 2, n - 1
                maxqjudge: If (q_matrix(i, j) > maxq .And. q_matrix(i, j) < (averq + 3*varianceq)) Then
                        maxq = q_matrix(i, j)
                End If maxqjudge
                End Do maxqcalculate_n
                End Do maxqcalculate_m
                minq = averq - 3*varianceq

                !3.calculate critical parameter qc
                bin_m: do i = 2, m - 1
                bin_n: do j = 2, n - 1
                ifnotoutrangedata: If (q_matrix(i, j) < maxq .And. q_matrix(i, j) > minq) Then
                        bincount = ceiling((q_matrix(i, j) - minq)/((maxq - minq)/binnumber))
                        q_hist(bincount) = q_hist(bincount) + 1
                End If ifnotoutrangedata
                end do bin_n
                end do bin_m
                findlocation: Do i = 1, binnumber
                ifismax: If (q_hist(i) == maxval(q_hist)) Then
                        qc = minq + (1.0*i - 0.5)*((maxq - minq)/binnumber)
                End If ifismax
                End Do findlocation
                Write (*, *) qc

                !4.calculate q_limit(the matrix used to group vortex)
                q_limit = 0
                qlimit_m: Do i = 2, m - 1
                qlimit_n: Do j = 2, n - 1
                qcondition: If (q_matrix(i, j) > qc) Then
                        q_limit(i, j) = q_matrix(i, j)
                End If qcondition
                End Do qlimit_n
                End Do qlimit_m

                !4.set the matrix used to groupvalue (which is q_limit if using q-Criterion Identify Method)
                groupjudgevalue = q_limit
        End Subroutine qcriterion

        !Group Vortex using groupjudgevalue matrix given by specific vortex identify method
        Subroutine groupallvortex(groupvalue, groupjudgevalue)
                Implicit None
                Integer, Intent(Out), Dimension(m, n)::groupvalue
                Real, Intent(In), Dimension(m, n)::groupjudgevalue
                Integer::i, j
                Integer::groupnumber
                groupnumber = 1
                setgroupvalue_m: do i = 1, m
                setgroupvalue_n: Do j = 1, n
                !set only if has groupjudgevalue and hasn't been set
                judgevalue: If (groupjudgevalue(i, j) > 0 .And. groupvalue(i, j) == 0) Then
                        Call setvortexgroupnumber(i, j, groupvalue, groupnumber, groupjudgevalue)
                        groupnumber = groupnumber + 1
                End If judgevalue
                End Do setgroupvalue_n
                End Do setgroupvalue_m
        End Subroutine groupallvortex

        !Recursive setvortexgroupnumber subroutine
        Recursive Subroutine setvortexgroupnumber(i, j, groupvalue, groupnumber, groupjudgevalue)
                Implicit None
                Integer, Intent(In)::i, j, groupnumber
                Integer, Intent(Out), Dimension(m, n)::groupvalue
                Real, Intent(In), Dimension(m, n)::groupjudgevalue
                groupvalue(i, j) = groupnumber
                notatmargin: If (i /= 1 .And. j /= 1 .And. i /= m .And. j /= n) Then
                        iadd: If (groupjudgevalue(i + 1, j) > 0 .And. groupvalue(i + 1, j) == 0) Then
                                Call setvortexgroupnumber(i + 1, j, groupvalue, groupnumber, groupjudgevalue)
                        End If iadd
                        iminus: If (groupjudgevalue(i - 1, j) > 0 .And. groupvalue(i - 1, j) == 0) Then
                                Call setvortexgroupnumber(i - 1, j, groupvalue, groupnumber, groupjudgevalue)
                        End If iminus
                        jadd: If (groupjudgevalue(i, j + 1) > 0 .And. groupvalue(i, j + 1) == 0) Then
                                Call setvortexgroupnumber(i, j + 1, groupvalue, groupnumber, groupjudgevalue)
                        End If jadd
                        jminus: If (groupjudgevalue(i, j - 1) > 0 .And. groupvalue(i, j - 1) == 0) Then
                                Call setvortexgroupnumber(i, j - 1, groupvalue, groupnumber, groupjudgevalue)
                        End If jminus
                End If notatmargin
        End Subroutine setvortexgroupnumber

        !the subroutine definition to record center of vortex, total vorticity and maximumvorticity
        Subroutine vortexcenter(outputvortexarray, groupvalue, vorticity, x, y)
                Implicit None
                Integer::vortex_amount
                Integer, Intent(In), Dimension(:, :)::groupvalue
                Type(vortex), Intent(InOut), Dimension(:)::outputvortexarray                      !an array of type vortex
                Real, Intent(In), Dimension(:, :)::vorticity, x, y
                Integer::i, j
                Integer, Dimension(:, :), Allocatable::maxvorticityij
                Real, Allocatable, Dimension(:, :)::maxvorticity       !record the maximux vorticity of each vortex (to determine position)
                vortex_amount = size(outputvortexarray)
                Allocate (maxvorticity(vortex_amount, 4))
                Allocate (maxvorticityij(vortex_amount, 2))
                maxvorticity = 0
                maxvorticityij = 0

                !judge the vortex center and total voticity and record to outputarray
                setallvortex_m: Do i = 1, m
                setallvortex_n: Do j = 1, n
                ifisinvortex: If (groupvalue(i, j) /= 0) Then
                        maxvorticity(groupvalue(i, j), 3) =&
                                &maxvorticity(groupvalue(i, j), 3) + vorticity(i, j)
                        centerjudge: If (abs(vorticity(i, j)) > abs(maxvorticity(groupvalue(i, j), 4))) Then
                                maxvorticity(groupvalue(i, j), 1) = x(i, j)
                                maxvorticity(groupvalue(i, j), 2) = y(i, j)
                                maxvorticityij(groupvalue(i, j), 1) = i
                                maxvorticityij(groupvalue(i, j), 2) = j
                                maxvorticity(groupvalue(i, j), 4) = vorticity(i, j)
                        End If centerjudge
                End If ifisinvortex
                End Do setallvortex_n
                End Do setallvortex_m

                !interpolation to get more accurate vortexcenter
                Call accurate_center(maxvorticity, maxvorticityij, x, y, vorticity)

                !record maxvorticity to outputvortexarray
                recordtoarray: Do i = 1, vortex_amount
                outputvortexarray(i)%xc = maxvorticity(i, 1)
                outputvortexarray(i)%yc = maxvorticity(i, 2)
                outputvortexarray(i)%vorticity = maxvorticity(i, 3)
                outputvortexarray(i)%maximumvorticity = maxvorticity(i, 4)
                End Do recordtoarray
        End Subroutine vortexcenter

        !the subroutine used in vortexcenter to interpolate to get more accurate vortexcenter
        Subroutine accurate_center(maxvorticity, maxvorticityij, x, y, vorticity)
                Implicit None
                Real, Intent(InOut), Dimension(:, :)::maxvorticity
                Integer, Intent(In), Dimension(:, :)::maxvorticityij
                Real, Intent(In), Dimension(:, :)::x, y, vorticity
                Integer::i, j, k
                Real::dx, dy, zx1, zx2, zx3, zy1, zy2, zy3, xcenter, ycenter
                dx = x(2, 1) - x(1, 1)
                dy = y(1, 2) - y(1, 1)
                forallvortex: Do k = 1, size(maxvorticity(:, 1))
                i = maxvorticityij(k, 1)
                j = maxvorticityij(k, 2)
                zx1 = vorticity(i - 1, j)
                zx2 = vorticity(i, j)
                zx3 = vorticity(i + 1, j)
                zy1 = vorticity(i, j - 1)
                zy2 = vorticity(i, j)
                zy3 = vorticity(i, j + 1)
                ifmaxinmiddle_x: If (zx2 > zx1 .And. zx2 > zx3) Then
                        xcenter = x(i, j) + dx*(zx1 - zx3)/(2*(zx1 - 2*zx2 + zx3))
                        maxvorticity(k, 1) = xcenter
                End If ifmaxinmiddle_x
                ifmaxinmiddle_y: If (zy2 > zy1 .And. zy2 > zy3) Then
                        ycenter = y(i, j) + dy*(zy1 - zy3)/(2*(zy1 - 2*zy2 + zy3))
                        maxvorticity(k, 2) = ycenter
                End If ifmaxinmiddle_y
                End Do forallvortex
        End Subroutine accurate_center

        !the subroutine to record the area of vortex
        Subroutine vortexarea(outputvortexarray, groupvalue)
                Implicit None
                Integer::vortex_amount
                Integer, Intent(In), Dimension(:, :)::groupvalue
                Type(vortex), Intent(InOut), Dimension(:)::outputvortexarray    !an array of type vortex
                Integer, Dimension(:), Allocatable::vortexpoints                !the points number of a vortex
                Real::pointsarea                                                !the area for a point
                Integer::i, j
                vortex_amount = size(outputvortexarray)
                Allocate (vortexpoints(vortex_amount))
                vortexpoints = 0
                pointsarea = (x_max - x_min)*(y_max - y_min)/(m*n)

                !calculate the points of vortex i using groupvalue matrix
                areaxi: Do i = 1, m
                areayj: Do j = 1, n
                areainavortex: If (groupvalue(i, j) /= 0) Then
                        vortexpoints(groupvalue(i, j)) = vortexpoints(groupvalue(i, j)) + 1
                End If areainavortex
                End Do areayj
                End Do areaxi

                !record area to outputvortexarray
                recordareatoarray: Do i = 1, vortex_amount
                outputvortexarray(i)%area = vortexpoints(i)*pointsarea
                End Do recordareatoarray
        End Subroutine vortexarea

        !the subroutine definition to record time, neardistance and direction of vortex
        Subroutine vortexinfo1(outputvortexarray, dataname)
                Implicit None
                Type(vortex), Intent(InOut), Dimension(:)::outputvortexarray                      !an array of type vortex
                Character(len=*), Intent(In)::dataname
                Integer::time, i, vortex_amount
                Real::formerdistance, latterdistance
                vortex_amount = size(outputvortexarray)
                !read time (frame number) from dataname
                Read (dataname((len(dataname) - 8):(len(dataname) - 4)), *) time

                !set the time and direction of vortex
                setvortexforeach: Do i = 1, vortex_amount
                outputvortexarray(i)%time = time
                setvortexdirection: If (outputvortexarray(i)%vorticity > 0) Then
                        outputvortexarray(i)%direction = 1
                Else
                        outputvortexarray(i)%direction = -1
                End If setvortexdirection
                End Do setvortexforeach

                nearvortex: Do i = 1, vortex_amount
                outputvortexarray(i)%neardistance = sqrt(eff_radius**2*3.1415/vortex_amount)
                End Do nearvortex

                !                !set the nearest distance of vortex
                !                formerdistance = 0
                !                latterdistance = 0
                !                nearvortex: Do i = 1, vortex_amount
                !                        ifnotfirstorlastvortex: If (i /= 1 .And. i /= vortex_amount) Then
                !                                formerdistance = sqrt((outputvortexarray(i)%xc - outputvortexarray(i - 1)%xc)**2&
                !                                        &+ (outputvortexarray(i)%yc - outputvortexarray(i - 1)%yc)**2)
                !                                latterdistance = sqrt((outputvortexarray(i)%xc - outputvortexarray(i + 1)%xc)**2&
                !                                        &+ (outputvortexarray(i)%yc - outputvortexarray(i + 1)%yc)**2)
                !                                formerorlatter: If (formerdistance > latterdistance) Then
                !                                        outputvortexarray(i)%neardistance = latterdistance
                !                                Else
                !                                        outputvortexarray(i)%neardistance = formerdistance
                !                                End If formerorlatter
                !                        Else If (i == 1) Then
                !                                latterdistance = sqrt((outputvortexarray(i)%xc - outputvortexarray(i + 1)%xc)**2&
                !                                        &+ (outputvortexarray(i)%yc - outputvortexarray(i + 1)%yc)**2)
                !                                outputvortexarray(i)%neardistance = latterdistance
                !                        Else
                !                                formerdistance = sqrt((outputvortexarray(i)%xc - outputvortexarray(i - 1)%xc)**2&
                !                                        &+ (outputvortexarray(i)%yc - outputvortexarray(i - 1)%yc)**2)
                !                                outputvortexarray(i)%neardistance = formerdistance
                !                        End If ifnotfirstorlastvortex
                !                End Do nearvortex
        End Subroutine vortexinfo1

        !the subroutine to calculate vortex fliter fliterthreshold in find_flitervalue
        Subroutine flitervalue_liutexomega(outputvortexarray, flitervalue)
                Implicit None
                Type(vortex), Dimension(:), Intent(InOut), Allocatable::outputvortexarray
                Real, Intent(Out)::flitervalue
                Integer, Parameter::bins = 10
                Integer, Dimension(bins)::histgramlist
                Integer::i, binposition, flitercount, j
                Real::meanmaximumvorticity

                !calculate average maximumvorticity
                meanmaximumvorticity = 0
                totalmaximumvorticity: Do i = 1, size(outputvortexarray)
                meanmaximumvorticity = meanmaximumvorticity + abs(outputvortexarray(i)%maximumvorticity)
                End Do totalmaximumvorticity
                meanmaximumvorticity = meanmaximumvorticity/size(outputvortexarray)

                !histrogram list and find most probable maximumvorticity
                histgramlist = 0
                gethistrogramlist: Do i = 1, size(outputvortexarray)
                binposition = ceiling(abs(outputvortexarray(i)%maximumvorticity)/(2*meanmaximumvorticity/bins))
                notoutofrange: If (.Not. binposition > bins) Then
                        histgramlist(binposition) = histgramlist(binposition) + 1
                End If notoutofrange
                End Do gethistrogramlist
                findflitervalue: Do i = 1, bins
                foundsuccessful: If (histgramlist(i) == maxval(histgramlist)) Then
                        flitervalue = (1.0*i - 0.5)*(2*meanmaximumvorticity/bins)
                End If foundsuccessful
                End Do findflitervalue
        End Subroutine flitervalue_liutexomega

        !the subroutine of vortex fliter: drop turbulence plume with maximumvorticity less than most probable maximumvorticity
        Subroutine vortexfliter_liutexomega(outputvortexarray, flitervalue)
                Implicit None
                Type(vortex), Dimension(:), Intent(InOut), Allocatable::outputvortexarray
                Type(vortex), Dimension(:), Allocatable::temp_outputvortexarray
                Real, Intent(In)::flitervalue
                Integer::i, flitercount, j

                !write flitered data to temp_outputvortexarray
                flitercount = 0
                j = 0
                countokvortex: Do i = 1, size(outputvortexarray)
                iffliterok: If (abs(outputvortexarray(i)%maximumvorticity) > flitervalue) Then
                        flitercount = flitercount + 1
                End If iffliterok
                End Do countokvortex
                Allocate (temp_outputvortexarray(flitercount))
                writeokvortex_temp: Do i = 1, size(outputvortexarray)
                iffliterok2: If (abs(outputvortexarray(i)%maximumvorticity) > flitervalue) Then
                        j = j + 1
                        temp_outputvortexarray(j) = outputvortexarray(i)
                End If iffliterok2
                End Do writeokvortex_temp

                !substitute outputvortexarray with temp_outputvortexarray
                Deallocate (outputvortexarray)
                Allocate (outputvortexarray(flitercount))
                outputvortexarray = temp_outputvortexarray
        End Subroutine vortexfliter_liutexomega

        !The subroutine for setting vortexdata in outputvortexarray to outputvortexarray_outside, the outer vortexarray should large
        !       enough.
        Subroutine setdatatoarray(outputvortexarray_outside, outputvortexarray)
                Implicit None
                Type(vortex), Dimension(:), Intent(Out)::outputvortexarray_outside
                Type(vortex), Dimension(:), Intent(In)::outputvortexarray
                Integer::i
                i = 1
                writevortex: Do
                if_not_vortex: If (outputvortexarray(i)%direction == 0) Then
                        Exit writevortex
                Else
                        outputvortexarray_outside(i) = outputvortexarray(i)
                        i = i + 1
                End If if_not_vortex
                End Do writevortex
        End Subroutine setdatatoarray

        !The subroutine for writing vortexdata in outputvortexarray to file outputfilename, the dataname is written on first line
        Subroutine writedata(dataname, outputfilename, outputvortexarray)
                Implicit None
                Integer::ierror, i
                Type(vortex), Dimension(:), Intent(In)::outputvortexarray
                Character(len=*), Intent(In)::outputfilename
                Character(len=*), Intent(In)::dataname
                i = 1
                Open (Unit=10, File=outputfilename, Status='New', Iostat=ierror)
                !First and second line is the information of vortex data
                Write (10, "(A27,100a)") 'The data of Vortex of file ', dataname
                Write (10, "(130a)") 'xc(m) yc(m) vorticity(m*s) direction time(frame) nearest distance(m) maximumvorticity(1/s)&
                        & area(m^2)'
                writevortex: Do i = 1, size(outputvortexarray)
                Write (10, "(2F14.7,F20.8,2I7,F14.7,F20.8,F20.10)") outputvortexarray(i)%xc, outputvortexarray(i)%yc,&
                        &outputvortexarray(i)%vorticity, outputvortexarray(i)%direction,&
                        &outputvortexarray(i)%time, outputvortexarray(i)%neardistance,&
                        &outputvortexarray(i)%maximumvorticity, outputvortexarray(i)%area
                End Do writevortex
                Close (Unit=10)
        End Subroutine writedata
End Module find_vortex

!Test Program Using Subroutine findvortex

!Program testfindvortex
!        Use find_vortex
!        Implicit None
!        Type(vortex), Dimension(1000)::testa
!        Call findvortex('B00001.txt', outputvortexarray_outside=testa, outputfilename='vortexdata.txt')
!End Program testfindvortex
