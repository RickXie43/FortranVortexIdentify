#源文件所在目录
SOURCE:=source

#并行程序编译
vortex_mpi_cmd : b1.o b2.o b3.o b4.o b5.o b6.o fvm.o f.o
	mpiifort b1.o b2.o b3.o b4.o b5.o b6.o fvm.o f.o -o vortex_mpi_cmd

f.o : $(SOURCE)/main_findvortex_MPI_CMD.f90
	mpiifort -c $(SOURCE)/main_findvortex_MPI_CMD.f90 -o f.o -I.

fvm.o : $(SOURCE)/find_vortex_module.f90 $(SOURCE)/parameters.F90
ifdef DATADIR
	./$(SOURCE)/setnxny.sh $(DATADIR)
endif
	mpiifort -fpp -c $(SOURCE)/find_vortex_module.f90 -o fvm.o -I.

b3.o : $(SOURCE)/bspline_module.f90 b5.o
	mpiifort -c ./$(SOURCE)/bspline_module.f90 -o b3.o -I.

b5.o : $(SOURCE)/bspline_defc_module.F90 b6.o
	mpiifort -c ./$(SOURCE)/bspline_defc_module.F90 -o b5.o -I.

b6.o : $(SOURCE)/bspline_oo_module.f90 b2.o
	mpiifort -c ./$(SOURCE)/bspline_oo_module.f90 -o b6.o -I.

b2.o : $(SOURCE)/bspline_sub_module.f90 b1.o
	mpiifort -c ./$(SOURCE)/bspline_sub_module.f90 -o b2.o -I.

b1.o : $(SOURCE)/bspline_blas_module.F90 b4.o
	mpiifort -c ./$(SOURCE)/bspline_blas_module.F90 -o b1.o -I.

b4.o : $(SOURCE)/bspline_kinds_module.F90
	mpiifort -c ./$(SOURCE)/bspline_kinds_module.F90 -o b4.o -I.

#清理编译文件
clean :
	rm -f *.o *.mod vortex_mpi_cmd test vortex_single

#编译测试文件
test : b1.o b2.o b3.o b4.o b5.o b6.o fvm.o test.o
	ifort b1.o b2.o b3.o b4.o b5.o b6.o fvm.o test.o -o test

test.o : $(SOURCE)/main_test.f90
	ifort -c ./$(SOURCE)/main_test.f90 -o test.o -I.

#单核程序编译
single_processor : b1_single.o b2_single.o b3_single.o b4_single.o b5_single.o b6_single.o fvm_single.o single.o
	ifort b1_single.o b2_single.o b3_single.o b4_single.o b5_single.o b6_single.o fvm_single.o single.o -o vortex_single

single.o : $(SOURCE)/main_findvortex.f90
	ifort -c ./$(SOURCE)/main_findvortex.f90 -o single.o -I.

fvm_single.o : $(SOURCE)/find_vortex_module.f90 $(SOURCE)/parameters.F90
ifdef DATADIR
	./$(SOURCE)/setnxny.sh $(DATADIR)
endif
	ifort -fpp -c $(SOURCE)/find_vortex_module.f90 -o fvm_single.o -I.

b3_single.o : $(SOURCE)/bspline_module.f90 b5_single.o
	ifort -c ./$(SOURCE)/bspline_module.f90 -o b3_single.o -I.

b5_single.o : $(SOURCE)/bspline_defc_module.F90 b6_single.o
	ifort -c ./$(SOURCE)/bspline_defc_module.F90 -o b5_single.o -I.

b6_single.o : $(SOURCE)/bspline_oo_module.f90 b2_single.o
	ifort -c ./$(SOURCE)/bspline_oo_module.f90 -o b6_single.o -I.

b2_single.o : $(SOURCE)/bspline_sub_module.f90 b1_single.o
	ifort -c ./$(SOURCE)/bspline_sub_module.f90 -o b2_single.o -I.

b1_single.o : $(SOURCE)/bspline_blas_module.F90 b4_single.o
	ifort -c ./$(SOURCE)/bspline_blas_module.F90 -o b1_single.o -I.

b4_single.o : $(SOURCE)/bspline_kinds_module.F90
	ifort -c ./$(SOURCE)/bspline_kinds_module.F90 -o b4_single.o -I.
