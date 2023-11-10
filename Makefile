vortex_mpi_cmd : b1.o b2.o b3.o b4.o b5.o b6.o fvm.o f.o
	mpiifort b1.o b2.o b3.o b4.o b5.o b6.o fvm.o f.o -o vortex_mpi_cmd

f.o : main_findvortex_MPI_CMD.f90
	mpiifort -c main_findvortex_MPI_CMD.f90 -o f.o -I.
fvm.o : find_vortex_module.f90 parameters.F90
ifdef DATADIR
	./setnxny.sh $(DATADIR)
endif
	mpiifort -fpp -c find_vortex_module.f90 -o fvm.o -I.
b3.o : bspline_module.f90 b5.o
	mpiifort -c ./bspline_module.f90 -o b3.o -I.
b5.o : bspline_defc_module.F90 b6.o
	mpiifort -c ./bspline_defc_module.F90 -o b5.o -I.
b6.o : bspline_oo_module.f90 b2.o
	mpiifort -c ./bspline_oo_module.f90 -o b6.o -I.
b2.o : bspline_sub_module.f90 b1.o
	mpiifort -c ./bspline_sub_module.f90 -o b2.o -I.
b1.o : bspline_blas_module.F90 b4.o
	mpiifort -c ./bspline_blas_module.F90 -o b1.o -I.
b4.o : bspline_kinds_module.F90
	mpiifort -c ./bspline_kinds_module.F90 -o b4.o -I.
clean : 
	rm -f *.o *.mod vortex_mpi_cmd test vortex_single

test : b1.o b2.o b3.o b4.o b5.o b6.o fvm.o test.o
	ifort b1.o b2.o b3.o b4.o b5.o b6.o fvm.o test.o -o test
test.o : main_test.f90
	ifort -c ./main_test.f90 -o test.o -I.

single_processor : b1.o b2.o b3.o b4.o b5.o b6.o fvm.o single.o
	ifort b1.o b2.o b3.o b4.o b5.o b6.o fvm.o single.o -o vortex_single
single.o : main_findvortex.f90
	ifort -c ./main_findvortex.f90 -o single.o -I.

