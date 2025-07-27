OBJ=eQCD.f90 modsave.o gauleg.o modcal.o initial.o phypoint.o fun.o nfdx_com.o lblf.o Intquarkloop.o thresholdfunc.o derivs.o odeint.o dtVdiff1.o polyakoveq.o intlin.o selfeq.o gapeq.o newt.o 
F77=gfortran

exe: $(OBJ)
	$(F77) -o eqcd -O3 -fopenmp -falign-functions=16 $(OBJ) $(LIBS)

modsave.o: modsave.f90
	$(F77)  -c -O3 -falign-functions=16 modsave.f90

modcal.o: modcal.f90
	$(F77)  -c -O3 -falign-functions=16 modcal.f90

Intquarkloop.o: Intquarkloop.f90
	$(F77)  -c -O3 -fopenmp -falign-functions=16 Intquarkloop.f90

derivs.o: derivs.f90
	$(F77)  -c -O3 -falign-functions=16 derivs.f90

dtVdiff1.o: dtVdiff1.f90
	$(F77)  -c -O3 -falign-functions=16 dtVdiff1.f90

fun.o: fun.f90
	$(F77)  -c -O3 -falign-functions=16 fun.f90

gapeq.o: gapeq.f90
	$(F77)  -c -O3 -falign-functions=16 gapeq.f90

gauleg.o: gauleg.f90
	$(F77)  -c -O3 -falign-functions=16 gauleg.f90

initial.o: initial.f90
	$(F77)  -c -O3 -falign-functions=16 initial.f90

intlin.o: intlin.f90
	$(F77)  -c -O3 -falign-functions=16 intlin.f90

lblf.o: lblf.f90
	$(F77)  -c -O3 -falign-functions=16 lblf.f90

newt.o: newt.f90
	$(F77)  -c -O3 -falign-functions=16 newt.f90

nfdx_com.o: nfdx_com.f90
	$(F77)  -c -O3 -falign-functions=16 nfdx_com.f90

odeint.o: odeint.f90
	$(F77)  -c -O3 -falign-functions=16 odeint.f90

phypoint.o: phypoint.f90
	$(F77)  -c -O3 -falign-functions=16 phypoint.f90

polyakoveq.o: polyakoveq.f90
	$(F77)  -c -O3 -fopenmp -falign-functions=16 polyakoveq.f90

selfeq.o: selfeq.f90
	$(F77)  -c -O3 -falign-functions=16 selfeq.f90

thresholdfunc.o: thresholdfunc.f90
	$(F77)  -c -O3 -falign-functions=16 thresholdfunc.f90
