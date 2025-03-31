OBJ=QM.f90 chebft.o fityukawa.o derivs.o dtVdiff1.o dtVdiff2.o FRG.o fun.o gapEq.o GAULEG.o initial.o intLin.o massbose2.o nbose.o newt.o nfdx.o nfplnfmi.o odeint.o phypoint2.o PolyakovEq.o selfEQ.o vInf.o 
F77=gfortran

exe: $(OBJ)
	$(F77) -o qm2p1 -O2 $(OBJ) $(LIBS)

chebft.o: chebft.f90
	$(F77)  -c -O2 chebft.f90

derivs.o: derivs.f90
	$(F77)  -c -O2 derivs.f90

dtVdiff1.o: dtVdiff1.f90
	$(F77)  -c -O2 dtVdiff1.f90

dtVdiff2.o: dtVdiff2.f90
	$(F77)  -c -O2 dtVdiff2.f90

FRG.o: FRG.f90
	$(F77)  -c -O2 FRG.f90

fun.o: fun.f90
	$(F77)  -c -O2 fun.f90

fityukawa.o: fityukawa.f90
	$(F77)  -c -O2 fityukawa.f90

gapEq.o: gapEq.f90
	$(F77)  -c -O2 gapEq.f90

GAULEG.o: GAULEG.f90
	$(F77)  -c -O2 GAULEG.f90

initial.o: initial.f90
	$(F77)  -c -O2 initial.f90

intLin.o: intLin.f90
	$(F77)  -c -O2 intLin.f90

massbose2.o: massbose2.f90
	$(F77)  -c -O2 massbose2.f90

nbose.o: nbose.f90
	$(F77)  -c -O2 nbose.f90

newt.o: newt.f90
	$(F77)  -c -O2 newt.f90

nfdx.o: nfdx.f90
	$(F77)  -c -O2 nfdx.f90

nfplnfmi.o: nfplnfmi.f90
	$(F77)  -c -O2 nfplnfmi.f90

odeint.o: odeint.f90
	$(F77)  -c -O2 odeint.f90

phypoint2.o: phypoint2.f90
	$(F77)  -c -O2 phypoint2.f90

PolyakovEq.o: PolyakovEq.f90
	$(F77)  -c -O2 PolyakovEq.f90

selfEQ.o: selfEQ.f90
	$(F77)  -c -O2 selfEQ.f90

vInf.o: vInf.f90
	$(F77)  -c -O2 vInf.f90
