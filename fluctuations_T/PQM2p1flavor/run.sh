#!/bin/bash
#source /opt/intel/oneapi/setvars.sh intel64
rm exam/eqcd
rm exam/*.out
rm exam/*.o
rm exam/*.mod
rm exam/BUFFER/*
cd exam
make
rm *.o
rm *.mod
cd ..
rm -rf data
rm -rf dataB/*
#rm -rf dataS/*
#rm -rf dataBS/*
#rm -rf data_sca/*
rm -rf Tem*
for a in {1..1} 
do
	mkdir TemB$a
#	mkdir TemS$a
#	mkdir TemBS$a
	cp -r exam/qm2p1 TemB$a
#	cp -r exam/qm2p1 TemS$a
#	cp -r exam/qm2p1 TemBS$a
	mkdir TemB$a/BUFFER
#	mkdir TemS$a/BUFFER
#	mkdir TemBS$a/BUFFER
	echo -e "$a\n1" >TemB$a/M1.DAT
#	echo -e "$a\n2" >TemS$a/M1.DAT
#	echo -e "$a\n3" >TemBS$a/M1.DAT
	echo -e "#!/bin/bash\ncd TemB$a\n ./qm2p1" >TemB$a/run.sh
#	echo -e "#!/bin/bash\ncd TemS$a\n ./qm2p1" >TemS$a/run.sh
#	echo -e "#!/bin/bash\ncd TemBS$a\n ./qm2p1" >TemBS$a/run.sh
done
sh exam/run2.sh &
wait
#for a in {1..1} 
#do
#cp TemB$a/BUFFER/VTOTAL.DAT dataB/V$a.dat
#cp TemS$a/BUFFER/VTOTAL.DAT dataS/V$a.dat
#cp TemBS$a/BUFFER/VTOTAL.DAT dataBS/V$a.dat
#done