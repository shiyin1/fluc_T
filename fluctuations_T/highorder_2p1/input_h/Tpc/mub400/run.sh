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
rm -rf dataB/*
rm -rf Tem*
for a in {1..1}
do
	mkdir TemB$a
	cp -r exam/qm2p1 TemB$a
	mkdir TemB$a/BUFFER
	echo -e "$a\n1" >TemB$a/M1.DAT
	echo -e "#!/bin/bash\ncd TemB$a\n ./qm2p1" >TemB$a/run.sh
done
sh TemB1/run.sh 
wait