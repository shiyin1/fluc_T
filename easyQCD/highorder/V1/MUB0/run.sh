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
for a in {1..21} 
do
	mkdir Tem$a
	cp -r exam/eqcd Tem$a
	mkdir Tem$a/BUFFER
	echo -e "$a\n1" >Tem$a/M1.DAT
	echo -e "#!/bin/bash\ncd Tem$a\n ./eqcd" >Tem$a/run.sh
done
sh Tem1/run.sh &
sh Tem2/run.sh &
sh Tem3/run.sh &
sh Tem4/run.sh &
sh Tem5/run.sh &
sh Tem6/run.sh &
sh Tem7/run.sh &
sh Tem8/run.sh &
sh Tem9/run.sh &
sh Tem10/run.sh &
sh Tem11/run.sh &
sh Tem12/run.sh &
sh Tem13/run.sh &
sh Tem14/run.sh &
sh Tem15/run.sh &
sh Tem16/run.sh &
sh Tem17/run.sh &
sh Tem18/run.sh &
sh Tem19/run.sh &
sh Tem20/run.sh &
sh Tem21/run.sh &
wait
for a in {1..21} 
do
mv Tem$a/BUFFER/VTOTAL.DAT data/V$a.dat
done