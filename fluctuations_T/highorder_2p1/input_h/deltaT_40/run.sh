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
sh TemB1/run.sh &
#sh TemB2/run.sh &
#sh TemB3/run.sh &
#sh TemB4/run.sh &
#sh TemB5/run.sh &
#sh TemB6/run.sh &
#sh TemB7/run.sh &
#sh TemB8/run.sh &
#sh TemB9/run.sh &
#sh TemB10/run.sh &
#sh TemB11/run.sh &
#sh TemB12/run.sh &
#sh TemB13/run.sh &
#sh TemB14/run.sh &
#sh TemB15/run.sh &
#sh TemB16/run.sh &
#sh TemB17/run.sh &
#sh TemB18/run.sh &
#sh TemB19/run.sh &
#sh TemB20/run.sh &
#sh TemB21/run.sh &
wait
#for a in {1..21} 
#do
#mv TemB$a/BUFFER/VTOTAL.DAT dataB/V$a.dat
#done