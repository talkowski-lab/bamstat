CFLAGS = -I /source/samtools-0.1.18 -L /source/samtools-0.1.18 -lbam -lz

bamstat: bamstat.cpp alignment_stats.o orientation.o samstring.o running_mean.o running_median.o
	g++ -g -o bamstat bamstat.cpp alignment_stats.o orientation.o samstring.o running_mean.o running_median.o $(CFLAGS)

alignment_stats.o: alignment_stats.cpp orientation.cpp samstring.cpp running_median.cpp running_mean.cpp
	g++ -g -c alignment_stats.cpp orientation.cpp samstring.cpp running_median.cpp running_mean.cpp  -I /usr/local/src/samtools-0.1.18

running_mean.o: running_mean.cpp
	g++ -c running_mean.cpp

running_median.o: running_median.cpp
	g++ -c running_median.cpp

orientation.o: orientation.cpp
	g++ -c orientation.cpp $(CFLAGS)

samstring.o: samstring.cpp
	g++ -c samstring.cpp $(CFLAGS)

clean:
	rm *.o
	rm bamstat
