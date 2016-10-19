all: quantification

quantification:
	 g++ -std=gnu++11 -I/usr/local/include/bamtools -o mapping-from-bam-file mapping-from-bam-file.cpp -lbamtools
