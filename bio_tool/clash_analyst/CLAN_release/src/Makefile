CPP=g++
CFLAGS = -std=c++0x -fopenmp -O2
OBJS = bio_alphabet.o bwt_frag_merger.o bwt.o bwt_search.o concatenator.o divsufsort.o kmer_unitcoder.o loader.o minimizer_sort.o sssort.o trsort.o utils.o mapping_driver.o output_printer.o read_merge.o dimer_folding.o
all: $(OBJS) main_index.cc main_search.cc main_output.cc main_annotate.cc
	$(CPP) $(CFLAGS) main_index.cc $(OBJS) -o clan_index
	$(CPP) $(CFLAGS) main_search.cc $(OBJS) -o clan_search
	$(CPP) $(CFLAGS) main_output.cc $(OBJS) -o clan_output
	$(CPP) $(CFLAGS) main_annotate.cc $(OBJS) -o clan_annotate
clan_index:	$(OBJS) main_index.cc parameters.h
	$(CPP) $(CFLAGS) main_index.cc $(OBJS) -o clan_index
clan_search: $(OBJS) main_search.cc parameters.h misc.h
	$(CPP) $(CFLAGS) main_search.cc $(OBJS) -o clan_search
clan_output: $(OBJS) main_output.cc parameters.h output_driver.h
	$(CPP) $(CFLAGS) main_output.cc $(OBJS) -o clan_output
clan_annotate: $(OBJS) main_annotate.cc parameters.h annotate_driver.h misc.h
	$(CPP) $(CFLAGS) main_annotate.cc $(OBJS) -o clan_annotate
dimer_folding.o: dimer_folding.cc dimer_folding.h free_energy.h
	$(CPP) $(CFLAGS) -c dimer_folding.cc -o dimer_folding.o
read_merge.o: read_merge.cc read_merge.h bwt_frag_merger.h
	$(CPP) $(CFLAGS) -c read_merge.cc -o read_merge.o
output_printer.o: output_printer.cc mapping_driver.h parameters.h bwt.o bwt_search.o bwt_frag_merger.o
	$(CPP) $(CFLAGS) -c output_printer.cc -o output_printer.o
mapping_driver.o: mapping_driver.cc mapping_driver.h parameters.h bwt.o bwt_search.o bwt_frag_merger.o
	$(CPP) $(CFLAGS) -c mapping_driver.cc -o mapping_driver.o
bio_alphabet.o: bio_alphabet.cc
	$(CPP) $(CFLAGS) -c bio_alphabet.cc -o bio_alphabet.o
bwt.o: bwt.cc divsufsort.o bio_alphabet.o concatenator.o
	$(CPP) $(CFLAGS) -c bwt.cc -o bwt.o
bwt_search.o: bwt_search.cc bwt.o seed_extension.h
	$(CPP) $(CFLAGS) -c bwt_search.cc -o bwt_search.o
bwt_frag_merger.o: bwt_frag_merger.cc bwt.o bwt_search.o
	$(CPP) $(CFLAGS) -c bwt_frag_merger.cc -o bwt_frag_merger.o
concatenator.o: concatenator.cc
	$(CPP) $(CFLAGS) -c concatenator.cc -o concatenator.o
kmer_unitcoder.o: kmer_unitcoder.cc bio_alphabet.o
	$(CPP) $(CFLAGS) -c kmer_unitcoder.cc -o kmer_unitcoder.o
loader.o: loader.cc bio_alphabet.o
	$(CPP) $(CFLAGS) -c loader.cc -o loader.o
minimizer_sort.o: minimizer_sort.cc kmer_unitcoder.o
	$(CPP) $(CFLAGS) -c minimizer_sort.cc -o minimizer_sort.o
divsufsort.o: divsufsort.cc 
	$(CPP) $(CFLAGS) -c divsufsort.cc -o divsufsort.o
sssort.o: sssort.cc divsufsort.o
	$(CPP) $(CFLAGS) -c sssort.cc -o sssort.o
trsort.o: trsort.cc divsufsort.o
	$(CPP) $(CFLAGS) -c trsort.cc -o trsort.o
utils.o: utils.cc divsufsort.o
	$(CPP) $(CFLAGS) -c utils.cc -o utils.o
clean:
	rm *.o
new:
	touch *.cc
