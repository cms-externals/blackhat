FC=gfortran
BH_PATH=$(prefix)


VPATH = ${srcdir}

fortran_example.o:	$(srcdir)/fortran_example.f
	$(FC) -c -g $(srcdir)/fortran_example.f -o $@

cpp_example.o:	$(srcdir)/cpp_example.cpp
	$(CXX) -c -g $(srcdir)/cpp_example.cpp -o $@



fortran_example:	fortran_example	
	$(FC) -g -o $@ fortran_example.o -lc -lstdc++ $$($(BH_PATH)/bin/blackhat-config --libs)  

cpp_example:	cpp_example.o	
	$(CXX) -o $@ cpp_example.o $$($(BH_PATH)/bin/blackhat-config --libs)  
	
	
	
EXTRA_dist = contract_file.lh BHLH_example.py order_example.lh 
	
examplesdir = $(pkgdatadir)/examples

dist_examples_DATA = fortran_example.f cpp_example.cpp BHLH_example.py contract_file.lh Sherpa_W0_example.dat Sherpa_W2_example.dat
