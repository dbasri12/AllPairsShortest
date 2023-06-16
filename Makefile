ifdef DEBUG
MACRO = -DPRINT_DEBUG
endif

#compiler setup
CXX = g++
CXXFLAGS = -std=c++14 -O3 -pthread $(MACRO)

MPICXX = mpic++
MPI_CXXFLAGS = -std=c++14 -O3 $(MACRO)

COMMON= core/utils.h core/cxxopts.h core/get_time.h core/graph.h core/quick_sort.h
SERIAL= SSSP
PARALLEL= SSSP_parallel
MPI = SSSP_MPI
ALL= $(SERIAL) $(PARALLEL) $(MPI)


all : $(ALL)

$(SERIAL): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(PARALLEL): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(MPI): %: %.cpp
	$(MPICXX) $(CXXFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
