CC	      = gcc
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
BAS		= basicmodel.o io.o utilities.o distributions.o likelihoods.o variants.o diagnostics.o thresholds.o threshold_data.o checkpointing.o process_likelihoods.o find_trans_networks.o find_adjacent.o convert_likelihoods_alt.o convert_checkpoint_alt.o find_orders.o calc_timings.o

basic: $(BAS)
	$(CC) $(CC_FLAGS) $(BAS) -o run_basic  $(LD_FLAGS)
basicmodel.o: basicmodel.cpp
	$(CC) $(CC_FLAGS) -c basicmodel.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp
distributions.o: distributions.cpp
	 $(CC) $(CC_FLAGS) -c distributions.cpp
likelihoods.o: likelihoods.cpp
	$(CC) $(CC_FLAGS) -c likelihoods.cpp
variants.o: variants.cpp
	$(CC) $(CC_FLAGS) -c variants.cpp
diagnostics.o: diagnostics.cpp
	$(CC) $(CC_FLAGS) -c diagnostics.cpp
thresholds.o: thresholds.cpp
	$(CC) $(CC_FLAGS) -c thresholds.cpp
threshold_data.o: threshold_data.cpp
	$(CC) $(CC_FLAGS) -c threshold_data.cpp
checkpointing.o: checkpointing.cpp
	$(CC) $(CC_FLAGS) -c checkpointing.cpp
process_likelihoods.o: process_likelihoods.cpp
	$(CC) $(CC_FLAGS) -c process_likelihoods.cpp
find_trans_networks.o: find_trans_networks.cpp
	$(CC) $(CC_FLAGS) -c find_trans_networks.cpp
find_adjacent.o: find_adjacent.cpp
	$(CC) $(CC_FLAGS) -c find_adjacent.cpp
convert_likelihoods_alt.o: convert_likelihoods_alt.cpp
	$(CC) $(CC_FLAGS) -c convert_likelihoods_alt.cpp
convert_checkpoint_alt.o: convert_checkpoint_alt.cpp
	$(CC) $(CC_FLAGS) -c convert_checkpoint_alt.cpp
find_orders.o: find_orders.cpp
	$(CC) $(CC_FLAGS) -c find_orders.cpp
calc_timings.o: calc_timings.cpp
	$(CC) $(CC_FLAGS) -c calc_timings.cpp

