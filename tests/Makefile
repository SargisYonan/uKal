cc=cc
cc_options=-Wall -Wextra -Werror -std=c99 -Os

ukal_src_dir=../src

ukal_options=-DUKAL_MAX_STATE_VECTOR_SIZE=5 -DUKAL_MAX_MEASUREMENT_VECTOR_SIZE=5

ulapack_options=-DULAPACK_MATRIX_ENTRY_TYPE=double -DULAPACK_INDEX_TYPE=uint64_t \
			 -DULAPACK_SIGNED_INDEX_TYPE=int64_t -DULAPACK_USE_STATIC_ALLOC=1 \
			 -DULAPACK_MAX_MATRIX_N_ROWS=5 -DULAPACK_MAX_MATRIX_N_COLS=5 \
			 -DULAPACK_INITIALIZE_MEMORY=1 -DULAPACK_INVERSE_VIA_LU=1 -DULAPACK_MINIMUM_THRESHOLD_TOLERANCE=.0001

ulapack_src_dir=$(ukal_src_dir)/uLAPack/src

exec=ukal_tests

all : $(exec)

$(exec): unit_test.c ukal.o
	$(cc) $(cc_options) -I$(ukal_src_dir) -I$(ulapack_src_dir) $(ukal_options) $(ulapack_options) ukal.o ulapack.o -lm unit_test.c -o $(exec)

ukal.o: ulapack.o $(ukal_src_dir)/ukal.c $(ukal_src_dir)/ukal.c $(ukal_src_dir)/ukal_type.h
	$(cc) $(cc_options) -I$(ukal_src_dir) -I$(ulapack_src_dir) -c $(ukal_src_dir)/ukal.c $(ukal_options) $(ulapack_options) -o ukal.o

ulapack.o:
	$(cc) $(cc_options) -I$(ulapack_src_dir) -c $(ulapack_src_dir)/ulapack.c  $(ulapack_options) -o ulapack.o

clean:
	rm ukal.o ulapack.o $(exec)