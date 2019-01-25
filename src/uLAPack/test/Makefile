cc=cc
cc_options=-std=c99 -Wall -Wextra -Werror -Os

static_test_exec=ulapack_test_static
static_test_options=-DULAPACK_USE_STATIC_ALLOC=1 -DULAPACK_MAX_MATRIX_N_COLS=40 -DULAPACK_MAX_MATRIX_N_ROWS=40

dynamic_test_exec=ulapack_test_dynamic
dynamic_test_options=-DULAPACK_USE_DYNAMIC_ALLOC=1

src=../src

all : $(static_test_exec) $(dynamic_test_exec)

test : $(static_test_exec) $(dynamic_test_exec)
	./$(static_test_exec)
	./$(dynamic_test_exec)

$(static_test_exec) : static_unit_tests.o ulapack_static.o
	$(cc) $(cc_options) unit_tests.o ulapack.o -lm -o $(static_test_exec)

$(dynamic_test_exec) : dynamic_unit_tests.o ulapack_dynamic.o
	$(cc) $(cc_options) unit_tests.o ulapack.o -lm -o $(dynamic_test_exec)

dynamic_unit_tests.o : unit_tests.c ulapack_dynamic.o
	$(cc) $(cc_options) -I$(src) -c unit_tests.c $(dynamic_test_options)

static_unit_tests.o : unit_tests.c ulapack_static.o
	$(cc) $(cc_options) -I$(src) -c unit_tests.c $(static_test_options)

ulapack_static.o : $(src)/ulapack.c ulapack_headers
	$(cc) $(cc_options) -c $(src)/ulapack.c $(static_test_options)

ulapack_dynamic.o : $(src)/ulapack.c ulapack_headers
	$(cc) $(cc_options) -c $(src)/ulapack.c $(dynamic_test_options)

ulapack_headers : $(src)/ulapack.h $(src)/ulapack_options.h $(src)/ulapack_type.h

clean :
	rm *.o $(dynamic_test_exec) $(static_test_exec)