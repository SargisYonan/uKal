cc=cc
cc_options=-Wall -Werror -Wextra -Os

ulap_src_dir=../src
ulap_opt=-DULAPACK_USE_STATIC_ALLOC=1 -DULAPACK_MAX_MATRIX_N_ROWS=20 -DULAPACK_MAX_MATRIX_N_COLS=20
ulap_headers=$(ulap_src_dir)/ulapack.h $(ulap_src_dir)/ulapack_options.h

all: pca polyfit

pca : ulapack.o polyfit.c
	cc -I$(ulap_src_dir) $(cc_options) $(ulap_opt) ulapack.o pca.c -o pca

polyfit : ulapack.o polyfit.c
	cc -I$(ulap_src_dir) $(cc_options) $(ulap_opt) ulapack.o polyfit.c -o polyfit

ulapack.o : $(ulap_src_dir)/ulapack.c $(ulap_headers)
	cc -I$(ulap_src_dir) $(cc_options) $(ulap_opt) -c $(ulap_src_dir)/ulapack.c -o ulapack.o

clean :
	rm *.o polyfit pca