
SRCFILES = collision.o finite_volume.o grid.o integrals.o quadrature.o simulation.o util.o

# uncomment to use Intel compiler
# CC = icc

# set these with appropriate libs for your system
EXTRALIBS = -lm -lfftw3
# EXTRALIBS = -lfftw3

# compiler options
CCOPTS = -Wall -m64 -O3 -DNDEBUG -Drestrict_=__restrict__


all: proj_main proj_test

proj_main: proj_main_hom proj_main_inhom

proj_test: proj_Cd_test proj_Cc_test proj_coll_step_test proj_I1_test proj_I2_test proj_I3_test proj_I4_test proj_lgwt_test proj_sim_step_test proj_slope_lim_test proj_slope_lim_test2


proj_main_hom: ${SRCFILES} main_hom.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o main_hom $?

proj_main_inhom: ${SRCFILES} main_inhom.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o main_inhom $?


proj_Cd_test: ${SRCFILES} Cd_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o Cd_test $?

proj_Cc_test: ${SRCFILES} Cc_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o Cc_test $?

proj_coll_step_test: ${SRCFILES} coll_step_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o coll_step_test $?

proj_I1_test: ${SRCFILES} I1_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o I1_test $?

proj_I2_test: ${SRCFILES} I2_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o I2_test $?

proj_I3_test: ${SRCFILES} I3_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o I3_test $?

proj_I4_test: ${SRCFILES} I4_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o I4_test $?

proj_lgwt_test: ${SRCFILES} lgwt_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o lgwt_test $?

proj_sim_step_test: ${SRCFILES} sim_step_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o sim_step_test $?

proj_slope_lim_test: ${SRCFILES} slope_lim_test.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o slope_lim_test $?

proj_slope_lim_test2: ${SRCFILES} slope_lim_test2.o
	${CC} ${CCOPTS} ${EXTRALIBS} -o slope_lim_test2 $?


Cd_test.o: ../test/Cd_test.c
	${CC} -c -I../include ${CCOPTS} $<

Cc_test.o: ../test/Cc_test.c
	${CC} -c -I../include ${CCOPTS} $<

coll_step_test.o: ../test/coll_step_test.c
	${CC} -c -I../include ${CCOPTS} $<

I1_test.o: ../test/I1_test.c
	${CC} -c -I../include ${CCOPTS} $<

I2_test.o: ../test/I2_test.c
	${CC} -c -I../include ${CCOPTS} $<

I3_test.o: ../test/I3_test.c
	${CC} -c -I../include ${CCOPTS} $<

I4_test.o: ../test/I4_test.c
	${CC} -c -I../include ${CCOPTS} $<

lgwt_test.o: ../test/lgwt_test.c
	${CC} -c -I../include ${CCOPTS} $<

sim_step_test.o: ../test/sim_step_test.c
	${CC} -c -I../include ${CCOPTS} $<

slope_lim_test.o: ../test/slope_lim_test.c
	${CC} -c -I../include ${CCOPTS} $<

slope_lim_test2.o: ../test/slope_lim_test2.c
	${CC} -c -I../include ${CCOPTS} $<


%.o: ../src/%.c
	${CC} -c -I../include ${CCOPTS} $<


clean:
	rm *.o
