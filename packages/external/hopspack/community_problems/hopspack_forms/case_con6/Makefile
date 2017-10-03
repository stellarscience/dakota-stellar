
# Compiler name
F77    = g77
INCLUDE_DIR = -I. -I/usr/include
LIB_DIR = . 
#LIBS = -llammpi++ -llammpio -lpmpi -llamf77mpi -lmpi -llam -lnsl 


#*******************************************************************

PGM    = EvalFn
OBJS   = $(PGM).o \
         location_constraint.o linear_constraint.o\
         head_constraint.o\
	 costac.o costao.o\
         feval.o writewel.o readhead.o

EXE    = evalCON6

# Rules
all: $(EXE)

$(EXE): $(OBJS)
	@echo
	$(F77) $(INCLUDE_DIR) $(FFLAGS) $(OBJS)  -L$(LIB_DIR) $(LIBS) -o $@
#	strip $@

clean:
	rm -f $(EXE) *.o

test:
	rm -f A1-1.hed A1-1.mfn A1-1.out A1-1.wel
	rm -f input.1 test_output.1
	./evaljob_con6.sh test_input.1 test_output.1 1

.f.o:
	$(F77) $(INCLUDE_DIR) $(FFLAGS) $(PG) -c $<
