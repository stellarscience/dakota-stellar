EXE                 = basic_lib.exe

COMPILATOR          = g++

COMPILATOR_OPTIONS  = -O2 -ansi 
L1                  = $(NOMAD_HOME)/lib/nomad.a
LIBS                = $(L1) -lm
INCLUDE             = -I$(NOMAD_HOME)/src -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = basic_lib.o


ifndef NOMAD_HOME
define ECHO_NOMAD
	@echo Please set NOMAD_HOME environment variable!
	@false
endef
endif


$(EXE): $(L1) $(OBJS)
	$(ECHO_NOMAD)
	@echo "   building the parallel (ptrhreads) version ..."
	@echo "   exe file : "$(EXE)
	@$(COMPILATOR) -o  $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)
	@strip  $(EXE)

basic_lib.o: basic_lib.cpp
	$(ECHO_NOMAD)
	@$(COMPILE) basic_lib.cpp

$(L1) : ;
	$(ECHO_NOMAD)


all: $(EXE) 

clean: ;
	@echo "   cleaning obj files"
	@rm -f $(OBJS) 

del: ;
	@echo "   cleaning trash files"
	@rm -f core *~
	@echo "   cleaning obj files"
	@rm -f $(OBJS)
	@echo "   cleaning exe file"
	@rm -f $(EXE)


