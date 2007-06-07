
# The source file where the main() function is

SOURCEMAIN = ukb_aw.cc

# Library files

SRC = w2syn.c mcrGraph.c disambGraph.c csentence.c kGraph.c

# Don't change anything below
#DEBUG = 1
#PROF = 1

INCLUDE_DIR = 
LIBS = -lboost_filesystem -lboost_program_options


ifdef DEBUG
OPTFLAGS = -g
else
OPTFLAGS = -O2
endif

ifdef PROF
PROFFLAGS = -pg
else
PROFFLAGS =
endif

CCOPTIONS = -Wall $(OPTFLAGS)
MEMBERS = $(SRC:.c=.o)
EXEC  = $(basename $(notdir $(SOURCEMAIN)))

all: $(EXEC)

%.o : %.cc
	g++ -c $(CCOPTIONS) $(PROFFLAGS) -o $@  $(INCLUDE_DIR) $< 

$(EXEC): $(TARGET) $(MEMBERS) $(SOURCEMAIN)
	gcc $(CCOPTIONS) $(PROFFLAGS) -o $(EXEC) $(SOURCEMAIN) $(MEMBERS) $(INCLUDE_DIR) $(LIBDIR) $(LIBS)

.PHONY : all clean

clean:
	find . -type f -name '*.o' | xargs rm -f
	rm -f $(EXEC)

