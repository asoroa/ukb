
# The source file where the main() function is

SOURCEMAIN = ukb_aw.cc

# Library files

SRC = common.cc globalVars.cc configFile.cc fileElem.cc w2syn.cc mcrGraph.cc disambGraph.cc csentence.cc kGraph.cc

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
MEMBERS = $(SRC:.cc=.o)
EXEC  = $(basename $(notdir $(SOURCEMAIN)))

all: $(EXEC)

%.o : %.cc
	g++ -c $(CCOPTIONS) $(PROFFLAGS) -o $@  $(INCLUDE_DIR) $< 

$(EXEC): $(TARGET) $(MEMBERS) $(SOURCEMAIN)
	gcc $(CCOPTIONS) $(PROFFLAGS) -o $(EXEC) $(SOURCEMAIN) $(MEMBERS) $(INCLUDE_DIR) $(LIBDIR) $(LIBS)

info_ukb: $(TARGET) $(MEMBERS) info_ukb.cc
	gcc $(CCOPTIONS) $(PROFFLAGS) -o info_ukb info_ukb.cc $(MEMBERS) $(INCLUDE_DIR) $(LIBDIR) $(LIBS)

create_mcrbin: $(TARGET) $(MEMBERS) create_mcrbin.cc
	gcc $(CCOPTIONS) $(PROFFLAGS) -o create_mcrbin create_mcrbin.cc $(MEMBERS) $(INCLUDE_DIR) $(LIBDIR) $(LIBS)

ukb_sentences: $(TARGET) $(MEMBERS) ukb_sentences.cc
	gcc $(CCOPTIONS) $(PROFFLAGS) -o ukb_sentences ukb_sentences.cc $(MEMBERS) $(INCLUDE_DIR) $(LIBDIR) $(LIBS)

.PHONY : all clean

clean:
	find . -type f -name '*.o' | xargs rm -f
	rm -f $(EXEC)

