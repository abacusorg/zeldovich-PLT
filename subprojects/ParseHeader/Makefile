ROOT_DIR := ..
include $(ROOT_DIR)/common.mk

LEX = flex
YACC = bison

# the library name
LIBRARY = libparseheader.a

# these are the object files which will make up the object archive
CXXFILES = stringutil.cc phParser.tab.cc phScanner.cc ParseHeader.cc phDriver.cc HeaderStream.cc
OFILES = ${CXXFILES:.cc=.o}

$(LIBRARY): $(LIBRARY)($(OFILES))

.INTERMEDIATE: $(OFILES)

# the example executable
phexe: tryit.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^

HeaderStream.o: HeaderStream.cc ParseHeader.hh

ParseHeader.o: ParseHeader.cc ParseHeader.hh phDriver.hh

phParser.tab.o: phParser.tab.cc phParser.tab.hh

phScanner.o: phScanner.cc phParser.tab.hh phDriver.hh

phScanner.cc phScanner.hh: phScanner.ll phParser.tab.hh phDriver.hh
	$(LEX) $^

# TODO: when we are ready to require bison >= 3.3, we can remove -Wno-deprecated
phParser.tab.cc phParser.tab.hh: phParser.yy phDriver.hh
	$(YACC) -Wno-deprecated $<

phDriver.o: phDriver.cc phDriver.hh phScanner.hh phParser.tab.hh

ParseHeader.hh: phDriver.hh

.PHONY: clean distclean
clean:
	$(RM) phDriver.o phParser.tab.o \
              phScanner.o phScanner.cc phScanner.hh \
              phParser.tab.cc phParser.tab.hh stack.hh location.hh position.hh \
              ParseHeader.o \
              tryit.o \
              phexe testit \
              libparseheader.a
	$(RM) $(wildcard *.o) $(wildcard *~) $(wildcard *.d)

distclean: clean
