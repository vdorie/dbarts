BART_INC=../$(INST_DIR)/include/dbarts
HEADERS=-I../$(INST_DIR)/include -I../include

PKG_CPPFLAGS=$(HEADERS)
ALL_CPPFLAGS=$(R_XTRA_CPPFLAGS) $(PKG_CPPFLAGS) $(CPPFLAGS)

LOCAL_SOURCES = bartFit.cpp birthDeathRule.cpp changeRule.cpp functions.cpp \
          likelihood.cpp node.cpp parameterPrior.cpp state.cpp \
          swapRule.cpp tree.cpp treePrior.cpp
LOCAL_OBJECTS = bartFit.o birthDeathRule.o changeRule.o functions.o \
          likelihood.o node.o parameterPrior.o state.o \
          swapRule.o tree.o treePrior.o

all : libdbarts.a

libdbarts.a : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) libdbarts.a $(LOCAL_OBJECTS)

clean :
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(BART_INC)/bartFit.hpp : $(BART_INC)/types.hpp $(BART_INC)/control.hpp $(BART_INC)/data.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp $(BART_INC)/state.hpp
$(BART_INC)/control.hpp :
$(BART_INC)/data.hpp : $(BART_INC)/types.hpp
$(BART_INC)/model.hpp :
$(BART_INC)/scratch.hpp :
$(BART_INC)/results.hpp :
$(BART_INC)/types.hpp :

birthDeathRule.hpp : 
changeRule.hpp : 
functions.hpp : $(BART_INC)/types.hpp
likelihood.hpp :
node.hpp : $(BART_INC)/types.hpp
prior.hpp :
swapRule.hpp : 
tree.hpp : node.hpp

bartFit.o : bartFit.cpp $(BART_INC)/bartFit.hpp $(BART_INC)/results.hpp functions.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

birthDeathRule.o : birthDeathRule.cpp birthDeathRule.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

changeRule.o : changeRule.cpp changeRule.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp $(BART_INC)/types.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

functions.o : functions.cpp functions.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp $(BART_INC)/types.hpp birthDeathRule.hpp changeRule.hpp node.hpp swapRule.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

likelihood.o : likelihood.cpp likelihood.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/state.hpp node.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

node.o : node.cpp node.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/data.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp functions.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

parameterPrior.o : parameterPrior.cpp $(BART_INC)/model.hpp $(BART_INC)/control.hpp node.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

state.o : state.cpp $(BART_INC)/state.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/control.hpp functions.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

swapRule.o : swapRule.cpp swapRule.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp $(BART_INC)/types.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

tree.o : tree.cpp tree.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/data.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp $(BART_INC)/state.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

treePrior.o : treePrior.cpp $(BART_INC)/model.hpp $(BART_INC)/data.hpp $(BART_INC)/scratch.hpp $(BART_INC)/types.hpp functions.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@
