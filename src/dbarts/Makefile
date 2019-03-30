BART_INC=../$(INST_DIR)/include/dbarts

PKG_CPPFLAGS=-I../$(INST_DIR)/include -I../include $(PTHREAD_FLAGS) $(DBARTS_FLAGS)

LOCAL_SOURCES = bartFit.cpp binaryIO.cpp birthDeathRule.cpp changeRule.cpp functions.cpp \
  likelihood.cpp node.cpp parameterPrior.cpp state.cpp \
  swapRule.cpp tree.cpp treePrior.cpp
LOCAL_OBJECTS = bartFit.o binaryIO.o birthDeathRule.o changeRule.o functions.o \
  likelihood.o node.o parameterPrior.o state.o \
  swapRule.o tree.o treePrior.o

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../dbarts.a $(LOCAL_OBJECTS)

clean :
	-rm -f config.hpp
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

binaryIO.hpp : 
birthDeathRule.hpp : 
changeRule.hpp : 
functions.hpp : $(BART_INC)/types.hpp
likelihood.hpp :
node.hpp : $(BART_INC)/types.hpp
prior.hpp :
swapRule.hpp : 
tree.hpp : node.hpp

bartFit.o : bartFit.cpp $(BART_INC)/bartFit.hpp $(BART_INC)/results.hpp binaryIO.hpp functions.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c bartFit.cpp -o bartFit.o

binaryIO.o : binaryIO.cpp binaryIO.hpp $(BART_INC)/control.hpp $(BART_INC)/data.hpp $(BART_INC)/model.hpp $(BART_INC)/state.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c binaryIO.cpp -o binaryIO.o

birthDeathRule.o : birthDeathRule.cpp birthDeathRule.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c birthDeathRule.cpp -o birthDeathRule.o

changeRule.o : changeRule.cpp changeRule.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/types.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c changeRule.cpp -o changeRule.o

functions.o : functions.cpp functions.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/types.hpp birthDeathRule.hpp changeRule.hpp node.hpp swapRule.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c functions.cpp -o functions.o

likelihood.o : likelihood.cpp likelihood.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/state.hpp node.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c likelihood.cpp -o likelihood.o

node.o : node.cpp node.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/data.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp functions.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c node.cpp -o node.o

parameterPrior.o : parameterPrior.cpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/control.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c parameterPrior.cpp -o parameterPrior.o

state.o : state.cpp $(BART_INC)/state.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/control.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c state.cpp -o state.o

swapRule.o : swapRule.cpp swapRule.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/model.hpp $(BART_INC)/types.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c swapRule.cpp -o swapRule.o

tree.o : tree.cpp tree.hpp $(BART_INC)/bartFit.hpp $(BART_INC)/data.hpp $(BART_INC)/model.hpp $(BART_INC)/scratch.hpp $(BART_INC)/state.hpp functions.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c tree.cpp -o tree.o

treePrior.o : treePrior.cpp $(BART_INC)/model.hpp $(BART_INC)/data.hpp $(BART_INC)/types.hpp functions.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c treePrior.cpp -o treePrior.o

