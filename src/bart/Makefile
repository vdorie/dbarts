BART_INC=../$(INST_DIR)/include
HEADERS=-I$(BART_INC) -I../include

PKG_CPPFLAGS=$(HEADERS)
ALL_CPPFLAGS=$(R_XTRA_CPPFLAGS) $(PKG_CPPFLAGS) $(CPPFLAGS)

LOCAL_SOURCES = bartFit.cpp birthDeathRule.cpp changeRule.cpp functions.cpp \
          likelihood.cpp node.cpp parameterPrior.cpp \
          swapRule.cpp tree.cpp treePrior.cpp
LOCAL_OBJECTS = bartFit.o birthDeathRule.o changeRule.o functions.o \
          likelihood.o node.o parameterPrior.o \
          swapRule.o tree.o treePrior.o

all : libbart.a

libbart.a : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) libbart.a $(LOCAL_OBJECTS)

clean :
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(BART_INC)/bart/bartFit.hpp : $(BART_INC)/bart/types.hpp $(BART_INC)/bart/control.hpp $(BART_INC)/bart/data.hpp $(BART_INC)/bart/model.hpp
$(BART_INC)/bart/control.hpp :
$(BART_INC)/bart/data.hpp :
$(BART_INC)/bart/model.hpp :
$(BART_INC)/bart/results.hpp :
$(BART_INC)/bart/types.hpp :

birthDeathRule.hpp : 
changeRule.hpp : 
functions.hpp : $(BART_INC)/bart/types.hpp
likelihood.hpp :
node.hpp : $(BART_INC)/bart/types.hpp
prior.hpp :
swapRule.hpp : 
tree.hpp : node.hpp

bartFit.o : bartFit.cpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/results.hpp functions.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

birthDeathRule.o : birthDeathRule.cpp birthDeathRule.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

changeRule.o : changeRule.cpp changeRule.hpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/model.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

functions.o : functions.cpp functions.hpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/model.hpp $(BART_INC)/bart/types.hpp birthDeathRule.hpp changeRule.hpp node.hpp swapRule.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

likelihood.o : likelihood.cpp likelihood.hpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/model.hpp node.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

node.o : node.cpp node.hpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/data.hpp $(BART_INC)/bart/model.hpp functions.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

parameterPrior.o : parameterPrior.cpp $(BART_INC)/bart/model.hpp $(BART_INC)/bart/control.hpp node.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

swapRule.o : swapRule.cpp swapRule.hpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/model.hpp $(BART_INC)/bart/types.hpp functions.hpp likelihood.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

tree.o : tree.cpp tree.hpp $(BART_INC)/bart/bartFit.hpp $(BART_INC)/bart/data.hpp $(BART_INC)/bart/model.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@

treePrior.o : treePrior.cpp $(BART_INC)/bart/model.hpp $(BART_INC)/bart/data.hpp $(BART_INC)/bart/types.hpp functions.hpp node.hpp tree.hpp
	$(CXX) $(ALL_CPPFLAGS) $(CXXFLAGS) -c $< -o $@
