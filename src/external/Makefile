INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR) $(PTHREAD_CPPFLAGS)

LOCAL_SOURCES=io.c linearAlgebra.c randomBase.c randomNorm.c random.c
LOCAL_OBJECTS=io.o linearAlgebra.o randomBase.o randomNorm.o random.o

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../external.a $(LOCAL_OBJECTS)

clean :
	@-rm -f config.h
	@-rm -f *.o
	@-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/external/io.h :
$(INCLUDE_DIR)/external/linearAlgebra.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/external/R.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/external/random.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/external/Rinternals.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/external/stats.h : $(INCLUDE_DIR)/misc/stddef.h

io.o : io.c $(INCLUDE_DIR)/external/io.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c io.c -o io.o

linearAlgebra.o : linearAlgebra.c $(INCLUDE_DIR)/external/linearAlgebra.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c linearAlgebra.c -o linearAlgebra.o

random.o : random.c $(INCLUDE_DIR)/external/random.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c random.c -o random.o

randomBase.o : randomBase.c $(INCLUDE_DIR)/external/random.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c randomBase.c -o randomBase.o

randomNorm.o : randomNorm.c $(INCLUDE_DIR)/external/random.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c randomNorm.c -o randomNorm.o

