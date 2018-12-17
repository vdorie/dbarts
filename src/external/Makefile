INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR) $(PTHREAD_FLAGS)

LOCAL_SOURCES=adaptiveRadixTree.c binaryIO.c blockingThreadManager.c io.c hierarchicalThreadManager.c linearAlgebra.c moments.c randomBase.c randomNorm.c random.c string.c thread.c
LOCAL_OBJECTS=adaptiveRadixTree.o binaryIO.o blockingThreadManager.o io.o hierarchicalThreadManager.o linearAlgebra.o moments.o randomBase.o randomNorm.o random.o string.o thread.o

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../external.a $(LOCAL_OBJECTS)

clean :
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/external/adaptiveRadixTree.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/alloca.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/binaryIO.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/io.h :
$(INCLUDE_DIR)/external/linearAlgebra.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/random.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/stats_mt.h : $(INCLUDE_DIR)/external/stddef.h $(INCLUDE_DIR)/external/thread.h
$(INCLUDE_DIR)/external/stats.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/stddef.h :
$(INCLUDE_DIR)/external/string.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/thread.h : $(INCLUDE_DIR)/external/stddef.h

adaptiveRadixTree.o : adaptiveRadixTree.c $(INCLUDE_DIR)/external/adaptiveRadixTree.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c adaptiveRadixTree.c -o adaptiveRadixTree.o

binaryIO.o : binaryIO.c $(INCLUDE_DIR)/external/binaryIO.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c binaryIO.c -o binaryIO.o

blockingThreadManager.o : blockingThreadManager.c $(INCLUDE_DIR)/external/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c blockingThreadManager.c -o blockingThreadManager.o

io.o : io.c $(INCLUDE_DIR)/external/io.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c io.c -o io.o

hierarchicalThreadManager.o : hierarchicalThreadManager.c $(INCLUDE_DIR)/external/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c hierarchicalThreadManager.c -o hierarchicalThreadManager.o

linearAlgebra.o : linearAlgebra.c $(INCLUDE_DIR)/external/linearAlgebra.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c linearAlgebra.c -o linearAlgebra.o

moments.o : moments.c $(INCLUDE_DIR)/external/stats.h $(INCLUDE_DIR)/external/stats_mt.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c moments.c -o moments.o

randomBase.o : randomBase.c $(INCLUDE_DIR)/external/random.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c randomBase.c -o randomBase.o

randomNorm.o : randomNorm.c $(INCLUDE_DIR)/external/random.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c randomNorm.c -o randomNorm.o

random.o : random.c $(INCLUDE_DIR)/external/random.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c random.c -o random.o

string.o : string.c $(INCLUDE_DIR)/external/string.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c string.c -o string.o

thread.o : thread.c pthread.c $(INCLUDE_DIR)/external/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c thread.c -o thread.o

pthread.c :

