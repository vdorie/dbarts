INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR) $(PTHREAD_FLAGS)

LOCAL_SOURCES=adaptiveRadixTree.c binaryIO.c blockingThreadManager.c hierarchicalThreadManager.c linearAlgebra.c moments.c partition.c simd.c string.c thread.c
LOCAL_OBJECTS=adaptiveRadixTree.o binaryIO.o blockingThreadManager.o hierarchicalThreadManager.o linearAlgebra.o moments.o partition.o simd.o string.o thread.o

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../misc.a $(LOCAL_OBJECTS)

clean :
	-rm -f config.h
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/misc/adaptiveRadixTree.h : $(INCLUDE_DIR)/misc/stddef.h $(INCLUDE_DIR)/misc/simd.h
$(INCLUDE_DIR)/misc/alloca.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/binaryIO.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/intrinsic.h :
$(INCLUDE_DIR)/misc/linearAlgebra.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/partition.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/simd.h :
$(INCLUDE_DIR)/misc/stats.h : $(INCLUDE_DIR)/misc/stddef.h $(INCLUDE_DIR)/misc/thread.h $(INCLUDE_DIR)/misc/simd.h
$(INCLUDE_DIR)/misc/stddef.h :
$(INCLUDE_DIR)/misc/string.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/thread.h : $(INCLUDE_DIR)/misc/stddef.h

adaptiveRadixTree.o : adaptiveRadixTree.c $(INCLUDE_DIR)/misc/adaptiveRadixTree.h
	$(CC) $(ALL_CPPFLAGS) $(SIMD_FLAGS) $(CFLAGS) -c adaptiveRadixTree.c -o adaptiveRadixTree.o

binaryIO.o : binaryIO.c $(INCLUDE_DIR)/misc/binaryIO.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c binaryIO.c -o binaryIO.o

blockingThreadManager.o : blockingThreadManager.c $(INCLUDE_DIR)/misc/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c blockingThreadManager.c -o blockingThreadManager.o

hierarchicalThreadManager.o : hierarchicalThreadManager.c $(INCLUDE_DIR)/misc/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c hierarchicalThreadManager.c -o hierarchicalThreadManager.o

linearAlgebra.o : linearAlgebra.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	$(CC) $(ALL_CPPFLAGS) $(SIMD_FLAGS) $(CFLAGS) -c linearAlgebra.c -o linearAlgebra.o

moments.o : moments.c $(INCLUDE_DIR)/misc/stats.h
	$(CC) $(ALL_CPPFLAGS) $(SIMD_FLAGS) $(CFLAGS) -c moments.c -o moments.o

partition.o : partition.c $(INCLUDE_DIR)/misc/partition.h partition_body.c
	$(CC) $(ALL_CPPFLAGS) $(SIMD_FLAGS) $(CFLAGS) -c partition.c -o partition.o

simd.o : simd.c $(INCLUDE_DIR)/misc/simd.h
	$(CC) $(ALL_CPPFLAGS) $(SIMD_FLAGS) $(CFLAGS) -c simd.c -o simd.o

string.o : string.c $(INCLUDE_DIR)/misc/string.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c string.c -o string.o

thread.o : thread.c $(INCLUDE_DIR)/misc/thread.h pthread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c thread.c -o thread.o

partition_impl.c : 

