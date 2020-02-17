BOOL,true = true
BOOL,1 = true
BOOL,false = false
BOOL,0 = false
BOOL, = true
OR,true,true  = true
OR,true,false = true
OR,false,true = true
OR,false,false = false

have.sse2.test, = false
have.sse2 = $(BOOL,$(have.sse2.test,$(SSE2_FLAG)))

have.sse4.1.test, = false
have.sse4.1 = $(BOOL,$(have.sse4.1.test,$(SSE4_1_FLAG)))

have.avx.test, = false
have.avx = $(BOOL,$(have.avx.test,$(AVX_FLAG)))

have.avx2.test, = false
have.avx2 = $(BOOL,$(have.avx2.test,$(AVX2_FLAG)))

INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR) $(PTHREAD_CPPFLAGS)

simd.source.sse2.true = linearAlgebra_sse2.c moments_sse2.c partition_sse2.c
simd.source.sse2.false =
simd.source.sse2 = $(simd.source.sse2.$(have.sse2))

simd.source.sse4.1.true = partition_sse4_1.c
simd.source.sse4.1.false =
simd.source.sse4.1 = $(simd.source.sse2) $(simd.source.sse4.1.$(have.sse4.1))

simd.source.avx.true = linearAlgebra_avx.c
simd.source.avx.false =
simd.source.avx = $(simd.source.sse4.1) $(simd.source.avx.$(have.avx))

simd.source.avx2.true = partition_avx2.c
simd.source.avx2.false = 

simd.sources = $(simd.source.avx) $(simd.source.avx2.$(have.avx2))
simd.objects = $(simd.sources:.c=.o)

LOCAL_SOURCES=adaptiveRadixTree.c binaryIO.c blockingThreadManager.c hierarchicalThreadManager.c linearAlgebra.c moments.c partition.c simd.c string.c thread.c $(simd.sources)
LOCAL_OBJECTS=adaptiveRadixTree.o binaryIO.o blockingThreadManager.o hierarchicalThreadManager.o linearAlgebra.o moments.o partition.o simd.o string.o thread.o $(simd.objects)

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../misc.a $(LOCAL_OBJECTS)

clean :
	-rm -f config.h
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/misc/adaptiveRadixTree.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/alloca.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/binaryIO.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/intrinsic.h :
$(INCLUDE_DIR)/misc/linearAlgebra.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/partition.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/simd.h :
$(INCLUDE_DIR)/misc/stats.h : $(INCLUDE_DIR)/misc/stddef.h $(INCLUDE_DIR)/misc/thread.h
$(INCLUDE_DIR)/misc/stddef.h :
$(INCLUDE_DIR)/misc/string.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/thread.h : $(INCLUDE_DIR)/misc/stddef.h

adaptiveRadixTree.o : adaptiveRadixTree.c $(INCLUDE_DIR)/misc/adaptiveRadixTree.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c adaptiveRadixTree.c -o adaptiveRadixTree.o

binaryIO.o : binaryIO.c $(INCLUDE_DIR)/misc/binaryIO.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c binaryIO.c -o binaryIO.o

blockingThreadManager.o : blockingThreadManager.c $(INCLUDE_DIR)/misc/thread.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c blockingThreadManager.c -o blockingThreadManager.o

hierarchicalThreadManager.o : hierarchicalThreadManager.c $(INCLUDE_DIR)/misc/thread.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c hierarchicalThreadManager.c -o hierarchicalThreadManager.o

linearAlgebra.o : linearAlgebra.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c linearAlgebra.c -o linearAlgebra.o

linearAlgebra_sse2.true = $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(SSE2_FLAG) -c linearAlgebra_sse2.c -o linearAlgebra_sse2.o
linearAlgebra_sse2.false = touch linearAlgebra_sse2.o
compile.linearAlgebra_sse2 = $(linearAlgebra_sse2.$(have.sse2))

linearAlgebra_sse2.o : linearAlgebra_sse2.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	echo $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c linearAlgebra_sse2.c -o linearAlgebra_sse2.c
	@$(compile.linearAlgebra_sse2)

linearAlgebra_avx.true = $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(AVX_FLAG) -c linearAlgebra_avx.c -o linearAlgebra_avx.o
linearAlgebra_avx.false = touch linearAlgebra_avx.o
compile.linearAlgebra_avx = $(linearAlgebra_avx.$(have.avx))

linearAlgebra_avx.o : linearAlgebra_avx.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	echo $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c linearAlgebra_avx.c -o linearAlgebra_avx.c
	@$(compile.linearAlgebra_avx)

moments.o : moments.c $(INCLUDE_DIR)/misc/stats.h $(INCLUDE_DIR)/misc/simd.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c moments.c -o moments.o

moments_sse2.true = $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(SSE2_FLAG) -c moments_sse2.c -o moments_sse2.o
moments_sse2.false = touch moments_sse2.o
compile.moments_sse2 = $(moments_sse2.$(have.sse2))

moments_sse2.o : moments_sse2.c $(INCLUDE_DIR)/misc/stats.h
	echo $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c moments_sse2.c -o moments_sse2.o
	@$(compile.moments_sse2)

partition.o : partition.c $(INCLUDE_DIR)/misc/partition.h partition_body.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c partition.c -o partition.o

partition_sse2.true = $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(SSE2_FLAG) -c partition_sse2.c -o partition_sse2.o
partition_sse2.false = touch partition_sse2.o
compile.partition_sse2 = $(partition_sse2.$(have.sse2))

partition_sse2.o : partition_sse2.c $(INCLUDE_DIR)/misc/partition.h partition_body.c
	echo $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c partition_sse2.c -o partition_sse2.o
	@$(compile.partition_sse2)

partition_sse4.1.true = $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(SSE4_1_FLAG) -c partition_sse4_1.c -o partition_sse4_1.o
partition_sse4.1.false = touch partition_sse4_1.o
compile.partition_sse4.1 = $(partition_sse4.1.$(have.sse4.1))

partition_sse4_1.o : partition_sse4_1.c $(INCLUDE_DIR)/misc/partition.h partition_body.c
	echo $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c partition_sse4_1.c -o partition_sse4_1.o
	@$(compile.partition_sse4.1)

partition_avx2.true = $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(AVX2_FLAG) -c partition_avx2.c -o partition_avx2.o
partition_avx2.false = touch partition_avx2.o
compile.partition_avx2 = $(partition_avx2.$(have.avx2))

partition_avx2.o : partition_avx2.c $(INCLUDE_DIR)/misc/partition.h partition_body.c
	echo $(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c partition_avx2.c -o partition_avx2.o
	@$(compile.partition_avx2)

simd.o : simd.c $(INCLUDE_DIR)/misc/simd.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c simd.c -o simd.o

string.o : string.c $(INCLUDE_DIR)/misc/string.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c string.c -o string.o

thread.o : thread.c $(INCLUDE_DIR)/misc/thread.h pthread.h
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c thread.c -o thread.o

partition_body.c : 

