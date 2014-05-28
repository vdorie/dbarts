INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR)
ALL_CPPFLAGS=$(R_XTRA_CPPFLAGS) $(PKG_CPPFLAGS) $(CPPFLAGS)

LOCAL_SOURCES = io.c linearAlgebra.c moments.c stats.c thread.c
LOCAL_OBJECTS = io.o linearAlgebra.o moments.o stats.o thread.o

all : libexternal.a

libexternal.a : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) libexternal.a $(LOCAL_OBJECTS)

clean :
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/external/alloca.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/io.h :
$(INCLUDE_DIR)/external/linearAlgebra.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/stats_mt.h : $(INCLUDE_DIR)/external/stddef.h $(INCLUDE_DIR)/external/thread.h
$(INCLUDE_DIR)/external/stats.h : $(INCLUDE_DIR)/external/stddef.h
$(INCLUDE_DIR)/external/stddef.h :
$(INCLUDE_DIR)/external/thread.h : $(INCLUDE_DIR)/external/stddef.h


linearAlgebra.o : linearAlgebra.c $(INCLUDE_DIR)/external/linearAlgebra.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c $< -o $@

io.o : io.c $(INCLUDE_DIR)/external/io.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c $< -o $@

moments.o : moments.c $(INCLUDE_DIR)/external/stats.h $(INCLUDE_DIR)/external/stats_mt.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c $< -o $@

stats.o : stats.c $(INCLUDE_DIR)/external/stats.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c $< -o $@

thread.o : thread.c pthread.c $(INCLUDE_DIR)/external/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c $< -o $@

pthread.c :
