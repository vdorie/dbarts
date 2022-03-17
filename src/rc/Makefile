INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR) $(PTHREAD_CPPFLAGS)

LOCAL_SOURCES=bounds.c util.c
LOCAL_OBJECTS=bounds.o util.o

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../rc.a $(LOCAL_OBJECTS)

clean :
	@-rm -f config.h
	@-rm -f *.o
	@-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/rc/bounds.h :
$(INCLUDE_DIR)/rc/util.h : $(INCLUDE_DIR)/misc/stddef.h

bounds.o : bounds.c $(INCLUDE_DIR)/rc/bounds.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c bounds.c -o bounds.o

util.o : util.c $(INCLUDE_DIR)/rc/util.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c util.c -o util.o
