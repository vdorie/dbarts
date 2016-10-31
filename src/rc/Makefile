INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR)
ALL_CPPFLAGS=$(R_XTRA_CPPFLAGS) $(PKG_CPPFLAGS) $(CPPFLAGS)

LOCAL_SOURCES=bounds.c util.c
LOCAL_OBJECTS=bounds.o util.o

all : librc.a

librc.a : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) librc.a $(LOCAL_OBJECTS)

clean :
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/rc/bounds.h :
$(INCLUDE_DIR)/rc/util.h :

bounds.o : bounds.c $(INCLUDE_DIR)/rc/bounds.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c bounds.c -o bounds.o

util.o : util.c $(INCLUDE_DIR)/rc/util.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c util.c -o util.o
