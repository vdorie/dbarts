# Read after src/misc/Makefile by the standalone recipe in this directory's
# Makefile: misc's sources come from src/misc, its objects stay in the build
# directory the sub-make runs in.
vpath %.c $(MISC_SRC)
vpath %.h $(MISC_SRC)

# Every object rebuilds when a configuration header it can read changes: the
# checked-in ones, and configure's own config.h, which a source's quoted
# include finds first in a tree configure has run in.
$(LOCAL_OBJECTS): $(SA_CONFIG)/misc/config.h $(SA_CONFIG)/misc/types.h \
                  $(wildcard $(MISC_SRC)/config.h)
