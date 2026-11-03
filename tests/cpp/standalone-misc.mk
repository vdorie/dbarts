# Read ahead of src/misc/Makefile by the standalone recipe in this directory's
# Makefile: misc's sources come from src/misc, its objects stay in the build
# directory the sub-make runs in.
vpath %.c $(MISC_SRC)
vpath %.h $(MISC_SRC)
