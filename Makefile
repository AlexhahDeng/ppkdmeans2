# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Ubu20/ppkdmeans2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Ubu20/ppkdmeans2

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/c/Ubu20/ppkdmeans2/CMakeFiles /mnt/c/Ubu20/ppkdmeans2/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /mnt/c/Ubu20/ppkdmeans2/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named ppkdmeans2

# Build rule for target.
ppkdmeans2: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ppkdmeans2
.PHONY : ppkdmeans2

# fast build rule for target.
ppkdmeans2/fast:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/build
.PHONY : ppkdmeans2/fast

cloud_one.o: cloud_one.cpp.o

.PHONY : cloud_one.o

# target to build an object file
cloud_one.cpp.o:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/cloud_one.cpp.o
.PHONY : cloud_one.cpp.o

cloud_one.i: cloud_one.cpp.i

.PHONY : cloud_one.i

# target to preprocess a source file
cloud_one.cpp.i:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/cloud_one.cpp.i
.PHONY : cloud_one.cpp.i

cloud_one.s: cloud_one.cpp.s

.PHONY : cloud_one.s

# target to generate assembly for a file
cloud_one.cpp.s:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/cloud_one.cpp.s
.PHONY : cloud_one.cpp.s

cloud_two.o: cloud_two.cpp.o

.PHONY : cloud_two.o

# target to build an object file
cloud_two.cpp.o:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/cloud_two.cpp.o
.PHONY : cloud_two.cpp.o

cloud_two.i: cloud_two.cpp.i

.PHONY : cloud_two.i

# target to preprocess a source file
cloud_two.cpp.i:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/cloud_two.cpp.i
.PHONY : cloud_two.cpp.i

cloud_two.s: cloud_two.cpp.s

.PHONY : cloud_two.s

# target to generate assembly for a file
cloud_two.cpp.s:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/cloud_two.cpp.s
.PHONY : cloud_two.cpp.s

comparator.o: comparator.cpp.o

.PHONY : comparator.o

# target to build an object file
comparator.cpp.o:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/comparator.cpp.o
.PHONY : comparator.cpp.o

comparator.i: comparator.cpp.i

.PHONY : comparator.i

# target to preprocess a source file
comparator.cpp.i:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/comparator.cpp.i
.PHONY : comparator.cpp.i

comparator.s: comparator.cpp.s

.PHONY : comparator.s

# target to generate assembly for a file
comparator.cpp.s:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/comparator.cpp.s
.PHONY : comparator.cpp.s

ppkdmeans2.o: ppkdmeans2.cpp.o

.PHONY : ppkdmeans2.o

# target to build an object file
ppkdmeans2.cpp.o:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/ppkdmeans2.cpp.o
.PHONY : ppkdmeans2.cpp.o

ppkdmeans2.i: ppkdmeans2.cpp.i

.PHONY : ppkdmeans2.i

# target to preprocess a source file
ppkdmeans2.cpp.i:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/ppkdmeans2.cpp.i
.PHONY : ppkdmeans2.cpp.i

ppkdmeans2.s: ppkdmeans2.cpp.s

.PHONY : ppkdmeans2.s

# target to generate assembly for a file
ppkdmeans2.cpp.s:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/ppkdmeans2.cpp.s
.PHONY : ppkdmeans2.cpp.s

tools.o: tools.cpp.o

.PHONY : tools.o

# target to build an object file
tools.cpp.o:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/tools.cpp.o
.PHONY : tools.cpp.o

tools.i: tools.cpp.i

.PHONY : tools.i

# target to preprocess a source file
tools.cpp.i:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/tools.cpp.i
.PHONY : tools.cpp.i

tools.s: tools.cpp.s

.PHONY : tools.s

# target to generate assembly for a file
tools.cpp.s:
	$(MAKE) -f CMakeFiles/ppkdmeans2.dir/build.make CMakeFiles/ppkdmeans2.dir/tools.cpp.s
.PHONY : tools.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... ppkdmeans2"
	@echo "... cloud_one.o"
	@echo "... cloud_one.i"
	@echo "... cloud_one.s"
	@echo "... cloud_two.o"
	@echo "... cloud_two.i"
	@echo "... cloud_two.s"
	@echo "... comparator.o"
	@echo "... comparator.i"
	@echo "... comparator.s"
	@echo "... ppkdmeans2.o"
	@echo "... ppkdmeans2.i"
	@echo "... ppkdmeans2.s"
	@echo "... tools.o"
	@echo "... tools.i"
	@echo "... tools.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
