# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /cygdrive/c/Users/Coyote/.CLion2017.3/system/cygwin_cmake/bin/cmake.exe

# The command to remove a file.
RM = /cygdrive/c/Users/Coyote/.CLion2017.3/system/cygwin_cmake/bin/cmake.exe -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cygdrive/d/paraProg/project/fuckclion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cygdrive/d/paraProg/project/fuckclion/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/fuckclion.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/fuckclion.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fuckclion.dir/flags.make

CMakeFiles/fuckclion.dir/serialordon.cpp.o: CMakeFiles/fuckclion.dir/flags.make
CMakeFiles/fuckclion.dir/serialordon.cpp.o: ../serialordon.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/d/paraProg/project/fuckclion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fuckclion.dir/serialordon.cpp.o"
	D:/progs/cygwin64/bin/g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fuckclion.dir/serialordon.cpp.o -c /cygdrive/d/paraProg/project/fuckclion/serialordon.cpp

CMakeFiles/fuckclion.dir/serialordon.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fuckclion.dir/serialordon.cpp.i"
	D:/progs/cygwin64/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/d/paraProg/project/fuckclion/serialordon.cpp > CMakeFiles/fuckclion.dir/serialordon.cpp.i

CMakeFiles/fuckclion.dir/serialordon.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fuckclion.dir/serialordon.cpp.s"
	D:/progs/cygwin64/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/d/paraProg/project/fuckclion/serialordon.cpp -o CMakeFiles/fuckclion.dir/serialordon.cpp.s

CMakeFiles/fuckclion.dir/serialordon.cpp.o.requires:

.PHONY : CMakeFiles/fuckclion.dir/serialordon.cpp.o.requires

CMakeFiles/fuckclion.dir/serialordon.cpp.o.provides: CMakeFiles/fuckclion.dir/serialordon.cpp.o.requires
	$(MAKE) -f CMakeFiles/fuckclion.dir/build.make CMakeFiles/fuckclion.dir/serialordon.cpp.o.provides.build
.PHONY : CMakeFiles/fuckclion.dir/serialordon.cpp.o.provides

CMakeFiles/fuckclion.dir/serialordon.cpp.o.provides.build: CMakeFiles/fuckclion.dir/serialordon.cpp.o


CMakeFiles/fuckclion.dir/commonordon.cpp.o: CMakeFiles/fuckclion.dir/flags.make
CMakeFiles/fuckclion.dir/commonordon.cpp.o: ../commonordon.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/d/paraProg/project/fuckclion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/fuckclion.dir/commonordon.cpp.o"
	D:/progs/cygwin64/bin/g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fuckclion.dir/commonordon.cpp.o -c /cygdrive/d/paraProg/project/fuckclion/commonordon.cpp

CMakeFiles/fuckclion.dir/commonordon.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fuckclion.dir/commonordon.cpp.i"
	D:/progs/cygwin64/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/d/paraProg/project/fuckclion/commonordon.cpp > CMakeFiles/fuckclion.dir/commonordon.cpp.i

CMakeFiles/fuckclion.dir/commonordon.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fuckclion.dir/commonordon.cpp.s"
	D:/progs/cygwin64/bin/g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/d/paraProg/project/fuckclion/commonordon.cpp -o CMakeFiles/fuckclion.dir/commonordon.cpp.s

CMakeFiles/fuckclion.dir/commonordon.cpp.o.requires:

.PHONY : CMakeFiles/fuckclion.dir/commonordon.cpp.o.requires

CMakeFiles/fuckclion.dir/commonordon.cpp.o.provides: CMakeFiles/fuckclion.dir/commonordon.cpp.o.requires
	$(MAKE) -f CMakeFiles/fuckclion.dir/build.make CMakeFiles/fuckclion.dir/commonordon.cpp.o.provides.build
.PHONY : CMakeFiles/fuckclion.dir/commonordon.cpp.o.provides

CMakeFiles/fuckclion.dir/commonordon.cpp.o.provides.build: CMakeFiles/fuckclion.dir/commonordon.cpp.o


# Object files for target fuckclion
fuckclion_OBJECTS = \
"CMakeFiles/fuckclion.dir/serialordon.cpp.o" \
"CMakeFiles/fuckclion.dir/commonordon.cpp.o"

# External object files for target fuckclion
fuckclion_EXTERNAL_OBJECTS =

fuckclion.exe: CMakeFiles/fuckclion.dir/serialordon.cpp.o
fuckclion.exe: CMakeFiles/fuckclion.dir/commonordon.cpp.o
fuckclion.exe: CMakeFiles/fuckclion.dir/build.make
fuckclion.exe: CMakeFiles/fuckclion.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cygdrive/d/paraProg/project/fuckclion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable fuckclion.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fuckclion.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fuckclion.dir/build: fuckclion.exe

.PHONY : CMakeFiles/fuckclion.dir/build

CMakeFiles/fuckclion.dir/requires: CMakeFiles/fuckclion.dir/serialordon.cpp.o.requires
CMakeFiles/fuckclion.dir/requires: CMakeFiles/fuckclion.dir/commonordon.cpp.o.requires

.PHONY : CMakeFiles/fuckclion.dir/requires

CMakeFiles/fuckclion.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fuckclion.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fuckclion.dir/clean

CMakeFiles/fuckclion.dir/depend:
	cd /cygdrive/d/paraProg/project/fuckclion/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cygdrive/d/paraProg/project/fuckclion /cygdrive/d/paraProg/project/fuckclion /cygdrive/d/paraProg/project/fuckclion/cmake-build-debug /cygdrive/d/paraProg/project/fuckclion/cmake-build-debug /cygdrive/d/paraProg/project/fuckclion/cmake-build-debug/CMakeFiles/fuckclion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fuckclion.dir/depend

