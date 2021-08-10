# Set the compiler to be used
COMPILER = g++-11

# Set the paths and extension type
SRC_PATH = src
BUILD_PATH = build
SRC_PATH_EXECS = mains
SRC_EXT = cpp

# Find all cpp files in the src and main directories
SOURCES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)')
SOURCES_EXECS = $(shell find $(SRC_PATH_EXECS) -name '*.$(SRC_EXT)')

# Set the object file names, with the source directory stripped
# from the path, and the build path prepended in its place
OBJECTS = $(SOURCES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)
OBJECTS_EXECS = $(SOURCES_EXECS:$(SRC_PATH_EXECS)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)

# Names of the executables
PROGRAM_NAMES = $(SOURCES_EXECS:$(SRC_PATH_EXECS)/%.cpp=%)

# Compile flags
COMPILE_FLAGS = -std=c++20 -O3 -Wall -Wextra -march=native -fopenmp
LIB_FLAGS = -ltbb -fopenmp

# Set the directories that need to be included
INCLUDES = -I include/ -I /usr/local/include -isystem /usr/local/Cellar/eigen/3.3.9/include/eigen3 -I /usr/local/Cellar/tbb/2021.2.0/include

# Libraries that need to be included by this project
LIBS = -L /usr/local/Cellar/tbb/2021.2.0/lib

.PHONY: default_target
default_target: release

.PHONY: release
release: dirs
	@$(MAKE) all

.PHONY: dirs
dirs:
	@echo "Creating directories"
	@mkdir -p $(dir $(OBJECTS))
	@mkdir -p $(dir $(OBJECTS_EXECS))

.PHONY: clean
clean:
	@echo "Deleting $(PROGRAM_NAMES) symlink"
	@$(RM) $(PROGRAM_NAMES) 
	@echo "Deleting directories"
	@$(RM) -r $(BUILD_PATH)

.PHONY: all
all: $(BUILD_PATH)/$(PROGRAM_NAMES) 

# Creation of the executable
$(BUILD_PATH)/$(PROGRAM_NAMES): %: $(OBJECTS) $(OBJECTS_EXECS)
	@echo "Linking: $(BUILD_PATH)/$(notdir $@)"
	$(COMPILER) $(COMPILE_FLAGS) $(LIBS) $(LIB_FLAGS) $(BUILD_PATH)/$(notdir $@).o $(OBJECTS) -o $(BUILD_PATH)/$(notdir $@) 

# Source file rules
$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	@echo "Compiling: $< -> $@"
	$(COMPILER) $(COMPILE_FLAGS) $(INCLUDES) -c $< -o $@

# Main file rules
$(BUILD_PATH)/%.o: $(SRC_PATH_EXECS)/%.$(SRC_EXT)
	@echo "Compiling: $< -> $@"
	$(COMPILER) $(COMPILE_FLAGS) $(INCLUDES) -c $< -o $@
