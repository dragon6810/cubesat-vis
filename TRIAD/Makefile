CC := gcc -c
CC_FLAGS := -Wall -O3 -g -std=gnu11 -fsanitize=address

LD := gcc
LD_FLAGS := -fsanitize=address

BIN_DIR := bin
OBJ_DIR := obj
SRC_DIR := src

BIN_FILE := $(BIN_DIR)/triad
SRC_FILE := $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILE := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILE))
CC_FLAGS := $(CC_FLAGS) -I$(SRC_DIR) -I../common/inc/ -Ilib/inc
LD_FLAGS := $(LD_FLAGS)

LIBS_LINUX := -lm -lGL -lGLU -lGLEW -lglfw
LIBS_MACOS := -lm -framework OpenGL -lGLEW -lglfw
LIBS := ../common/bin/common.a

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
    LIBS += $(LIBS_LINUX)
else ifeq ($(UNAME_S), Darwin)
    LIBS += $(LIBS_MACOS)
endif

.PHONY: all mkdirs clean

all: mkdirs $(BIN_FILE)

mkdirs:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)

$(BIN_FILE): $(OBJ_FILE) ../common/bin/common.a
	$(LD) -o $(BIN_FILE) $(OBJ_FILE) $(LD_FLAGS) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CC_FLAGS) -c $< -o $@

clean:
	@rm -r $(BIN_DIR)
	@rm -r $(OBJ_DIR)
