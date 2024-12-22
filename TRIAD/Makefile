CC := gcc -c
CC_FLAGS := -Wall -g -std=gnu11

LD := gcc
LD_FLAGS :=

BIN_DIR := bin
OBJ_DIR := obj
SRC_DIR := src

BIN_FILE := $(BIN_DIR)/triad
SRC_FILE := $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILE := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILE))
CC_FLAGS := $(CC_FLAGS) -I$(SRC_DIR) -I../common/inc/
LD_FLAGS := $(LD_FLAGS) ../common/bin/common.a -lm

.PHONY: all mkdirs clean

all: mkdirs $(BIN_FILE)

mkdirs:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)

$(BIN_FILE): $(OBJ_FILE)
	$(LD) -o $(BIN_FILE) $(OBJ_FILE) $(LD_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CC_FLAGS) -c $< -o $@

clean:
	@rm -r $(BIN_DIR)
	@rm -r $(OBJ_DIR)
