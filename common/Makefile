CC := gcc -c
CC_FLAGS := -Wall -g -std=gnu11

LD := ar
LD_FLAGS := rcs
     
BIN_DIR := bin
OBJ_DIR := obj
INC_DIR := inc
SRC_DIR := src
 
BIN_FILE := $(BIN_DIR)/common.a
SRC_FILE := $(shell find $(SRC_DIR) -name '*.c')
OBJ_FILE := $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILE))
CC_FLAGS := $(CC_FLAGS) -I$(SRC_DIR) -I$(INC_DIR)
 
.PHONY: all mkdirs clean
 
all: mkdirs $(BIN_FILE)
 
mkdirs:
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)
 
$(BIN_FILE): $(OBJ_FILE)
	$(LD) $(LD_FLAGS) $(BIN_FILE) $(OBJ_FILE)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CC_FLAGS) -c $< -o $@

clean:
	@rm -r $(BIN_DIR)
	@rm -r $(OBJ_DIR)

