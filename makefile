BIN_FOLDER  = bin
SRC_FOLDER  = src
OBJ_FOLDER  = obj
TEST_FOLDER = test
DOC_FOLDER  = doc

TEST_SCRIPT = $(SRC_FOLDER)/test_all.sh

TEST        = test_2d test_3d test_6d test_es
TEST_EXE    = $(addprefix $(BIN_FOLDER)/,$(addsuffix .exe, $(basename $(TEST))))

OBJ         = math_lib
OBJ_LIB     = $(addprefix $(OBJ_FOLDER)/,$(addsuffix .o, $(basename $(OBJ))))

all: dirs
all: $(OBJ_LIB)
all: $(TEST_EXE)

dirs:
	@mkdir -p $(BIN_FOLDER)
	@mkdir -p $(OBJ_FOLDER)
	@mkdir -p $(TEST_FOLDER)

$(BIN_FOLDER)/%.exe: $(SRC_FOLDER)/%.c
	$(CC) -o $@ $(OBJ_LIB) $< -lm

$(OBJ_FOLDER)/%.o: $(SRC_FOLDER)/%.c $(SRC_FOLDER)/%.h
	$(CC) -c -o $@ $<

doc: Doxyfile $(SRC_FOLDER)/math_lib.h
	mkdir -p $(DOC_FOLDER); \
	doxygen Doxyfile; \
#	cd $(DOC_FOLDER)/latex; \
#	$(MAKE) 

test: $(TEST_EXE) $(OBJ_LIB)
	@./$(BIN_FOLDER)/test_2d.exe > $(TEST_FOLDER)/test_2d.log
	@./$(BIN_FOLDER)/test_3d.exe > $(TEST_FOLDER)/test_3d.log
	@./$(BIN_FOLDER)/test_6d.exe > $(TEST_FOLDER)/test_6d.log
	@./$(BIN_FOLDER)/test_es.exe > $(TEST_FOLDER)/test_es.log
	@echo "Running tests..."
	@$(TEST_SCRIPT)

clean:
	rm -rf $(OBJ_LIB) $(TEST_EXE) $(DOC_FOLDER)

