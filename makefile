BIN_FOLDER  = bin
SRC_FOLDER  = src
OBJ_FOLDER  = obj
TEST_FOLDER = test
DOC_FOLDERS = doc/html doc/latex

TEST_SCRIPT = scripts/test_all.sh

TEST        = test_2d test_3d test_6d test_es
TEST_EXE    = $(addprefix $(BIN_FOLDER)/,$(addsuffix .exe, $(basename $(TEST))))

OBJ         = libbbmutils
OBJ_LIB     = $(addprefix $(OBJ_FOLDER)/,$(addsuffix .o, $(basename $(OBJ))))

DOXYFILE    = Doxyfile
DOXY_HEADER = doc/doxy_header.tex

all: dirs $(OBJ_LIB) $(TEST_EXE)

dirs:
	@mkdir -p $(BIN_FOLDER)
	@mkdir -p $(OBJ_FOLDER)
	@mkdir -p $(TEST_FOLDER)
	@mkdir -p $(DOC_FOLDERS)

doc: dirs $(DOXY_HEADER) $(SRC_FOLDER)/* README.md $(DOXYFILE)
	doxygen $(DOXYFILE) ;\
	cd doc/latex ;\
	$(MAKE) 

$(DOXY_HEADER): dirs
	[ -f $(DOXY_HEADER) ] || curl -o $(DOXY_HEADER) https://raw.githubusercontent.com/physycom/templates/master/doxy_header.tex

test: dirs $(OBJ_LIB) $(TEST_EXE) 
	@./$(BIN_FOLDER)/test_2d.exe > $(TEST_FOLDER)/test_2d.log
	@./$(BIN_FOLDER)/test_3d.exe > $(TEST_FOLDER)/test_3d.log
	@./$(BIN_FOLDER)/test_6d.exe > $(TEST_FOLDER)/test_6d.log
	@./$(BIN_FOLDER)/test_es.exe > $(TEST_FOLDER)/test_es.log
	@echo "Running tests..."
	@$(TEST_SCRIPT)

$(BIN_FOLDER)/%.exe: $(SRC_FOLDER)/%.c dirs 
	$(CC) -I.. -o $@ $(OBJ_LIB) $< -lm

$(OBJ_FOLDER)/lib%.o: $(SRC_FOLDER)/%.c %.h dirs 
	$(CC) -I.. -c -o $@ $<

clean:
	rm -f $(OBJ_LIB) $(TEST_EXE) $(DOXY_HEADER)

cleanall:
	rm -rf $(BIN_FOLDER) $(OBJ_FOLDER) $(DOC_FOLDERS)
