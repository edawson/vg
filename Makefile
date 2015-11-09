DEP_DIR:=./deps
SRC_DIR:=src
BIN_DIR:=bin
OBJ_DIR:=obj
LIB_DIR:=lib
INC_DIR:=include
CPP_DIR:=cpp

EXE:=vg

CXX:=g++
CXXFLAGS:=-O0 -msse4 -fopenmp -std=c++11

CWD:=$(shell pwd)

LD_INCLUDE_FLAGS:=-I$(INC_DIR) -I. -Isrc
LD_LIB_FLAGS:= -L$(LIB_DIR) -lvcflib -lgssw -lprotobuf -lhts -lpthread -ljansson -lncurses -lrocksdb -lsnappy -lz -lbz2 -lgcsa2 -lxg -lsdsl -ldivsufsort -ldivsufsort64 #-lvg

OBJ:=$(OBJ_DIR)/gssw_aligner.o $(OBJ_DIR)/vg.o cpp/vg.pb.o $(OBJ_DIR)/main.o $(OBJ_DIR)/index.o $(OBJ_DIR)/mapper.o $(OBJ_DIR)/region.o $(OBJ_DIR)/progress_bar.o $(OBJ_DIR)/vg_set.o $(OBJ_DIR)/utility.o $(OBJ_DIR)/path.o $(OBJ_DIR)/alignment.o $(OBJ_DIR)/edit.o $(OBJ_DIR)/sha1.o $(OBJ_DIR)/json2pb.o $(OBJ_DIR)/entropy.o $(OBJ_DIR)/pileup.o $(OBJ_DIR)/caller.o

PROTOBUF_DIR:=deps/protobuf
SDSL_DIR:=deps/sdsl-lite
SNAPPY_DIR:=deps/snappy
ROCKSDB_DIR:=deps/rocksdb
GCSA2_DIR:=deps/gcsa2
PROGRESS_BAR_DIR:=deps/progress_bar
FASTAHACK_DIR:=deps/fastahack
HTSLIB_DIR:=deps/htslib
VCFLIB_DIR:=deps/vcflib
XG_DIR:=deps/xg
GSSW_DIR:=deps/gssw
SPARSEHASH_DIR:=deps/sparsehash
SHA1_DIR:=deps/sha1

$(shell ./source_me.sh)

.PHONY: clean get-deps set-path

all: $(BIN_DIR)/vg $(LIB_DIR)/libvg.a


$(BIN_DIR)/vg: deps $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(LIB_DIR)/libvg.a: $(BIN_DIR)/vg
	ar rs $(LIB_DIR)/libvg.a


deps: $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libgssw.a $(LIB_DIR)/libgcsa2.a $(LIB_DIR)/libsnappy.a $(LIB_DIR)/libvcflib.a sparsehash $(OBJ_DIR)/sha1.o rocksdb htslib $(LIB_DIR)/libxg.a

$(LIB_DIR)/libprotobuf.a:
	cd $(PROTOBUF_DIR) && ./autogen.sh && ./configure --prefix="$(CWD)" && make -j 8 && make install

$(LIB_DIR)/libsdsl.a:
	cd $(SDSL_DIR) && ./install.sh $(CWD)

$(LIB_DIR)/libsnappy.a:
	cd $(SNAPPY_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) && $(MAKE) && $(MAKE) install

rocksdb:
	cd $(ROCKSDB_DIR) && $(MAKE) static_lib && mv librocksdb.a $(CWD)/lib && cp -r include/* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libgcsa2.a: $(LIB_DIR)/libsdsl.a
	cd $(GCSA2_DIR) && $(MAKE) && mv libgcsa2.a $(CWD)/$(LIB_DIR) && cp *.h* $(CWD)/$(INC_DIR)/
	touch $(LIB_DIR)/libgcsa2.a

progress_bar:
	cd $(PROGRESS_BAR_DIR) && $(MAKE) && cp progress_bar.o $(CWD)/$(OBJ_DIR) && cp *.h* $(CWD)/$(INC_DIR)

$(OBJ_DIR)/Fasta.o:
	cd $(FASTAHACK_DIR) && make && mv Fasta.o $(CWD)/$(OBJ_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

htslib:
	cd $(HTSLIB_DIR) && $(MAKE) lib-static && mv libhts.a $(CWD)/$(LIB_DIR) && cp *.h $(CWD)/$(INC_DIR) && cp -r htslib/ $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libxg.a: $(LIB_DIR)/libsdsl.a $(LIB_DIR)/libprotobuf.a
	cd $(XG_DIR) && $(MAKE) all && cp obj/xg.o $(CWD)/$(OBJ_DIR)/ && cp lib/libxg.a $(CWD)/$(LIB_DIR)/ && cp src/*.hpp $(CWD)/$(INC_DIR)/ && cp include/* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libvcflib.a: pre
	cd $(VCFLIB_DIR) && $(MAKE) && cp lib/* $(CWD)/$(LIB_DIR)/ && cp include/* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libgssw.a: pre
	cd $(GSSW_DIR) && $(MAKE) && cp lib/* $(CWD)/$(LIB_DIR)/ && cp obj/* $(CWD)/$(OBJ_DIR) && cp src/*.h $(CWD)/$(INC_DIR)

$(INC_DIR)/lru_cache.h:
	cd $(DEP_DIR)/lru_cache && $(MAKE) && cp *.h* $(CWD)/$(INC_DIR)/

sparsehash:
	cd $(SPARSEHASH_DIR) && ./autogen.sh && ./configure --prefix=$(CWD) && $(MAKE) && $(MAKE) install

$(OBJ_DIR)/sha1.o: $(SHA1_DIR)/sha1.cpp $(SHA1_DIR)/sha1.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) && cp $(SHA1_DIR)/*.h* $(CWD)/$(INC_DIR)/

###################################
## VG source code compilation begins here
####################################

$(CPP_DIR)/vg.pb.cc: $(CPP_DIR)/vg.pb.h pre
$(CPP_DIR)/vg.pb.h: $(LIB_DIR)/libprotobuf.a pre
	protoc $(SRC_DIR)/vg.proto --proto_path=$(SRC_DIR) --cpp_out=cpp
	cp $@ $(INC_DIR)

$(OBJ_DIR)/vg.o: $(SRC_DIR)/vg.cpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libvcflib.a $(FASTAHACK_DIR)/Fasta.o $(LIB_DIR)/libgssw.a sparsehash $(INC_DIR)/lru_cache.h $(INC_DIR)/stream.hpp $(LIB_DIR)/libprotobuf.a $(LIB_DIR)/libsdsl.a progress_bar
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/gssw_aligner.o: $(SRC_DIR)/gssw_aligner.cpp $(SRC_DIR)/gssw_aligner.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libgssw.a $(LIB_DIR)/libprotobuf.a sparsehash $(LIB_DIR)/libvcflib.a
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/vg_set.o: $(SRC_DIR)/vg_set.cpp $(SRC_DIR)/vg_set.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/index.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libgssw.a $(LIB_DIR)/libprotobuf.a sparsehash $(LIB_DIR)/libsdsl.a
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/mapper.o: $(SRC_DIR)/mapper.cpp $(SRC_DIR)/mapper.hpp $(CPP_DIR)/vg.pb.h $(LIB_DIR)/libprotobuf.a sparsehash $(LIB_DIR)/libsdsl.a
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(LIB_DIR)/libvcflib.a $(OBJ_DIR)/Fasta.o $(LIB_DIR)/libgssw.a $(INC_DIR)/stream.hpp $(LIB_DIR)/libprotobuf.a sparsehash $(LIB_DIR)/libsdsl.a rocksdb $(CPP_DIR)/vg.pb.cc 
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/region.o: $(SRC_DIR)/region.cpp $(SRC_DIR)/region.hpp $(LIB_DIR)/libprotobuf.a sparsehash
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)


$(OBJ_DIR)/index.o: $(SRC_DIR)/index.cpp $(SRC_DIR)/index.hpp $(LIB_DIR)/libprotobuf.a sparsehash
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/utility.o: $(SRC_DIR)/utility.cpp $(SRC_DIR)/utility.hpp $(LIB_DIR)/libprotobuf.a sparsehash
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/path.o: $(SRC_DIR)/path.cpp $(SRC_DIR)/path.hpp $(LIB_DIR)/libprotobuf.a sparsehash
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) -Iinclude

$(OBJ_DIR)/edit.o: $(SRC_DIR)/edit.cpp $(SRC_DIR)/edit.hpp $(LIB_DIR)/libprotobuf.a
	$(CXX) $(CXXFLAGS) -c -o edit.o edit.cpp $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/alignment.o: $(SRC_DIR)/alignment.cpp $(SRC_DIR)/alignment.hpp htslib $(LIB_DIR)/libprotobuf.a  sparsehash  $(SRC_DIR)/edit.hpp $(SRC_DIR)/edit.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_INCLUDE_FLAGS)


$(OBJ_DIR)/json2pb.o: $(SRC_DIR)/json2pb.cpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/bin2ascii.h $(LIB_DIR)/libprotobuf.a
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/entropy.o: $(SRC_DIR)/entropy.cpp $(SRC_DIR)/entropy.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/pileup.o: $(SRC_DIR)/pileup.cpp $(SRC_DIR)/pileup.hpp $(CPP_DIR)/vg.pb.h $(INC_DIR)/stream.hpp $(SRC_DIR)/vg.hpp $(SRC_DIR)/json2pb.h $(LIB_DIR)/libprotobuf.a sparsehash
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/caller.o: $(SRC_DIR)/caller.cpp $(SRC_DIR)/caller.hpp $(CPP_DIR)/vg.pb.h $(SRC_DIR)/vg.hpp $(INC_DIR)/stream.hpp $(SRC_DIR)/json2pb.h $(SRC_DIR)/pileup.hpp $(LIB_DIR)/libprotobuf.a sparsehash
	$(CXX) $(CXXFLAGS) -c -o $(OBJ_DIR)/caller.o $(SRC_DIR)/caller.cpp $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

pre:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(CPP_DIR) ]; then mkdir -p $(CPP_DIR); fi

clean:
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r $(CPP_DIR)
	cd $(DEP_DIR) && cd protobuf && $(MAKE) clean
	cd $(DEP_DIR) && cd xg && $(MAKE) clean

## TODO vg source code
## TODO LRU_CACHE
## TODO bash-tap
## TODO sha1
