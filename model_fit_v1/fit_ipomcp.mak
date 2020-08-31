BOOST = "/usr/include/boost" #directory for boost library
LIB64 = "/usr/lib" #64 bit lib
STD64 = "/usr/lib/gcc/x86_64-linux-gnu/7.5.0" #64 bit std
#release variables
FLAGS = -O3 -DNDEBUG -DBOOST_SYSTEM_STATIC_LINK=1 -DBOOST_THREAD_USE_LIB=1 -DDATE_TIME_INLINE
LINKING = -L$(BOOST) -lboost_thread -lboost_random -lboost_system -L$(LIB64) -lpthread -lrt -L$(STD64) -lstdc++

SOURCES = main.cpp  global.cpp
INCLUDES = -I$(BOOST)

EXECUTABLE = fit_ipomcp
OBJ_DIR = objects

.PHONY: all clean
.SILENT:

OBJECTS = $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(basename $(notdir $(SOURCES)))))
DEPFILES = $(OBJECTS:.o=.d)

all: $(EXECUTABLE)

-include $(DEPFILES)

$(EXECUTABLE): $(OBJ_DIR) $(OBJECTS)
	@echo \>\> linking $(EXECUTABLE)
	g++ $(OBJECTS) $(FLAGS) -o $(EXECUTABLE) $(LINKING)


$(OBJ_DIR):
	mkdir $@

$(OBJECTS):
	g++ $(filter $(basename $(@F)).%, $(SOURCES)) $(FLAGS) $(INCLUDES) -MMD -MF $(@:.o=.d) -MT $@ -c -o $@

clean:
	del $(OBJ_DIR)
	del $(EXECUTABLE)
