
add_library(malloccountfiles STATIC malloc_count.c malloc_count.h )
set_target_properties(malloccountfiles PROPERTIES LINKER_LANGUAGE CXX)
