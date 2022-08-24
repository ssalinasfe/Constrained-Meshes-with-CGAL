file(GLOB malloc_count_SRC
     "*.h"
     "*.c"
)

add_library(malloccountfiles STATIC ${malloc_count_SRC} )
set_target_properties(malloccountfiles PROPERTIES LINKER_LANGUAGE CXX)
