file(GLOB cgal_files
     "*.h"
     "*.c"
)

add_library(cgalsource STATIC ${cgal_files} )
set_target_properties(cgalsource PROPERTIES LINKER_LANGUAGE CXX)
