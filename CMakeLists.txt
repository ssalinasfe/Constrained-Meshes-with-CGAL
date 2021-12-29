cmake_minimum_required(VERSION 3.1)
project(triangulation)
#CGAL_Qt5 is needed for the drawing and CGAL_Core is needed for this special Kernel.
find_package(CGAL REQUIRED COMPONENTS Qt5 Core)
if(CGAL_FOUND AND CGAL_Qt5_FOUND)
  #required to use basic_viewer
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
  #create the executable of the application
  add_executable(triangulation main.cpp)
  #link it with the required CGAL libraries
  target_link_libraries(triangulation CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Core)
else()
  message("ERROR: this program requires CGAL and CGAL_Qt5 and will not be compiled.")
endif()