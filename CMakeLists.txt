cmake_minimum_required(VERSION 3.3)
project (VSORTS)

set (OSG_ROOT "/opt/osg" CACHE PATH "OSG Root folder")
set (OSG_BUILD ${OSG_ROOT}/build/)
set (OSG_LIB ${OSG_BUILD}/lib/)

find_library (OPENTHREADS_LIBRARY NAME libOpenThreads.dylib PATHS ${OSG_LIB} NO_DEFAULT_PATH)
find_library (OSG_LIBRARY NAMES libosg.dylib PATHS ${OSG_LIB} NO_DEFAULT_PATH)
find_library (OSGDB_LIBRARY NAMES libosgDB.dylib PATHS ${OSG_LIB} NO_DEFAULT_PATH)
find_library (OSGGA_LIBRARY NAMES libosgGA.dylib PATHS ${OSG_LIB} NO_DEFAULT_PATH)
find_library (OSGVIEWER_LIBRARY NAMES libosgViewer.dylib PATHS ${OSG_LIB} NO_DEFAULT_PATH)

add_executable (vsort main.cpp)

set_property(TARGET vsort PROPERTY CXX_STANDARD 14)

target_include_directories (vsort PUBLIC ${OSG_ROOT}/include ${OSG_BUILD}/include)
target_link_libraries (
  vsort
  ${OPENTHREADS_LIBRARY}
  ${OSG_LIBRARY}
  ${OSGDB_LIBRARY}
  ${OSGGA_LIBRARY}
  ${OSGVIEWER_LIBRARY}
)
