PROJECT (MeshQuery)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FILE(GLOB LIB_SRC "*.cpp")

INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/include")
ADD_LIBRARY (meshquery  ${LIB_SRC})
SET(LIB_TYPE STATIC)	#set the lib type
