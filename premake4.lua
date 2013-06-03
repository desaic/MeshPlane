solution "Planar"
configurations{"Debug", "Release"}
project "Planar"
kind "ConsoleApp"
language "C++"
files{
  "src/*.cpp"
}
includedirs{"include"}
libdirs{"lib"}
links{"GLEW", "GLU", "glut", "GL"}
links{"mesh_query", "png", "pthread"}
buildoptions { "-std=c++0x"}
configuration {"Debug"}
targetdir "debug"
defines{"Debug"}
flags{"Symbols"}
              
configuration {"Release"}
targetdir "release"
defines{"NDebug"}
flags{"Optimize"}
              
