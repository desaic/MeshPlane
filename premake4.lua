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
configuration {"Debug"}
  defines{"Debug"}
  flags{"Symbols"}
              
configuration {"Release"}
  defines{"NDebug"}
  flags{"Optimize"}
              
