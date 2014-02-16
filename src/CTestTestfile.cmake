# CMake generated Testfile for 
# Source directory: /home/kisame/Work/Kafedra/kafedra-resheto/src
# Build directory: /home/kisame/Work/Kafedra/kafedra-resheto/src
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(MPILaunchQuad "/usr/bin/mpiexec" "-np" "4" "../bin//kafedra-resheto")
