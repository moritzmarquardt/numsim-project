## numsim setup (follow docs of numsim lecture)
- Install wsl with ubuntu
- install g++ and gdp: ```sudo apt-get install build-essential gdb```
- ```sudo apt install cmake```
- ```sudo apt install libvtk9-dev```
- ```sudo apt install qtbase5-dev qttools5-dev qttools5-dev-tools qt5-qmake libqt5opengl5-dev```
- ```sudo apt install openjdk-17-jdk```
- server jdk path to .bashrc
- lib path to .profile at the end via sudo
- paraview installation via https://www.paraview.org/download/ 


## To run the code:
- make sure the build/ folder is empty and exists
- run ```cmake ..``` from build/ folder: This should create the Makefile in build/ folder
- run ```make install``` from build/ folder: This compiles the code and creates the ./numsim executable in build/ folder
- run this executable ```./numsim``` from build/ folder:
this should create the outputs in build/out/


## Helpful links:
- [https://numsim-exercises.readthedocs.io/en/latest/exercise1/index.html](https://numsim-exercises.readthedocs.io/en/latest/exercise1/index.html)
- [https://ipvs.informatik.uni-stuttgart.de/SGS/numsim/ex1_doxygen/annotated.html](https://ipvs.informatik.uni-stuttgart.de/SGS/numsim/ex1_doxygen/annotated.html)