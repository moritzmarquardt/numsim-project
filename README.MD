## numsim setup (follow docs of numsim lecture)
- Install wsl with ubuntu
- install g++ and gdp: ```sudo apt-get install build-essential gdb```
- ```sudo apt install cmake```
- ```sudo apt install libvtk9-dev```
- ```sudo apt install qtbase5-dev qttools5-dev qttools5-dev-tools qt5-qmake libqt5opengl5-dev```
- ```sudo apt install openjdk-17-jdk```
- server jdk path to .bashrc : export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-17-openjdk-amd64/lib:/usr/lib/jvm/java-17-openjdk-amd64/lib/server
- lib path to .profile at the end via sudo nano (same path) -> might be redundant
- paraview installation via https://www.paraview.org/download/ 


## To run the code:
- make sure the build/ folder is empty and exists
- run ```cmake ..``` from build/ folder: This should create the Makefile in build/ folder
- run ```make install``` from build/ folder: This compiles the code and creates the ./numsim executable in build/ folder
- run this executable ```./numsim``` from build/ folder:
this should create the outputs in build/out/
- cmake nur machen, wenn sich die CMakeLists.txt ändert, also wenn man dateien hinzufügt oder entfernt
- wenn das nicht der fall ist einfach schnell machen: ```make install``` und dann ```./numsim``` oder direkt ```make install && ./numsim```


## Helpful links:
- [https://numsim-exercises.readthedocs.io/en/latest/exercise1/index.html](https://numsim-exercises.readthedocs.io/en/latest/exercise1/index.html)
- [https://ipvs.informatik.uni-stuttgart.de/SGS/numsim/ex1_doxygen/annotated.html](https://ipvs.informatik.uni-stuttgart.de/SGS/numsim/ex1_doxygen/annotated.html)
- Submission system: https://opendihu-ci.informatik.uni-stuttgart.de/numsim/