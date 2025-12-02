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

## To run parallel code locally with mpirun:
- use same convenience script
- run ``` mpirun -n 4 ./numsim "../src/lid_driven_cavity.txt"``` from build/ folder to run with 4 processes

## To run the code on the IPVS cluster:
- to be able to access the ipvslogin.informatik.uni-stuttgart.de server you have to not only be in the uni eduroam network via the normal vpn, but in the informatics network via their vpn.
- To compile the code on the cluster one has to use the module system there to load the required compilers and libraries.
simcl1:[numsim-project] >module load gcc/10.2 openmpi/3.1.6-gcc-10.2 vtk/9.0.1 cmake/3.18.2
simcl1:[numsim-project] >./set_up_and_cmake.sh

### get into the informatics network:
- one can use the normal cisco anyconnect vpn to connect to the informatics network
- but you have to download the informatics network certificate from [http://www.zdi.uni-stuttgart.de/zdi/vpn/anyconnect.html](http://www.zdi.uni-stuttgart.de/zdi/vpn/anyconnect.html) (from within the informatics network for example by accessing the informatics computer pool)
- then install the certificate in windows according to the instructions on the website.
- to be able to access the informatics computer pool and then the vpn in vpn.informatik.uni-stuttgart.de you have to use you informatics login which is different from your ipvs cluster login

### access the ipvs cluster:
- once you are in the informatics network you can ssh into the ipvs cluster via ```ssh -X username@ipvslogin.informatik.uni-stuttgart.de``` as described in https://numsim-exercises.readthedocs.io/en/latest/exercise1/general_remarks.html
- for this ssh you have to use the ipvs login
- to then use the simcl1 computer inside that cluster use ```ssh -X simcl1```
- then use slurm to submit jobs as described in [https://numsim-exercises.readthedocs.io/en/latest/exercise2/cluster.html](https://numsim-exercises.readthedocs.io/en/latest/exercise2/cluster.html) 


## Helpful links:
- [https://numsim-exercises.readthedocs.io/en/latest/exercise1/index.html](https://numsim-exercises.readthedocs.io/en/latest/exercise1/index.html)
- [https://ipvs.informatik.uni-stuttgart.de/SGS/numsim/ex1_doxygen/annotated.html](https://ipvs.informatik.uni-stuttgart.de/SGS/numsim/ex1_doxygen/annotated.html)
- Submission system: https://opendihu-ci.informatik.uni-stuttgart.de/numsim/