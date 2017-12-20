cp python_orocos_kdl/CMakeLists.ros.txt python_orocos_kdl/CMakeLists.txt 
mkdir build
cd build/
cmake ..
make
cd lib/
sudo cp PyKDL.so /usr/lib/python3/dist-packages/

python3 ../../test.py

cd ../..
sudo rm -rf build/
sudo rm python_orocos_kdl/CMakeLists.txt
