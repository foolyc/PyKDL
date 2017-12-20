As the original PyKDL is usually used in python2 with ros, so i got a naive python3 binding for kdl lib in ros with boost.python.

# Requirement

- python 3.5
- EIGEN3
- KDL in ros kinect or original KDL


# install
before install the PyKDL, ensure that your Eigen3 header file is in /usr/local/include

    sudo sh install.sh


# test

use demo

    python test.py


# tips

- just a basic bindings for KDL use in python3, default parameter/function overload/data type converting have not been processed.
- besides, Joint.None is changed to Joint.NoJoint.
- fk and ik can be used as a function as the module without any solver established in python level.