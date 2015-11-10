# JQuante #

JQuante is a baic Java QM toolkit that was inspired by PyQuante, written in Java. 

Having started off attempting to port the Python PyQuante code to Java, it was soon discovered that some of the PyQuante code had been ported by V. Ganesh in his [Metastudio](https://github.com/tovganesh/metastudio) project. V. Ganesh, re-licensed this code base to a MIT license, hence this can be used here.

The code aims to:

* Provide a basic set of Java QM methods
* Provide associated unit tests
* Packaged within a Maven framework

# Compiling from Source
This will compile and install the Jquante to the Debian based OS.
## [Debian Wheezy](http://www.debian.org/releases/wheezy/) / [Ubuntu Precise](http://releases.ubuntu.com/precise/) / [Ubuntu Trusty](http://releases.ubuntu.com/trusty/)

1) Install the needed packages and configure Java.

    sudo apt-get install git maven openjdk-7-jdk
    # Ensure java7 is selected
    sudo update-alternatives --config java
    sudo update-alternatives --config javac

2) Clone, build and install locally

    git clone https://bitbucket.org/mjw99/jquante.git
    cd jquante ; mvn clean install