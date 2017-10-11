# JQuante #

JQuante is a basic Java QM toolkit, written in Java, that is inspired by [PyQuante](http://pyquante.sourceforge.net/).

Having started to port the Python PyQuante code to Java, it was soon discovered that some of the PyQuante code had been ported by V. Ganesh in his [Metastudio](https://github.com/tovganesh/metastudio) project. V. Ganesh, re-licensed this code base to a [MIT](https://tldrlegal.com/license/mit-license) license, hence this can be used here.

The code aims to:

* Provide a basic set of Java QM methods.
* Provide associated unit tests.
* Packaged within a Maven framework and available on [Maven Central](http://search.maven.org/).

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

## [Mac OS 10.13](https://www.apple.com/macos/high-sierra/)

1) Install the needed packages and configure Java.
   
    a) Install Java SDK : http://www.oracle.com/technetwork/java/javase/downloads/index.html
    b) Install maven, follow instructions at: https://maven.apache.org/install.html

2) Clone, build and install locally

    git clone https://bitbucket.org/mjw99/jquante.git
    cd jquante ; mvn clean install

## [Windows](https://www.microsoft.com/en-in/windows/)

1) Install the needed packages and configure Java.
   
    a) Install Java SDK : http://www.oracle.com/technetwork/java/javase/downloads/index.html
    b) Install maven, follow instructions at: https://maven.apache.org/install.html
    c) Install git: https://git-scm.com/download/win
    d) Ensure Java path is set to JDK and not JRE

2) Clone, build and install locally

    git clone https://bitbucket.org/mjw99/jquante.git
    cd jquante ; mvn clean install


