# JQuante ![Java CI with Maven](https://github.com/mjw99/jquante/workflows/Java%20CI%20with%20Maven/badge.svg)

JQuante is a basic Java QM toolkit, written in Java, that is inspired by [PyQuante](http://pyquante.sourceforge.net/).

Having started to port the Python PyQuante code to Java, it was soon discovered that some of the PyQuante code had been ported by V. Ganesh in his [MeTA Studio](https://github.com/tovganesh/metastudio) project. V. Ganesh, re-licensed this code base to a [MIT](https://tldrlegal.com/license/mit-license) license, hence this can be used here.

The code aims to:
* Provide a basic set of Java QM methods.
* Provide associated unit tests.
* Packaged within a Maven framework and available on [Maven Central](http://search.maven.org/).

*Note, this code is in development and is very likely to contain errors.*

Please see [jquante-jmh](https://github.com/mjw99/jquante-jmh) for benchmarks.

# Compiling from Source
## [Ubuntu Noble](http://releases.ubuntu.com/noble/)
This will compile and install the Jquante to the Debian based OS.

1) Install the needed packages.
```bash
    sudo apt-get install git maven openjdk-17-jdk
```
2) Clone, build and install locally
```bash
    git clone https://github.com/mjw99/jquante.git
    cd jquante ; mvn clean install
```
