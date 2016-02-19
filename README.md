DelphesFCCAnalyzer
==================

Software tool to read &amp; analyze FCC-EDM file from FCCSW (Delphes) 

Installation on LxPlus SL machine:
----------------------------------

First, set FCCSW environment (i.e. paths to FCC-EDM and PODIO):
    cd ${FCCSW}
    . init.sh
    cd -

Second, set ROOTSYS variable:
    . setup.sh

At last, build the software using CMake
    mkdir build
    cd build/
    cmake ..
    make install

Installation procedure will install the executable to bin and link it to your ${HOME}/bin directory

To run the code, use help first to see all options:

    analyzeDelphesFCC -h

Options:
--------

    Usage: analyzeDelphesFCC [options]
    Program options:
     -h [ --help ]            Display help
     -i [ --input-file ] arg  Specify name of input FCC-EDM ROOT file
     -o [ --output-file ] arg Specify name of output ROOT file, where results of analysis will be saved)
     -n [ --n-events ] arg    Specify number of processed events, n<=number of events in the input file (optional)
