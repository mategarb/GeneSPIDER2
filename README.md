      _____                   _____ _____ _____ _____  ______ _____
     / ____|                 / ____|  __ \_   _|  __ \|  ____|  __ \
    | |  __  ___ _ __   ___ | (___ | |__) || | | |  | | |__  | |__) |
    | | |_ |/ _ \ '_ \ / _ \ \___ \|  ___/ | | | |  | |  __| |  _  /
    | |__| |  __/ | | |  __/ ____) | |    _| |_| |__| | |____| | \ \
     \_____|\___|_| |_|\___||_____/|_|   |_____|_____/|______|_|  \_\
     
<img src="https://sonnhammer-tutorials.bitbucket.io/images/gs_logo_small.png"/>
# Gene regulatory network inference benchmarking with controlled network and data properties #

This is the collection repository for the GeneSPIDER toolbox (Generation and Simulation Package for Informative Data ExploRation).
This toolbox is comprised of five modules:

* [datastruct](https://bitbucket.org/sonnhammergrni/genespider/src/master/+datastruct) containing functionality which is used for handling relevant data structure formats.
* [Methods](https://bitbucket.org/sonnhammergrni/genespider/src/master/+Methods) containing original scripts and wrappers for various inference methods accepting special `datastruct` objects for analysis.
* [analyse](https://bitbucket.org/sonnhammergrni/genespider/src/master/+analyse), for analysing data and inference methods.
* [gsUtilities](https://bitbucket.org/sonnhammergrni/genespider/src/master/+gsUtilities), including helper functions.
* [illustrate](https://bitbucket.org/sonnhammergrni/genespider/src/master/+illustrate) toolbox for helping with data visualization and exporting to other graphical tools/formats.

## Why GeneSPIDER? ##
Inference of gene regulatory networks (GRNs) is a central goal in systems biology.
It is therefore important to evaluate the accuracy of GRN inference methods in the light of network and data properties.
Although several packages are available for modelling, simulate, and analyse GRN inference, they offer limited control of network topology together with system dynamics, experimental design, data properties, and noise characteristics.
Independent control of these properties in simulations is key to drawing conclusions about which inference method to use in a given condition and what performance to expect from it, as well as to obtain properties representative of real biological systems.

## GeneSPIDER website ##

The official GeneSPIDER website includes updated and extended tutorials under [this link.](https://sonnhammer-tutorials.bitbucket.io/gene-spider.html)

### Earlier tutorials ###

Initial tutorials for GeneSPIDER are archived under [this link.](https://bitbucket.org/sonnhammergrni/genespider/src/master/ARCHIVE.md)

## Who do I talk to? ##

* For questions contact [Nordling](mailto:torbjorn.nordling@nordlinglab.org)
* How to cite [bibtex]:

        @article{Tjarnberg2017-GeneSPIDER,
        author ="Tj\"arnberg, Andreas and Morgan, Daniel C. and Studham, Matthew  and Nordling, T\"orbjorn E. M. and Sonnhammer, Erik L. L.",
        title  ="GeneSPIDER - gene regulatory network inference benchmarking with controlled network and data properties",
        journal  ="Mol. BioSyst.",
        year  ="2017",
        pages  ="-",
        publisher  ="The Royal Society of Chemistry",
        doi  ="10.1039/C7MB00058H",
        url  ="http://dx.doi.org/10.1039/C7MB00058H",
        }

