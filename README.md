      _____                   _____ _____ _____ _____  ______ _____
     / ____|                 / ____|  __ \_   _|  __ \|  ____|  __ \
    | |  __  ___ _ __   ___ | (___ | |__) || | | |  | | |__  | |__) |
    | | |_ |/ _ \ '_ \ / _ \ \___ \|  ___/ | | | |  | |  __| |  _  /
    | |__| |  __/ | | |  __/ ____) | |    _| |_| |__| | |____| | \ \
     \_____|\___|_| |_|\___||_____/|_|   |_____|_____/|______|_|  \_\

# Generation and Simulation Package for Informative Data ExploRation #

This is the umbrella repository for the GeneSPIDER toolbox.
This toolbox is comprised of five parts:

* [datastruct](https://bitbucket.org/sonnhammergrni/datastruct) repository containing classes which is used for handling relevant data structure formats.
* [Methods](https://bitbucket.org/sonnhammergrni/methods) repository containing original scripts and wrappers for various inference methods accepting special `datastruct` objects for analysis.
* [analyse](https://bitbucket.org/sonnhammergrni/analyse), for analysing data and inference methods.
* [gsUtilities](https://bitbucket.org/sonnhammergrni/gsutilities), include helper functions.
* [illustrate](https://bitbucket.org/sonnhammergrni/illustrate) toolbox for help with data visualization and exporting to other graphical tools/formats.

### What is this repository for? ###
Motivation:
Inference of gene regulatory networks (GRNs) is a central goal in systems biology.
It is therefore important to evaluate the accuracy of GRN inference methods in the light of network and data properties.
Although several packages are available to model, simulate, and analyse GRN inference, they offer limited control of network topology, system dynamics, experimental design, data properties, and noise characteristics.
Independent control of these properties in simulations is key to drawing conclusions about which inference method to use in a given condition and what performance to expect from it, as well as to obtain properties representative of real biological systems.

Results:
We present a Matlab package GeneSPIDER for generation and analysis of networks and data in a dynamical systems framework with focus on the ability to vary properties.
It supplies essential components that have been missing and wrappers to existing network inference methods in common use.
\gs contains tools for controlling and evaluating network topology (random, small-world, scale-free), stability of linear time-invariant systems, signal to noise ratio (SNR), and network Interampatteness.
Procedures for design of perturbation experiments, bootstrapping, analysis of linear dependence, sample selection, scaling of SNR, and performance evaluation are included.
The ability of \gs to independently control network and data properties in simulations, together with its tools to analyse these properties and the quality of inferred GRNs enables much more informative analysis of GRN inference performance than was previously possible.

Availability and Implementation: Source code freely available for download at https://bitbucket.org/sonnhammergrni/genespider, implemented in Matlab.

Contact: tn@kth.se

Supplementary information: online-only supplementary data available at the journal's web site.

### How do I get set up? ###

* *Alternative 1*: Download genespider-v1 from [here](https://bitbucket.org/sonnhammergrni/genespider/downloads)
  and place in e.g. `~/src/genespider`

* *Alternative 2*: Fetch this repository with the command

        git clone git@bitbucket.org:sonnhammergrni/genespider.git ~/src/genespider

    and run `cd ~/src/genespider`.

    To fetch the complete toolbox run

        git submodule init
        git submodule update


* Now you should have the complete toolbox availible by adding the path `~/src/genespider` to your MATLAB path with the command

        addpath('~/src/genespider')


* If fetched through `git` you can develop and keep track of changes for each submodule by checking out the master branch from each submodule with the command `git checkout master`.

* Dependencies:

    * [MATLAB](https://se.mathworks.com/products/matlab/), version including (2015a) and above is preferred.
    * [git](https://git-scm.com/) for easily keeping up to date with the GeneSPIDER toolbox and for fetching sub-module toolboxes.
    * [Glmnet](https://web.stanford.edu/~hastie/glmnet_matlab/) for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N. used with the wrapper script `Methods.Glmnet`, `Methods.Bolasso`
    * [JSONlab](http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave) for exporting to storage format .json or .ubj.
    * [xml4mat](https://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0) is nessessary for exporting to storage format xml.
    * [CVX](http://cvxr.com/cvx/): Matlab software for disciplined convex programming for some implemented Methods.

* Datasets are available [here](https://bitbucket.org/sonnhammergrni/gs-datasets).
* Networks are available [here](https://bitbucket.org/sonnhammergrni/gs-networks).

* To test that basic functionality after setup as detailed above open matlab and run the commands (require version 2015a at least):

        Net = datastruct.Network.fetch()
        Data = datastruct.Dataset.fetch()

    You should have the default network and a dataset loaded in matlab.

### Who do I talk to? ###

* For questions contact [Torbjörn](https://bitbucket.org/temn/)
* How to cite:

    GeneSPIDER -- generation and simulation package for informative data exploration.  
    Tjärnberg *et. al.*, *manuscript.*, (2015).
