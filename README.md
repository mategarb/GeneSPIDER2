
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
* [illustrate](https://bitbucket.org/sonnhammergrni/illustrate) toolbox for help with data visualization and exporting to other graphical tools/formats.
* [tools](https://bitbucket.org/sonnhammergrni/tools), include helper functions.

### What is this repository for? ###
  A wide range of tools are available to model, simulate and
  analyze gene regulatory networks (GRN). However, it is not
  always clear how or when to apply different modeling techniques
  to create a meaningful analysis. We present a tool that aids in
  the GRN evaluation and model creation pipeline in a dynamical
  systems framework. This toolbox, (GeneSPIDER) not only supplies
  essential components that have been missing but also aggregates
  approaches and field sources that are, or could be incorporated
  into common use.  The systems engineering perspective of \gs
  focuses on gene regulatory networks to help bring forth the
  necessary and useful analytical components to the dynamical
  systems themselves, both in terms of the data and systems.  The
  goal is to make it simpler to evaluate data quality and model
  hypotheses generated from that data. This package can be fetched from the
  online =git= repository [GeneSPIDER](https://bitbucket.org/sonnhammergrni/genespider).

* Version 0.1
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)
* [Learn Org-mode](http://orgmode.org/)

### How do I get set up? ###

* Fetch this repository with the command

        git clone git@bitbucket.org:sonnhammergrni/genespider.git ~/src/genespider

    or download from `https://bitbucket.org/sonnhammergrni/genespider` to `~/src/genespider`
    and run `cd ~/src/genespider`.

    To fetch the complete toolbox run

        git submodule init
        git submodule update

    Now you will have the complete toolbox availible by adding the path `~/src/genespider` to your MATLAB path with the command

        addpath('~/src/genespider')

* To develop and keep track of changes for each submodule you need to check out the master branch from each submodule with the command `git checkout master`.
* Dependencies:

    * [MATLAB](https://se.mathworks.com/products/matlab/)
    * [git](https://git-scm.com/) for easily keeping up to date with the GeneSPIDER toolbox.
    * [Glmnet](https://web.stanford.edu/~hastie/glmnet_matlab/) for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N. used with the wrapper script `Methods.Glmnet`, `Methods.Bolasso`
    * [JSONlab](http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave) for exporting to storage format .json or .ubj.
    * [xml4mat](https://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0) is nessessary for exporting to storage format xml.
    * [CVX](http://cvxr.com/cvx/): Matlab software for disciplined convex programming for some implemented Methods.

* Datasets are available [here]()

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* For questions contact [Andreas](https://bitbucket.org/xparx)
* Other community or team contact
* How to cite: Generation and Simulation Package for Informative Data ExploRation
