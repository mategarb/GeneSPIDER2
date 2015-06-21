
      _____                    _____ _____ _____ _____  ______ _____
     / ____|                  / ____|  __ \_   _|  __ \|  ____|  __ \
    | |  __  ___ _ __   ___  | (___ | |__) || | | |  | | |__  | |__) |
    | | |_ |/ _ \ '_ \ / _ \  \___ \|  ___/ | | | |  | |  __| |  _  /
    | |__| |  __/ | | |  __/  ____) | |    _| |_| |__| | |____| | \ \
     \_____|\___|_| |_|\___| |_____/|_|   |_____|_____/|______|_|  \_\

# Generation and Simulation Package for Informative Data ExploRation #

This is the umbrella repository for the GeneSPIDER toolbox.
The toolbox is comprised of four parts:

* [datastruct](https://bitbucket.org/sonnhammergrni/datastruct) repository containing classes which is used for handling relevant data structure formats.
* [Methods](https://bitbucket.org/sonnhammergrni/methods) repository containing original scripts and wrappers for various inference methods accepting special `datastruct` objects for analysis.
* [tools](https://bitbucket.org/sonnhammergrni/tools), for analysing data and inference methods.
* [illustrate](https://bitbucket.org/sonnhammergrni/illustrate) toolbox for help with visualizing and data and exporting to other graphical tools/formats.

### What is this repository for? ###
A wide range of available tools have been created to be able to simulate, model and analyze gene regulatory networks (GRN). It's not always clear how or when to apply different modeling techniques and approaches to create a useful analysis.  We present a tool that describes and aids in creating a pipeline for how to evaluate and model GRN networks in a dynamical systems framework.  This tool (GeneSpider) not only supply essential components that has previously been missing but also aggregates approaches and sources found in the field that are, or should be, widely used.  GeneSpider is focused on gene regulatory networks from a systems engineering and dynamical systems perspective and will help to bring fourth necessary and useful analysis components both in terms of the data suitable for inferring these kinds of systems, to the systems themselves. The goal is to make it simpler to evaluate data quality and model hypothesis generated from that data.

* Version 0.7
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Fetch this repository with the command

        git clone git@bitbucket.org:sonnhammergrni/genespider.git ~/src/genespider

    and run `cd ~/src/genespider`.

    To fetch the complete toolbox run

        git submodule init
        git submodule update

    Now you will have the complete toolbox availible by adding the path `~/src/genespider` to your MATLAB path with the command

        addpath('~/src/genespider')

* To develop and keep track of changes for each submodule you need to check out the master branch from each submodule with the command `git checkout master`.
* Dependencies:
  * [MATLAB](https://se.mathworks.com/products/matlab/)
  * [git](https://git-scm.com/) for easily keeping up to date with the GeneSpider toolbox.
  * [Glmnet](https://web.stanford.edu/~hastie/glmnet_matlab/) for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N. used with the wrapper script `Methods.Glmnet`, `Methods.Bolasso`
  * [xml4mat](https://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0) is nessessary for exporting to storage format xml.
* Datasets are available [here]()
* How to run tests

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* For questions contact [Andreas](https://bitbucket.org/xparx)
* Other community or team contact
* How to cite: Generation and Simulation Package for Informative Data ExploRation
