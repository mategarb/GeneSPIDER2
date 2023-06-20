      _____                   _____ _____ _____ _____  ______ _____
     / ____|                 / ____|  __ \_   _|  __ \|  ____|  __ \
    | |  __  ___ _ __   ___ | (___ | |__) || | | |  | | |__  | |__) |
    | | |_ |/ _ \ '_ \ / _ \ \___ \|  ___/ | | | |  | |  __| |  _  /
    | |__| |  __/ | | |  __/ ____) | |    _| |_| |__| | |____| | \ \
     \_____|\___|_| |_|\___||_____/|_|   |_____|_____/|______|_|  \_\
     
![genespider](https://sonnhammer-tutorials.bitbucket.io/images/gs_logo.png source **=600x400**)
<img src="https://sonnhammer-tutorials.bitbucket.io/images/gs_logo.png" alt="drawing" width="50%" height="50%"/>
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

### Results: ###
We here present the Matlab toolbox GeneSPIDER for generation and analysis of networks and data in a dynamical systems framework with focus on the ability to vary properties on data to mimic plausible real world settings.
It supplies essential components that have been missing to existing network inference methods in common use and wrappers for a selected number of inference methods.
GeneSPIDER contains tools for controlling and evaluating network topology (random, small-world, scale-free), stability of linear time-invariant systems, signal to noise ratio (SNR), and network Interampatteness, properties that has been shown to play a major role in the ability to infer a good network that can explain the data.
Procedures for design of perturbation experiments, bootstrapping, analysis of linear dependence, sample selection, scaling of SNR, and performance evaluation are included.
The ability of GeneSPIDER to independently control network and data properties in simulations, together with its tools to analyse these properties and the quality of inferred GRNs enables much more informative analysis of GRN inference performance than was previously possible.

Availability and Implementation: Source code freely available for download at https://bitbucket.org/sonnhammergrni/genespider, implemented in Matlab.

Contact: [Torbjörn Nordling](mailto:torbjorn.nordling@nordlinglab.org), [Erik sonnhammer](mailto:erik.sonnhammer@scilifelab.se)

Supplementary information: online-only supplementary data available at the journal's web site.

### How do I get set up? ###

* *Alternative 1*: Download genespider-v1.11 from [here](https://bitbucket.org/sonnhammergrni/genespider/downloads)
  and place in e.g. `~/src/genespider`

* *Alternative 2*: Fetch this repository with the command

        git clone https://bitbucket.org/sonnhammergrni/genespider.git ~/src/genespider

    and run `cd ~/src/genespider`.

* Now you should have the complete toolbox availible by adding the path `~/src/genespider` to your MATLAB path with the command

        addpath('~/src/genespider')

* Dependencies:

    * [MATLAB](https://se.mathworks.com/products/matlab/), version including (2015a) and above is preferred.
    * [git](https://git-scm.com/) for easily keeping up to date with the GeneSPIDER toolbox and for fetching sub-module toolboxes.
    * [Glmnet](https://web.stanford.edu/~hastie/glmnet_matlab/) for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N. used with the wrapper script `Methods.Glmnet`, `Methods.Bolasso`
    * [JSONlab](http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave) for exporting to storage format .json or .ubj.
    * [xml4mat](https://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0) is nessessary for exporting to storage format xml.
    * [CVX](http://cvxr.com/cvx/): Matlab software for disciplined convex programming for some implemented Methods.
    * [aracne](http://califano.c2b2.columbia.edu/aracne/): for the aracne wrapper. Configuration files needs to be set as per instructions for aracne defualt or per users specific use cases. The PATH to the aracne home directory needs to be set before MATLAB is started.

* Datasets are available [here](https://bitbucket.org/sonnhammergrni/gs-datasets).
* Networks are available [here](https://bitbucket.org/sonnhammergrni/gs-networks).

* To test that basic functionality after setup as detailed above open matlab and run the commands (require version 2015a at least):

        Net = datastruct.Network.fetch()
        Data = datastruct.Dataset.fetch()

    You should now have the default network and a dataset loaded in matlab.

Generating example data as used in results section
--------------------------------------------------

GeneSPIDER provides four toolboxes: `datastruct`, `analyse`, `Methods`, and
`gsUtilities`. Each toolbox is aimed at a specific function and their
usage is exemplified here. The data used in the examples below can be
downloaded from the online repository at
<https://bitbucket.org/sonnhammergrni/gs-networks>. The network is
[Nordling-D20100302-random-N10-L25-ID1446937](https://bitbucket.org/sonnhammergrni/gs-networks/raw/ec384db2750b5ef229d1c613e3dd04a5e3b634e2/random/N10/Nordling-D20100302-random-N10-L25-ID1446937.json)
and dataset is
[Nordling-ID1446937-D20150825-E15-SNR3291-IDY15968](https://bitbucket.org/sonnhammergrni/gs-datasets/raw/a9d9b00aaa5fa6f4059ba03fd0cb5ec8eb80f3f0/N10/Nordling-ID1446937-D20150825-E15-SNR3291-IDY15968.json).
Note that if the code below is used to generate a new network and
dataset, then they will differ from the presented ones due to the use of
random number generators to create the network and noise matrices. To get the exact same data used here issue the following commands

        Net = datastruct.Network.fetch('Nordling-D20100302-random-N10-L25-ID1446937.json')
        Data = datastruct.Dataset.fetch('Nordling-ID1446937-D20150825-N10-E15-SNR3291-IDY15968.json')

Any named network or data set in the repositories can be input to the fetch function and used as an example.

### Network generation

We start by generating a stable random network with `10` nodes and
sparsity `0.25`. The following code snippet demonstrate how to create a
`datastruct.Network` object with the above specifications.

    N = 10; S=0.25;
    A = datastruct.randomNet(N,S)-eye(N);
    A = datastruct.stabilize(A,'iaa','high');
    Net = datastruct.Network(A,'random');
    setname(Net,struct('creator','Nordling'));
    Net.description = ['This is a sparse network with 10 nodes,'...
        '10 negative self-loops and 15 randomly chosen'...
        'links generated by Nordling 2010-03-02.'...
        'The coefficients are chosen such that they form one'...
        'strong component and a stable dynamical system with'...
        'time constants in the range 0.089 to 12 and an'...
        'interampatteness level of 145 that is in-between'...
        'the estimated level of an E. coli (Gardner et al. 2003 Science)'...
        'and Yeast (Lorenz et al. 2009 PNAS) gene regulatory network.'...
        'The coefficients of the network have not been tuned to explain'...
        'any of the data sets in the mentioned articles.'];

`datastruct.stabilize` takes the random network and the desired IAA as
input parameters and stabilises the network by making the real part of
all eigenvalues negative while adjusting the IAA level. The `setname`
method is used to specify the fields of the `Network` object. The name
is automatically generated based on the network properties to ensure
that each one is unique.

The displayed output of the `Network` object is in this case:

    Net =

      10x10 Network array with properties:

                  network: 'Nordling-D20100302-random-N10-L25-ID1446937'
                        A: [10x10 double]
                        G: [10x10 double]
                    names: {'G1' 'G2' 'G3' 'G4' 'G5' 'G6' 'G7' 'G8' 'G9' 'G10'}
                     desc: 'This is a sparse network with 10 nodes, 10 negative
                            self-loop and 15 randomly chosen links generated by
                            Nordling 2010-03-02. The coefficients are chosen such
                            that they forms one strong component and a stable
                            dynamical system with time constants in the range 0.089
                            to 12 and an interampatteness level of 145 that is in
                            between the estimated level of an
                            E. coli (Gardner et al. 2003 Science) and
                            Yeast (Lorenz et al. 2009 PNAS) gene regulatory network.
                            The coefficients of the network have not been tuned to
                            explain any of the data sets in the mentioned articles.'

The displayed output shows the non-hidden properties of the `Network`
object. `network` is the name of the object, which contains the name of
the creator `Nordling`, the date of creation `D`, the type of network
`random`, the number of nodes, and the number of edges `L`. `A` is the
network matrix. `G` is the static gain matrix (inverse of `A`), which is
precomputed to save time when used in an inference algorithm. `names`
contains the name assigned to each node, which are generated
automatically if they are not specified. `desc` is a description of the
network. The Network class can handle sparse matrices.


### Data generation

We now use the generated network to simulate perturbation experiments to
obtain an expression dataset. The following code snippet simulates `N`
single gene perturbation experiments where each gene is perturbed one by
one followed by `N/2` experiments in which genes are perturbed randomly.

    SNR = 7;
    P = double([eye(N),full(logical(sprandn(N,round(N/2),0.2)))]);
    Y = Net.G*P;
    s = svd(Y);
    stdE = s(N)/(SNR*sqrt(chi2inv(1-analyse.Data.alpha,prod(size(P)))));
    E = stdE*randn(size(P));
    F = zeros(size(P));

We have created a perturbation matrix `P` and a corresponding response
matrix `Y`. The standard deviation has been selected such that the SNR
became 7 when it was used to generate the noise matrix `E`. We didn’t
use the input noise matrix `F` here, but it needs to be specified, so it
was set to zero. With this information, we build a data struct, which we
later use to populate the `Dataset` object.

    D(1).network = Net.network;
    D(1).E = E;
    D(1).F = F;
    D(1).Y = Y+D.E;
    D(1).P = P;
    D(1).lambda = [stdE^2,0];
    D(1).cvY = D.lambda(1)*eye(N);
    D(1).cvP = zeros(N);
    D(1).sdY = stdE*ones(size(D.P));
    D(1).sdP = zeros(size(D.P));

The two easiest ways to populate the `Dataset` object with generated
data is to either initialise it with the data and/or network or to use
the function `populate`. To initialise the `datastruct.Dataset` object
with data we do the following:

    Data = datastruct.Dataset(D,Net);
    setname(Data,struct('creator','Nordling'));
    data.description = ['This data set contains 15 simulated experiments with additive'...
        'white Gaussian noise with variance 0.00028 added to the response'...
        'in order to make the SNR 7 and the data partly informative for'...
        'network inference. The singular values of the response matrix'...
        'are in the range 0.77 to 1.2.'];

The displayed output of the `Dataset` object is in this case:

    Data =

    Dataset with properties:

          dataset: 'Nordling-ID1446937-D20150825-E15-SNR3291-IDY15968'
          network: 'Nordling-D20100302-random-N10-L25-ID1446937'
                P: [10x15 double]
                F: [10x15 double]
              cvP: [10x10 double]
              sdP: [10x15 double]
                Y: [10x15 double]
                E: [10x15 double]
              cvY: [10x10 double]
              sdY: [10x15 double]
           lambda: [0.00028399 0]
            SNR_L: 3.2912
            names: {'G01'  'G02'  'G03'  'G04'  'G05'  'G06'  'G07'  'G08'  'G09'  'G10'}
      description: 'This data set contains 15 simulated experiments with additive
                   white Gaussian noise with variance 0.00028 added to the response
                   in order to make the SNR 7 and the data partly informative for
                   network inference. The singular values of the response matrix
                   are in the range 0.77 to 1.2.'

It is important to be able to connect a dataset to a specific network if
the data was generated , hence the network name is reported in the
`Data` object.

We also provide normalisation procedures for the `Data` object that will normalise the expression matrix `Y`. Three different normalisation procedures are available, standard normalisation, min max range scaling and unit scaling. All methods works over rows or columns, depending on input, *e.g.*

For standard normalisation

    NewData = Data.std_normalize(2);
    sum(response(NewData),2)
    sum(response(NewData).^2,2)

should return zeros as sum over rows and the squared values should be 1 for each sample so the sum over rows should be = M.


For unit scaling

    NewData = Data.unit_length_scaling(2);
    sum(response(NewData).^2,2)

the squared values should sum to 1.

For range scaling

    NewData = Data.range_scaling(2);
    max(response(NewData),[],2)
    min(response(NewData),[],2)
    
the max and min of each row should be 1 and 0 respectively.

It should be noted that the noise estimates are currently not scaled according to the new data and should therefore not be used *as is* in subsequent calculations.

### Analysis

The `analysis` toolbox provides tools to analyse data, networks, and
benchmark results.

First we demonstrate how to load the correct network and dataset from
the online repository:

    v = version('-release');
    if str2num(v(1:end-1)) >= 2015
        disp('Fetching example data online')
        Net = datastruct.Network.fetch('Nordling-D20100302-random-N10-L25-ID1446937.json')
        Data = datastruct.Dataset.fetch('Nordling-ID1446937-D20150825-N10-E15-SNR3291-IDY15968.json')
    else
        disp('Older versions of MATLAB does not support fetching datasets online.')
    end

#### Network analysis:

To analyse the network we input it to the `analyse.Model` module:

    net_prop = analyse.Model(Net);
    disp(net_prop)

It produces the output:

    net_prop =

      Model with properties:

                  network: 'Nordling-D20100302-random-N10-L25-ID1446937'
         interampatteness: 144.6937
        NetworkComponents: 1
            AvgPathLength: 2.8778
                     tauG: 0.085032
                       CC: 0.1
                       DD: 1.5

Six measures are calculated. The interampatteness degree,
`interampatteness`, is the number reported by `cond(A)` in .
`NetworkComponents` is the number of strongly connected components, as
reported by the function `graphconncomp`. `AvgPathLength` is the average
path length of the graph of the network in question, as reported by
`graphallshortestpaths` in . `tauG` is the time constant of the system.
`CC` is the average Clustering coefficient, which can be interpreted as
the neighbourhood sparsity of each node in the network, not considering
the node itself. `DD` is the average degree distribution of the model.
The property `analyse.Model.type` can be set to `directed` (default) or
`undirected` depending on the network and the properties one wishes to
calculate. This is a persistent property, so the value will remain the
default one until it is changed.

Individual properties can also be calculated, all clustering
coefficients can be calculated by

    disp(['Clustering coefficients of the network ',Net.network])
    CCs = analyse.Model.clustering_coefficient(Net)


#### Data analysis:

To analyse the data we input the `Dataset` object to the `analyse.Data`
module:

    data_prop = analyse.Data(Data);
    disp(data_prop)

It will result in the following output:

    data_prop =

    Data with properties:

            dataset: 'Nordling-ID1446937-D20150825-E15-SNR3291-IDY15968'
       SNR_Phi_true: 7
      SNR_Phi_gauss: 3.2912
       SNR_phi_true: 10.991
      SNR_phi_gauss: 10.341

The SNRs reported here correspond the definitions in equations in Tjärnberg *et. al.*, *manuscript.* 2017 by default. However, the SNR is calculated for all `i` with the following two functions:

    disp('SNR estimate based on actual noise matrix E for each variable')
    SNRe = analyse.Data.calc_SNR_phi_true(Data);
    disp(SNRe)

    disp('SNR estimate based on variance estimate each variable')
    SNRl = analyse.Data.calc_SNR_phi_gauss(Data);
    disp(SNRl)

#### Performance evaluation:

To analyse the performance of an inference method we first need to
generate an output. This is accomplished easily thanks to the wrappers.
Each method has an associated wrapper that parses the data of the method
itself. To run the Glmnet implementation we execute:

    [estA,zetavec,zetaRange] = Methods.Glmnet(Data,'full');

The variable `zetavec` is the returned regularisation parameters that
was used within the algorithm. The option “full” will instruct the
method to try to generate the complete regularization path from full to
empty network with the `zeta` values scaled between 0 and 1. It should
be noted that not all methods can reliably do this. For those cases a
zetavec can be specified and supplied to the method. `zetaRange` gives
the scaling factors used for the parameters.

    zetavec = logspace(-6,0,100)
    estA = Methods.Glmnet(Data,zetavec);

and the method will use that vector of values to infere the networks.

To analyse the performance of the model, we input the network estimates
produced by the algorithm to the model comparison method:

    M = analyse.CompareModels(Net,estA);

The `max` operation can now be used to find the optimal performance for
each calculated measure:

    maxM = max(M);

Note that `maxM` will contain the maximum of all measures calculated in
`analyse.CompareModels`. If one wants to get all measures when a
specific measure is maximised, one should specify that as an input.

    max_MCC_M = max(M,'MCC');

This will return all applicable measures to that point.

The measures currently available are detailed in tables in supplementary Tjärnberg *et. al.*, *manuscript.*. `CompareModels` will calculate similarity of non-diagonal elements if the input gold standard
model is not square, assuming that the diagonal has been removed and
truncated along the second dimension.


### Who do I talk to? ###

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

