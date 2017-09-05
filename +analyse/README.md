# analyse is a tool for analysing networks and datasets as well as GRN inference methods performances

Add the parent directory to the path in MATLAB. This makes `+analyse` a MATLAB toolbox and ready for use.

To test in MATLAB >(2015a) (requaires the [+datastruct](https://bitbucket.org/sonnhammergrni/datastruct) and [+Methods](https://bitbucket.org/sonnhammergrni/methods) toolbox) run:

    % prepare data and method result for analysis
    Net = datastruct.Network.fetch()
    Data = datastruct.Dataset.fetch()
    [Aest,z] = Methods.lsco(Data,'full');

After preparing the data and results we can try the analysis

    PN = analyse.Model(Net)
    PD = analyse.Data(Data)
    MR = analyse.CompareModels(Net,Aest)

This should return a brief overview of the network and data and the methods results compared to the gold standard network `Net`.

The measure results, `MR` can be saved by specifying an output file, but we also want to fill in some data specific parameters that would help with analysing the result further in the future

    MR.zetavec = z;
    MR.addprop('method');
    MR.addprop('dataset');
    MR.addprop('network');
    MR.method = repmat({'lsco'},size(z));
    MR.dataset = repmat({Data.dataset},size(z));
    MR.network = repmat({Net.network},size(z));

    save(MR,'test.tsv')

This will save a table of all measures to the file `test.tsv`. dataset and network designators could be any string or number array the length of `z` to identify the specific analysis. If a return value is provided to save, then nothing will be saved and a matlab dataset will be returned to the user.
