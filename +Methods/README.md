# Methods is a toolbox for running methods for inferring gene regulatory networks (GRNs)


Add the parent directory to the path in MATLAB. This makes `+Methods` a MATLAB toolbox and ready for use.

To test run in MATLAB >(2015a):

    Net = datastruct.Network.fetch()
    Data = datastruct.Dataset.fetch()

Now when we have the data we can run the method through a wrapper:

    method_name = 'lsco';
    [Aest,z] = Methods.(method_name)(Data,'full');

`Aest` will now contain a 3D matrix where each layer `i` in `Aest(:,:,i)` will contain a network along the regularization path detailed in the parameter `z`.

This repository has a few external dependencies for running specific wrappers.
For more details see the [GeneSPIDER](https://bitbucket.org/sonnhammergrni/genespider/README.md).
