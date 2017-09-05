# datastruct is a toolbox for managing gene regulatory networks and datasets
Or more generally network models and response data/measurements from the network models.


clone the repository
```
git clone git@bitbucket.org:sonnhammergrni/datastruct.git +datastruct
```
or
```
git clone https://$USER@bitbucket.org/sonnhammergrni/datastruct.git +datastruct
```

where `$USER` is your bitbucket user account name.

Add the parent directory to the path in MATLAB. This makes `+datastruct` a MATLAB toolbox and ready for use.

To test run in MATLAB >(2015a):

    Net = datastruct.Network.fetch()
    Data = datastruct.Dataset.fetch()

Datasets and networks availible can be found in the [gs-datasets](https://bitbucket.org/sonnhammergrni/gs-datasets) and [gs-networks](https://bitbucket.org/sonnhammergrni/gs-networks) repositories. For thesting these can be easily loaded with their repository file-names right in the matlab terminal

    Data = datastruct.Dataset.fetch('name-of-dataset.json')
    Net = datastruct.Network.fetch([Data.network,'.json'])

and it will return a network and dataset pair.
