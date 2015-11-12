# datastruct is a tool for managing regulatory networks and datasets

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

To test run in MATLAB(2015a):

    Net = datastruct.Network.fetch()
    Data = datastruct.Dataset.fetch()