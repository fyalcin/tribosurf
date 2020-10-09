# Triboflow

Triboflow is a collection of workflows designed to build up crystalline interfaces consisting of nearly arbitrary materials, calculate tribological figures of merit for them using DFT, and store the results in a database.
Triboflow uses the [FireWorks](https://materialsproject.github.io/fireworks/index.html) framework and relies heavily on [pymatgen](https://pymatgen.org/), [atomate](https://atomate.org/), and the [Materialsproject](https://materialsproject.org/) in general.

The project is a collaboration between [Michael Wolloch](https://www.researchgate.net/profile/Michael_Wolloch) at the University of Vienna and the [Group of Prof. M. Clelia Righi](http://www.tribchem.it/) at the University of Bologna.


## Download and Installation
Here we assume that you install Triboflow within a virtual python environment. We use [conda](https://docs.conda.io/en/latest/) in this guide, but [virtualenv](https://virtualenv.pypa.io/en/latest/) or similar is also possible.

### Create a virtual environment for TriboFlow:
 1. Make sure that you have conda installed or download and install miniconda from [here](https://docs.conda.io/en/latest/miniconda.html).
 2. Update conda and then create a new python 3 environment named "TriboFlow" (or whatever you decide) by typing `conda update conda` followed by `conda create --name TriboFlow python=3`
 3. Switch to your new environment: `conda activate TriboFlow`
### Install and configure MonogoDB locally if you want your database to run also in conda:
This is **optional** and you should skip these steps if you **already have a database running** somewhere that you want to use and can access from the machine you are using. See also [this section](https://atomate.org/installation.html#mongodb) of the atomate installation instructions.
 1. Install the conda package (which is unfortunately a bit out of date) of MongoDB by simply running `conda install -c anaconda mongodb` while the TriboFlow environment is active.
 2. You should now have access to the mongo shell via the `mongo`  command and the daemon via the `mongod` command. Before we create the database and set up the daemon, we have to prepare a configuration file and decide where the database should be located. The location will be referred to as `<YourMongoPath>`, and you put the configuration file `mongod.conf` there. It should look somewhat similar to the yaml file below, but you can also tune it according to the options given [here](https://docs.mongodb.com/v4.0/reference/configuration-options/). Importand settings are `storage.dbPath` (where your database will be located), `processManagement.fork` (If true, the process will run in the background and not block your shell), `net.port` (use the standard port for MongoDB: 27017), and `security.authorization` (should be disabled in the beginning to set up the users). Note that you must use spaces and not tabs in YAML!

Example mongod.conf file:
 ~~~
storage:  
    dbPath: "<YourMongoPath>/data/db"  
    journal:  
        enabled: true  
        commitIntervalMs: 100  
    directoryPerDB: true  
    engine: "wiredTiger"  
    wiredTiger:  
        engineConfig:  
           cacheSizeGB: 4.0  
           journalCompressor: "snappy"  
           directoryForIndexes: false  
systemLog:
    verbosity: 0  
    quiet: false  
    traceAllExceptions: false  
    path: "<YourMongoPath>/mongod.log"  
    logAppend: true  
    logRotate: "reopen"  
    destination: "file"  
    component:  
        accessControl:  
            verbosity: 0  
        command:  
            verbosity: 0
 processManagement:  
    fork: true  
    pidFilePath: "<YourMongoPath>/pidfile"
net:  
    port: 27017  
    bindIp: localhost # Can also use a comma separated list  
    bindIpAll: false # If true can access from everywhere.  
    ipv6: true  
security:  
    authorization: disabled
~~~
3. Now create the folder for the database to be put in by typing:  
    `cd <YourMongoPath>; mkdir data; mkdir data/db` and then start the mongo daemon using the configuration file by executing: `mongod --config mongod.conf` which should result in something ending with : `child process started successfully, parent exiting` and produce a `mongod.log` and a `pidfile`. If you have some problems, check the logfile.

4. Now we will use the mongo shell to create two users for the databases. We start by going in the shell and switching to the admin database (which we will use for authorization):  
`mongo`
`use admin`
Now we will create an adminUser and a readOnlyUser by typing:
`db.createUser({user: "<RootUser>", pwd: "<RootPassword>", roles: [{role: "root", db: "admin" }]})` and `db.createUser({user: "<ReadUser>", pwd: "<ReadPassword>", roles: [{role:"readAnyDatabase", db: "admin" }]})`
Now exit the mongo shell: 
`exit`
