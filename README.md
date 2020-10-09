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
 4. Now install some fundamental packages with conda: `conda install numpy scipy matplotlib ipython vtk dnspython`

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
### Create a folder structure, download and install TriboFlow
1. Select a location where you want to have your TriboFlow files located. We will assume this is `<YourPath>`. Now create two subfolders: `<YourPath>/config` for your configuration files and `<YourPath>/pps` for your pseudopotentials.
2. In the same folder we will now download the TriboFlow files from gitlab by typing `git clone https://gitlab.com/triboteam/TriboFlow.git`. This will (for now) only work if your gitlab account has authorisation!
3. You should now see a folder `<YourPath>/TriboFlow`. `cd` into it and run `pip install .` to install TriboFlow and all the other packages that are required to run it into your active conda environment. 
3a. If there is an error that the directory is not installable, there might not yet be a setup.py file in the master branch. Checkout a development branch by typing `git checkout development` or `git checkout mwo_dev` to switch branches and run `pip install .` again.

### Configuring everything
Here we assume that the database is locally installed in the same conda environment, but the procedure is not much different if a cloud service like Atlas is used or if the database is hosted on a different server.

1. Change into your `<YourPath>/config` folder and write a `db.json` file so FireWorks can access your database (If your database is not local, the "host" has to change of course, e.g. for a Atlas DB you might have something like `"mongodb+srv://cluster0-4bevc.mongodb.net"` instead of `"localhost"`):
 ~~~
{  
	"host": "localhost",  
	"port": 27017,  
	"database": "FireWorks",  
	"collection": "tasks",  
	"admin_user": "<RootUser>",  
	"admin_password": "<RootPassword>",  
	"readonly_user": "<ReadUser>",  
	"readonly_password": "<ReadPassword>",  
	"aliases": {}  
	"authsource": "admin"  
}
~~~
2. Test this by running the following python script inside your `<YourPath>/config` folder:
`from atomate.vasp.database import VaspCalcDb`
`x = VaspCalcDb.from_db_file("db.json")`
`x.reset()`
`print("SUCCESS")`
All is well if you see "SUCCESS" printed on screen.
3. Now we setup the launchpad (A database where Fireworks are started from), using functionalities of FireWorks by typing `lpad init` and select the same database name, username (`<RootUser>`) and password as when you created the database. Keep the default (None) for the `ssl_ca_file` parameter, but enter `admin` as the `authsource` parameter as suggested by the prompt. You can check and modify the results by looking into `my_launchpad.yaml`. 
4. Test this by typing `lpad -l my_launchpad.yaml reset`
5. Now we will setup the local machine as a FireWorker, which is nothing else than a machine which can pull and execute calculations from the launchpad that are ready to run. It contains information about the database file, the command running VASP on this machine, and other information necessary to run VASP smoothly and efficiently on the local machine. For this we put a `my_fworker.yaml` file into the `<YourPath>/config` folder:
~~~
name: <WorkerName>
category: ''
query: '{}'
env:
    db_file: <YourPath>/config/db.json
    vasp_cmd: mpirun -n <YourCoreCount> <YourVaspCommand>
    scratch_dir: <YourScratchDir>
    vdw_kernel_dir: <YourVdwKernelFolder>
    incar_update:
        KPAR: <YourKparSetting>
        NCORE: <YourNcoreSetting>
~~~

