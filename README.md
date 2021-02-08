# Triboflow<a name="toc"></a>

Triboflow is a collection of workflows designed to build up crystalline interfaces consisting of nearly arbitrary materials, calculate tribological figures of merit for them using DFT, and store the results in a database.
Triboflow uses the [FireWorks](https://materialsproject.github.io/fireworks/index.html) framework and relies heavily on [pymatgen](https://pymatgen.org/), [atomate](https://atomate.org/), and the [Materialsproject](https://materialsproject.org/) in general.

The project is a collaboration between [Michael Wolloch](https://www.researchgate.net/profile/Michael_Wolloch) at the University of Vienna and the [Group of Prof. M. Clelia Righi](http://www.tribchem.it/) at the University of Bologna.

 1. [Download and Installation](#install)
	  1. [Create a virtual environment for TriboFlow](#virtualenv)
	  2. [ Install and configure MonogoDB locally if you want your database to run also in conda](#mongodb)
	  3. [Create a folder structure, download and install](#triboflowinstall)
	  4. [Configuring FireWorks](#configurefw)
	  5. [Configuring pymatgen](#configurepymatgen)
 2. [Testing your installation of FireWorks](#testing)
 3. ["Fixing" some issues in atomate](#fixingatomate) 
 4. [Using TriboFlow](#using)
	 1. [Running a workflow](#runningawf)
	 2. [Looking at the results](#results)


## Download and Installation <a name="install"></a>
Here we assume that you install Triboflow within a virtual python environment. We use [conda](https://docs.conda.io/en/latest/) in this guide, but [virtualenv](https://virtualenv.pypa.io/en/latest/) or similar is also possible.

### Create a virtual environment for TriboFlow <a name="virtualenv"></a>
 1. Make sure that you have conda installed or download and install miniconda from [here](https://docs.conda.io/en/latest/miniconda.html).
 2. Update conda and then create a new python 3 environment named "TriboFlow" (or whatever you decide) by typing `conda update conda` followed by `conda create --name TriboFlow python=3`
 3. Switch to your new environment: `conda activate TriboFlow`
 4. Now install some fundamental packages with conda: `conda install numpy scipy matplotlib ipython vtk dnspython maggma -c conda-forge`


[Back to top](#toc)
### Install and configure MonogoDB locally if you want your database to run also in conda<a name="mongodb"></a>
This is **optional** and you should skip these steps if you **already have a database running** somewhere that you want to use and can access from the machine you are using. See also [this section](https://atomate.org/installation.html#mongodb) of the atomate installation instructions.
 1. Install the conda package (use conda-forge for a newer version; 4.2.10 instead of 4.0.3) of MongoDB by simply running `conda install -c anaconda mongodb=4.2.10 -c conda-forge` while the TriboFlow environment is active.
 2. You should now have access to the mongo shell via the `mongo`  command and the daemon via the `mongod` command. Before we create the database and set up the daemon, we have to prepare a configuration file and decide where the database should be located. The location will be referred to as `<YourMongoPath>`, and you put the configuration file `mongod.conf` there. It should look somewhat similar to the yaml file below, but you can also tune it according to the options given [here](https://docs.mongodb.com/v4.0/reference/configuration-options/). Importand settings are `storage.dbPath` (where your database will be located), `processManagement.fork` (If true, the process will run in the background and not block your shell), `net.port` (use the standard port for MongoDB: 27017), and `security.authorization` (should be disabled in the beginning to set up the users). Note that you must use spaces and not tabs in YAML!

Example mongod.conf file:<a name="mongodconf"></a>
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

5. Stop the mongo daemon by executing `mongod --shutdown --dbpath <YoutMongoPath>/data/db` and activate authorization in the `mongod.conf` file by changing the last line to `authorization: enabled`. Now, your database is password protected and you have to the username and password you created before in the mongo shell. To do that, while in the mongo shell, switch to the admin database by executing `use admin` and then authenticate with `db.auth(<RootUser>, <RootPassword>)`. Then exit the mongo shell again.

6. It might be a good idea to define some aliases in your `.bashrc` or `.myaliases` files to start and stop the mongod daemon, so you can do it quickly if needed. E.g.:
```
alias mongo_start="mongod -f <YourMongoPath>/mongod.conf"
alias mongo_stop="mongod --shutdown --dbpath <YourMongoPath>/data/db"
```

[Back to top](#toc)
### Create a folder structure, download and install TriboFlow<a name="triboflowinstall"></a>
1. Select a location where you want to have your TriboFlow files located. We will assume this is `<YourPath>`. Now create two subfolders: `<YourPath>/config` for your configuration files and `<YourPath>/pps` for your pseudopotentials.
2. In the same folder we will now download the TriboFlow files from gitlab (using ssh to facilitate automatic login with ssh-keys) by typing `git clone git@gitlab.com:triboteam/TriboFlow.git`. This will (for now) only work if your gitlab account has authorisation!
3. You should now see a folder `<YourPath>/TriboFlow`. `cd` into it and run `pip install -e .` to install TriboFlow and all the other packages that are required to run it into your active conda environment. 
4. If there is an error that the directory is not installable, there might not yet be a setup.py file in the master branch. Checkout a development branch by typing `git checkout development` to switch the branch and run `pip install -e .` again.

[Back to top](#toc)
### Configuring FireWorks<a name="configurefw"></a>
Here we assume that the database is locally installed in the same conda environment, but the procedure is not much different if a cloud service like Atlas is used or if the database is hosted on a different server.

1. Change into your `<YourPath>/config` folder and write a `db.json` file so FireWorks can access your database (If your database is not local, the "host" has to change of course, e.g. for a Atlas DB you might have something like `"mongodb+srv://cluster0-4bevc.mongodb.net"` instead of `"localhost"`):
<a name="dbjson"></a>
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
3. Test this by running the following python script inside your `<YourPath>/config` folder:
`from atomate.vasp.database import VaspCalcDb`
`x = VaspCalcDb.from_db_file("db.json")`
`x.reset()`
`print("SUCCESS")`
All is well if you see "SUCCESS" printed on screen.
4. Now we setup the launchpad (A database where Fireworks are started from), using functionalities of FireWorks by typing `lpad init` and select the same database name, username (`<RootUser>`) and password as when you created the database. Keep the default (None) for the `ssl_ca_file` parameter, but enter `admin` as the `authsource` parameter as suggested by the prompt. You can check and modify the results by looking into `my_launchpad.yaml`. 
5. Test this by typing `lpad -l my_launchpad.yaml reset`
6. Now we will setup the local machine as a FireWorker, which is nothing else than a machine which can pull and execute calculations from the launchpad that are ready to run. It contains information about the database file, the command running VASP on this machine, and other information necessary to run VASP smoothly and efficiently on the local machine. For this we put a `my_fworker.yaml` file into the `<YourPath>/config` folder:
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
7. If the computer or cluster where you are installing TriboFlow has a job scheduler, you have to set up a file called `my_qadapter.yaml` in the `<YourPath>/config` directory. `<SchedulerType>` can be PBS/Torque, SLURM, SGE, or IBM LoadLeveler. `pre_rocket` and `post_rocket` are optional commands to be run in the job script before and after the workflow is executed, this is e.g. useful for loading and unloading modules. You probably will have to activate your conda environment here. The `--timeout` option tells the job to stop pulling new FireWorks from the Launchpad after `<sec>` number of seconds, which should of course be smaller than the walltime. E.g. if you have an allowed walltime of 72 hours, set the `--timeout` option to e.g. 172800 (2 days in seconds).  
It is important to note that this type of rocket_launch will only work if the compute nodes on your cluster can access the MongoDB database, which is a problem for many clusters, since only the login nodes have access to the internet and firewall rules are strict. Possible solutions are described [here](https://materialsproject.github.io/fireworks/offline_tutorial.html). The `my_qadapter.yaml` file might look something like this:
~~~
_fw_name: CommonAdapter  
_fw_q_type: <SchedulerType>
rocket_launch: rlaunch -c <YourPath>/config rapidfire --timeout <sec>
nodes: <NrOfNodes>
walltime: <hh:mm:ss>
queue: null  
account: null  
job_name: null  
pre_rocket: <Custom commands to load modules and conda etc.>
post_rocket: <Custom commands to unload modules etc.>
logdir: <YourPath>/logs
~~~
8.  To tell TriboFlow (more precisely FireWorks) where to find the configuration files you just created, we will write the final configuration file `FW_config.yaml` with the line `CONFIG_FILE_DIR: <YourPath>/config`. We will also set an environment variable to tell FireWorks where this file is. Make sure not to use spaces before or after the ‘=’ sign, since this might lead to problems and type
~~~
conda env config vars set FW_CONFIG_FILE=<YourPath>/config/FW_config.yaml  
~~~
9. You will have to restart your conda environment and afterward verify if this worked by e.g. typing conda env config vars list or simply `echo $FW_CONFIG_FILE.`

[Back to top](#toc)
### Configuring pymatgen<a name="configurepymatgen"></a>
When running VASP calculations FireWorks relies heavily on   [pymatgen](https://pymatgen.org/) and [Custodian](https://materialsproject.github.io/custodian/). Some configuration of pymatgen is required:

1. We assume that you have a folder `<EXTRACTED_VASP_POTCAR>` where you have the VASP pseudopotentials (delivered with the VASP package) already extracted. Type `pmg config -p <EXTRACTED_VASP_POTCAR> <YourPath>/pps` to put these potentials in an order where pymatgen can find them. The final file structure should look something like this (you maybe have to rename the directories in the pps folder):
~~~
pps
├── POT_GGA_PAW_PBE  
│   ├── POTCAR.Ac.gz  
│   ├── POTCAR.Ac_s.gz  
│   ├── POTCAR.Ag.gz  
│   └── ...  
├── POT_GGA_PAW_PW91  
│   ├── POTCAR.Ac.gz  
│   ├── POTCAR.Ac_s.gz  
│   ├── POTCAR.Ag.gz  
│   └── ...  
└── POT_LDA_PAW  
    ├── POTCAR.Ac.gz  
    ├── POTCAR.Ac_s.gz  
    ├── POTCAR.Ag.gz  
    └── ...
~~~
2. Now we have to set a config variable (it will be a file .`pmgrc.yaml` in your home folder) so pymatgen can find the potentials and add your default functional as well (this could also be PBE_54) if you have this potential family and did not rename the folders in the previous step):  
`pmg config --add PMG_VASP_PSP_DIR <YourPath>/pps`
`pmg config --add PMG_DEFAULT_FUNCTIONAL PBE`
3. For integration of the Materials Project REST API, you should register for free at the website [https://materialsproject.org/](https://materialsproject.org/) and get an API Key from your dashboard there. Put it in the configuration file:  
    `pmg config --add PMG_MAPI_KEY <Your_API_Key>`

[Back to top](#toc)
## Testing your installation of FireWorks<a name="testing"></a>
 - Try to load a workflow (for a simple structure relaxation of Si) to the launchpad with `atwf add -l vasp -p wf_structure_optimization -m mp-149` and check that it has been added with `lpad get_wflows`. This should result in something like:
~~~
[  
 {  
  "state": "READY",  
  "name": "Si—1",  
  "created_on": "2020-02-27T14:44:42.634000",  
  "states_list": "REA"  
 },  
]
~~~
 - Navigate to a scratch directory and run the workflow (without a scheduler) with `rlaunch rapidfire`
 - Afterwards you can run `lpad get_wflows` again to see that the state has changed from "READY" to "COMPLETED"
 - It would probably be a good idea to continue with the [tutorial of Atomate](https://atomate.org/running_workflows.html) if you are not familiar with it already, but you can also jump straight into TriboFlow in the next section.

[Back to top](#toc)

## "Fixing" some issues in atomate<a name="fixingatomate"></a>
Atomate supplies a buch of Fireworks and workflows that are used in TriboFlow. However, there are some bugs or maybe incomplete features, as it is quite commong for scientific software. We recommend to check the [corresponding issue on the Triboflow github](https://gitlab.com/triboteam/TriboFlow/-/issues/14) for current "fixes" to atomate. At the time of writing, the only thing to change is to slightly change the `OptimizeFW` and `StaticFW` Fireworks in `atomate.vasp.fireworks.core` to automatically copy the vdw_kernel.bindat of VASP to the execution directory if the vdw parameter in the [vasp input set](https://pymatgen.org/pymatgen.io.vasp.sets.html) passed to the Firework is not `None`. For that you have to add `vdw_kernel_dir=VDW_KERNEL_DIR,` in `def __init__` before `**kwargs` and copy the following lines before `t.append(RunVaspCustodian(...`, in the same manner that it is already done for `ScanOptimizeFW`:
```
# Copy the pre-compiled VdW kernel for VASP, if required
if vasp_input_set.vdw is not None:
    t.append(CopyFiles(from_dir=vdw_kernel_dir))
```

[Back to top](#toc)

## Using TriboFlow<a name="using"></a>
This section is not really a complete user manual and more of a quickstart guide. More documentation is going to follow once the package is getting ready to be released.
### Running a workflow<a name="runningawf"></a>
Main workflows are located in the `triboflow.workflows.main` module. To run a workflow one has to import it from there and pass the input parameters in a dictionary.
### The Heterogeneous_WF
This workflow converges the computational parameters and the lattice parameters of two materials, constructs slabs from them using the supplied Miller indices to define the surface direction and matches those slabs to an interface.
In the (near future), also the PES and PPES will be calculated alongside relevant tribological data.
The workflow is imported from the `triboflow.workflows.main` module as:
```
from triboflow.workflows.main import heterogeneous_wf
```
and an inputs dictionary has to be passed to it. An example can be found in `triboflow.tests.heterogeneous_WF`:
```
from fireworks import LaunchPad
from triboflow.workflows.main import heterogeneous_wf

inputs = {'material_1': {'formula': 'Fe',
                         'miller': '110',
                         'min_vacuum': 25,
                         'min_thickness': 10
                         },
          'material_2': {'formula': 'Cu',
                         'miller': '111',
                         'mp_id': 'mp-30',
                         'min_vacuum': 25,
                         'min_thickness': 10
                         },
          'computational_params':{'functional': 'PBE',
                                  'energy_tolerance': 0.001,
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'False'},
          'interface_params':{'interface_distance': 2.5,
                              'max_area': 500,
                              'r1r2_tol': 0.05
                              }
          }

WF = heterogeneous_wf(inputs)
lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
```
#### Inputs dictionary
The input dictionary has to include the  4 keys, which each has to contain another dictionary with some essential inputs:
 - `material_1  = {'formula': <Chemial formula> 'miller': <Miller indices as Str or list of Int>}`
 - `material_2 = {'formula': <Chemial formula> 'miller': <Miller indices as Str or list of Int>}`
 - `computational_params = {'use_vdw': <True or False>}`
 - `interface_params = {'max_area': <Max crossection of the matched cell>}`
 
It is pretty self explanatory what these input are for: The two materials need to be defined (note that just defining the formula will fetch the structure with the lowest energy that matches the formula! Provide an 'mp_id' key for a clear selection) and need a surface orientation. Van der Waals forces can be taken into account or disregarded, and a maximal tolerated cell cross section area (in Angstrom squared) needs to be given to avoid humongous cells during the interface matching process. All other inputs have default values, or are derived from the other ones if not especially provided.
A full list of possible inputs with types are:
1. material_1 (and material_2):
	- formula (str)
	- miller (list of int, or single str)
	- mp_id (str)
	- min_vacuum (float)
	- min_thickness (float)
2. computational_params:
	- use_vdw (bool, or str)
	- volume_tolerance (float)
	- energy_tolerance (float)
	- use_spin (bool, or str)
	- BM_tolerance (float)
	- functional (str)
3. interface_params:
	- max_area (float)
	- interface_distance (float)
	- max_mismatch (float)
	- max_angle_diff (float)
	- r1r2_tol (float)

[Back to top](#toc)

### Looking at the results<a name="results"></a>
The results of TribolFlow are saved in a separate MongoDB database, which nevertheless is hosted on the same server  (note that the `directoryPerDB: true` line in [`mongod.conf`](#mongodconf) assures that the results are stored in a different folder). This database of results is called "triboflow", in contrast to the "FireWorks" database that is used by FireWorks (see the [db.json](#dbjson) file). The results are stored in different [collections](https://docs.mongodb.com/manual/core/databases-and-collections/#databases), separated for the functional used (PBE, or SCAN) and then split between bulk, slab, and interfaces results. Data in the triboflow database can be queried of course directly from the [mongo shell](https://docs.mongodb.com/manual/mongo/), but it is probably more useful to use the python interface to MongoDB, [pymongo](https://pymongo.readthedocs.io/en/stable/). Some functions to aid with this are provided in the `triboflow.utils` module. To look quickly at single results or get a feel for how the data is structured in the database, it might be beneficial to install a GUI for MongoDB, like [Compass](https://www.mongodb.com/products/compass). (Note that you can use the web-GUI of Atlas when you are using this cloud based solution instead of a local installation of MongoDB.)

Also note that there is a web-GUI provided by FireWorks where you can check out the Workflows and Fireworks in your FireWorks database. Just type `lpad web_gui` for that.

[Back to top](#toc)

