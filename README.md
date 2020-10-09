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
