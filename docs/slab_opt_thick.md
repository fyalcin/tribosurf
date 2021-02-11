# Convergece slab thickness

## Introduction

**Goal**
This is a workflow to calculate the optimal thickness for a slab of a given material with a specific orientation.

**Context**
DFT simulations usually take a long time to complete, so it is essential to optimize the structural and computational parameters of the calculations. In the case of a slab it is fundamental to provide a correct estimate of its minimum thickness in order to avoid huge slowdown in vasp simulations and not to waste CPU time.
 To estimate it we propose to converge the surface energy or the lattice parameter with respect to the number of atomic layers. 

**Implementation**
At the moment, our code implements the convergence of a slab thickness via the evaluation of its surface energy. This is supposed to work properly only with bulks with mono-atomic basis presenting symmetric slab terminations, i.e. materials from the periodic tables.

## Code Implementatioon

It is possible to use our code both as a single workflow, which scan for the bulk which are saved in a given database and converge the slab thickness for a given miller index, or as a subworkflow called by other workflow as intermediate step.

### Flowchart of the workflows

A schematic representation of our subworkflow is reported in the figure above. Now, a description of the main Firetasks as well as the list of required and optional parameters is made:

![Flowchart](file:///home/gl/Work/WORKFLOW/triboflow/docs/slab_opt_thick.png)

#### Firetasks

1. **FT_SlabOptThick**: Firetask to start a subworkflow as a detour within another workflow. Check if a specific slab (identified by `mp_id`,` functional` and `miller`) is already present in the ` high_level` db and if it contains a <span style =" color: orange ">'opt_thickness'</span> type his dictionary. Otherwise, the corresponding bulk structure is retrieved from the `high_level` db and a subworkflow is started to converge the slab thickness, either via the surface energy or the lattice parameter (not implemented yet). The subworkflow is created through the static method of an external class, e.g. *SlabWF.conv_slabthick_surfene* and is composed of two Firetasks: *FT_StartThickConvo* and *FT_EndThickConvo*.<br/>
Input arguments are:
	- **required_params** = [<span style="color:orange">'mp_id'</span>, <span style="color:orange">'miller'</span>, <span style="color:orange">'functional'</span>]
	- **optional_params** = [<span style="color:orange">'db_file'</span>, <span style="color:orange">'low_level'</span>, <span style="color:orange">'high_level'</span>, <span style="color:orange">'relax_type'</span>, <span style="color:orange">'convo_kind'</span>, <span style="color:orange">'thick_start'</span>, <span style="color:orange">'thick_incr'</span>, <span style="color:orange">'nsteps'</span>, <span style="color:orange">'vacuum'</span>, <span style="color:orange">'bulk_name'</span>, <span style="color:orange">'slab_name'</span>] <br/>
Parameters description:
-- `mp_id`: mp-id of the structure from the MP database.
-- `miller`: miller indexes of the slab orientation.
-- `functional`: pseudopotential functional.
-- `db_file`: path to the location of the database.
-- `low_level`: collection name of the "low level" database. The intermediate calculations and raw data during the convergence process will be saved here.
-- `high_level`: collection name of the "high level" database. The final results concerning the slab optimal thickness and the surface energy will be saved here.
-- `relax_type`: type of relaxation to be performed during the simulation.
-- `convo_kind`: type of convolution. Allowed values are: <span style="color:orange">'surfene'</span>, <span style="color:orange">'alat'</span>.
-- `thick_start`: number of atomic layers for the slab to be used as starting point.
-- `thick_incr`: incremental number of atomic layers at each step.
-- `nsteps`: number of steps to be done to converge the slab thickness.
-- `vacuum`: Minimum vacuum to be used for creating the slabs.
-- `bulk_name`: Name of the bulk dictionary in the `high_level` db containing the data.
-- `slab_name`: Name of the slab dictionary to be saved in the `high_level` db at the end of the convergence.

2. **FT_StartThickConvo**: Firetask to start a subworkflow to calculate and converge the surface energy for a `structure` passed by input, which identified by means of `mp_id`, `miller`, and `functional`. A detour is started to run a self consistent calculation to relax the bulk along a particular orientation and the corresponding slabs, with incresing number of layers. The workflow is created through a static method of an external class, i.e. *SurfEneWF.surface_energy*, and is composed of three Firetasks: *FT_GenerateSlabs*, *FT_RelaxStructure*, and *FT_SurfaceEnergy* <br/>
Input arguments are:
	- **required_params** = [<span style="color:orange">'structure'</span>, <span style="color:orange">'mp_id'</span>, <span style="color:orange">'miller'</span>, <span style="color:orange">'functional'</span>]
	- **optional_params** = [<span style="color:orange">'db_file'</span>, <span style="color:orange">'collection'</span>, <span style="color:orange">'comp_params'</span>, <span style="color:orange">'convo_kind'</span>, <span style="color:orange">'relax_type'</span>, <span style="color:orange">'thick_start'</span>, <span style="color:orange">'thick_incr'</span>, <span style="color:orange">'nsteps'</span>,  <span style="color:orange">'vacuum'</span>, <span style="color:orange">'ext_index'</span>, <span style="color:orange">'slab_name'</span>, <span style="color:orange">'cluster_params'</span>] <br/>
	Parameters description:
-- `structure`: pymatgen structure of the bulk to be converged.
-- `mp_id`: mp-id of the structure from the MP database.
-- `miller`: miller indexes of the slab orientation.
-- `functional`: pseudopotential functional.
-- `db_file`: path to the location of the database.
-- `collection`: collection name of the database where the surface energy data will be stored.
-- `comp_params`: computational params to be used in VASP simulations.
-- `convo_kind`: type of convolution. Allowed values are: <span style="color:orange">'surfene'</span>, <span style="color:orange">'alat'</span>.
-- `relax_type`: type of relaxation to be performed during the simulation.
-- `thick_start`: number of atomic layers for the slab to be used as starting point.
-- `thick_incr`: incremental number of atomic layers at each step.
-- `nsteps`: number of steps to be done to converge the slab thickness.
-- `vacuum`: minimum vacuum to be used for creating the slabs.
-- `ext_index`: boolean exit index, if it is True return the first element from `SlabGenerator.get_slabs`.
-- `slab_name`: name of the slab dictionary to be saved in the `collection` at the end of the convergence.

3. **FT_GenerateSlabs**: Generate a slab or a list of slabs out of a given structure. Parameters that are taken into account to generate the possible different slabs are: `miller`, `thickness`, `vacuum`, `slab_name`. The slabs are generated with SlabGenerator and stored in the `collection` of the database. <br/>
Input arguments are:
	- **required_params** = [<span style="color:orange">'structure'</span>, <span style="color:orange">'mp_id'</span>, <span style="color:orange">'miller'</span>, <span style="color:orange">'functional'</span>]
	- **optional_params** = [<span style="color:orange">'db_file'</span>, <span style="color:orange">'collection'</span>, <span style="color:orange">'thickness'</span>,  <span style="color:orange">'vacuum'</span>, <span style="color:orange">'ext_index'</span>, <span style="color:orange">'symmetrize'</span>, <span style="color:orange">'slab_name'</span>] <br/>
	Parameters description:
-- `structure`: pymatgen structure of the bulk to be converged.
-- `mp_id`: mp-id of the structure from the MP database.
-- `miller`: miller indexes of the slab orientation. It can be just one (list) or more than one (list of lists).
-- `functional`: pseudopotential functional.
-- `db_file`: path to the location of the database.
-- `collection`: collection name of the database where the slabs will be stored.
-- `thickness`: thicknesses to be used to create the slab(s). It can be a single value or a list of values.
-- `vacuum`: minimum vacuum to be used for creating the slabs. It can be either a float or a list.
-- `ext_index`: boolean exit index, if it is True return the first element from `SlabGenerator.get_slabs`.
-- `slab_name`: name of the slabs dictionary to be saved in `collection` at the end of the convergence.

4. **FT_RelaxStructure**: Do a self consistent relaxation of a structure which is passed as input. The structure can be a bulk, a slab or an interfaces. The simulation is performed only if the calculation is not present in `collection` under the the relaxation it checks if the calculation has been already done previously calculated and stored on databases. It stores the final data in a dictionary on the given `collection` under a `tag` created within the Firetask. (It would be necessary to create the tag outside and pass it as required params, in order to have it in higher level FT and retrieve them in next FTs). <br/>
Input arguments are:
	- **required_params** = [<span style="color:orange">'mp_id'</span>, <span style="color:orange">'functional'</span>, <span style="color:orange">'struct_kind'</span>]
	- **optional_params** = [<span style="color:orange">'comp_params'</span>, <span style="color:orange">'miller'</span>, <span style="color:orange">'name'</span>, <span style="color:orange">'db_file'</span>, <span style="color:orange">'collection'</span>, <span style="color:orange">'relax_type'</span>] <br/>
	Parameters description:
-- `mp_id`: mp-id of the structure from the MP database.
-- `functional`: pseudopotential functional.
-- `struct_kind`: Type of structure to be relaxed. Allowed values: <span style="color:orange">'bulk'</span>, <span style="color:orange">'slab'</span>, <span style="color:orange">'interface'</span>.
-- `comp_params`: computational params to be used in VASP simulations.
-- `miller`: miller indexes of the slab orientation. Meaningful only when relaxing slabs and interfaces.
-- `name`: custom name used as dictionary key to search for the structures in `collection`. Useful if bulk and slabs have been saved with different names, really meaninful for interfaces.
-- `db_file`: path to the location of the database.
-- `collection`: collection name of the database where the slabs will be stored.
-- `relax_type`: type of relaxation to be performed during the simulation.

5. **FT_SurfaceEnergy**: Calculate the surface energy of a given number of slabs, compared to their bulk. A list of `tags` pointing to the bulk (first element) and the slabs should be passed as input. The energy is retrieved from the concluded calculations and stored in the provided collection. It ends the calculation of the surface energy subworkflow. <br/>

6. **FT_EndThickConvo**: It finds out the optimal surface energy and the corresponding number of atomic layers, by fitting the data retrieved from the provided collection. The final data (optimal thickness, surface energy, and relaxed slab if done) will be stored in the high_level collection provided as input argument. If the convergence is not reached, it calls the FT_StartThickConvo again with different values for the thickness to be analyzed. Otherwise it ends the optimization of the slab thickness subworkflow. <br/>

7. **(OPTIONAL)** **FT_PutStructInDB**: Retrieve relaxed structure(s) or data from `collection_1` of `db_file_1` and store it to `collection_2` of `db_file_2`. This is meant to be a generalization of the already existing *FT_MakeSlabInDB*. The functionality of this Firetask will be already included within  *FT_SurfaceEnergy*. <br/>
Input arguments are:
	- **required_params** = [<span style="color:orange">'mp_id'</span>, <span style="color:orange">'functional'</span>, <span style="color:orange">'tag'</span>, <span style="color:orange">'struct_kind'</span>]
	- **optional_params** = [<span style="color:orange">'miller'</span>, <span style="color:orange">'name'</span>, <span style="color:orange">'db_file'</span>, <span style="color:orange">'collection_1'</span>, <span style="color:orange">'collection_2'</span>, <span style="color:orange">'cluster_params'</span>]

#### Functions generating the workflows

These firetasks rely on two subworkflows:
- **conv_slabthick_surfene**: Function to start a subworkflow to calculate the optimal thickness for a given bulk structure along a given orientation. Set the computational and cluster parameters and open a subworkflow. <br/>
Input arguments are:
-- `structure`:
-- `mp_id`:
-- `miller`:
-- `functional`:
-- `comp_params`:
-- `db_file`:
-- `collection`:
-- `thick_start`:
-- `thick_incr`:
-- `nsteps`:
-- `vacuum`:
-- `relax_type`:
-- `ext_index`:
-- `slab_name`:
-- `spec`:
-- `cluster_params`:

- **surface_energy**: Function to start a subworkflow to calculate the surface energy of a slab of a list of slabs. Set the computational and cluster parameters and open a subworkflow. <br/>
Input arguments are:
-- `structure`:
-- `mp_id`:
-- `miller`:
-- `functional`:
-- `comp_params`:
-- `db_file`:
-- `collection`:
-- `thick_start`:
-- `thick_incr`:
-- `nsteps`:
-- `vacuum`:
-- `relax_type`:
-- `ext_index`:
-- `slab_name`:
-- `spec`:
-- `cluster_params`:

### Functioning
