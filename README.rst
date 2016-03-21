==========
MatMethods
==========

MatMethods is a library for Materials Science workflows. It is currently in development and **unsupported**.

To cite MatMethods, you can cite the following two papers::

    (1) (1) Jain, A.; Ong, S. P.; Chen, W.; Medasani, B.; Qu, X.; Kocher, M.; Brafman, M.; Petretto, G.; Rignanese, G.-M.; Hautier, G.; Gunter, D.; Persson, K. A. FireWorks: a dynamic workflow system designed for high-throughput applications, Concurr. Comput. Pract. Exp., 2015, 22, doi:10.1002/cpe.3505.

    (2) Ong, S. P.; Richards, W. D.; Jain, A.; Hautier, G.; Kocher, M.; Cholia, S.; Gunter, D.; Chevrier, V. L.; Persson, K. a.; Ceder, G. Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis, Comput. Mater. Sci., 2013, 68, 314–319, doi:10.1016/j.commatsci.2012.10.028.


========
Overview
========

MatMethods is still not ready for production. However, the following notes might be helpful.

============
Installation
============

Although MatMethods itself should be pip-installable, the actual process to get an infrastructure to run VASP at supercomputers might be a bit more involved (not difficult, just more steps). See "Detailed Installation Notes" at the end of this document for some notes on how to do this.

==============================
Testing the VASP functionality
==============================

In order to use the VASP functionality, make sure you set up VASP_PSP_DIR variable (see pymatgen docs). Also, make sure you have MongoDB running in order to execute all the tests.

To test the VASP functionality, run the unit tests in ``matmethods.vasp.tests``. These unit tests are designed to run without installing VASP. Some of them start with a VASP workflow but apply the ``make_fake_workflow`` method to replace calling the VASP executable with a "Faker" that verifies basic properties of the inputs and copies pre-stored output files to the current directory, thus simulating the execution of VASP.

The unit tests in matmethods/vasp/tests/test_vasp_workflows.py can be modified to actually run VASP by setting VASP_CMD to a String representing your VASP command.

Many tests have a DEBUG option that can sometimes help in finding problems. Sometimes you need to toggle DEBUG on/off a couple of times if you are doing this to make sure all the old data is actually cleared between debug runs; the tearDown() and setUp() methods are stil a bit finicky.

==========================
Learning to use MatMethods
==========================

If you are familiar with (i) VASP, (ii) pymatgen, (iii) custodian, and (iv) FireWorks, then most of MatMethods such be fairly straightforward. For example, the workflows in matmethods/vasp/examples/vasp_workflows.py should make sense to you. A good place to begin is by reviewing these workflows and confirming that this is the case.

There are only a couple of new concepts in MatMethods that you might need to familiarize yourself with, and they are described below.

The "env_chk", e.g. ``>>db_file<<``
===================================

One issue in coding workflows is what to do when different machines require different settings. For example, the path to the VASP executable or the path to a file containing database credentials might be located in different places on different machines. For users wanting to run on multiple machines, such parameters cannot be hard-coded. However, users that are running on a single machine, or those that are testing things out, might prefer to hard-code those parameters.

The ``env_chk`` functionality is a way to support both hard-coding of parameters as well as letting the machine (or more specifically, the FireWorker) set the parameter. Many of the FireTasks in MatMethods, e.g., ``RunVaspDirect``, state in the docs that they "support ``env_chk``" for a parameter such as ``vasp_cmd``. What this means is that you have two options for creating the FireTask:

Option 1 is to use something like ``my_task = RunVaspDirect(vasp_cmd="vasp")``. This behaves exactly as you would expect in regular Python, i.e., the string literal "vasp" set as the ``vasp_cmd`` parameter.

Option 2 is to use the ``env_chk`` notation which looks like this: ``my_task = RunVaspDirect(vasp_cmd=">>my_vasp_cmd<<")``. If ``env_chk`` parameters like `vasp_cmd`` are enclosed in the ``>><<`` symbols, it is interpreted that the user wants to get the values from the FireWorker's ``env`` value. That is, when executing the workflow, one must use a FireWorker that contains an env that looks like ``{"my_vasp_cmd": "mpirun -n 24 vasp"}``. Here, the ``my_vasp_cmd`` in the dictionary matches the ``>>my_vasp_cmd<<`` string in the env_chk. Thus, when VASP is executed via this FireWorker, it will execute the command ``mpirun -n 24 vasp``. Other FireWorkers might execute different VASP commands and can support this by setting a different value of the FireWorker ``env``. The workflow can be kept intact since the workflow is merely pointing to the ``my_vasp_cmd`` env variable and not setting the VASP command explicitly. There are more details about setting the FireWorker env variables in the FireWorks tutorials (in particular the Worker tutorial). The unit tests also use the env_chk feature to find the db configuration file. e.g., see the unit test: ``matmethods.vasp.tests.test_vasp_workflows.TestVaspWorkflows#test_single_Vasp_dbinsertion`` and you will have a flavor for how this works. Just remember that if you see something like this ``>>db_file<<``, when running your Workflow your FireWorker will need to set the env like this: ``FWorker(env={"db_file": "path/to/db.json"})`` and you will need to use that FireWorker when launching the jobs.

VaspLocs
========

If you are running multiple VASP jobs that depend on copying the outputs of previous jobs, one issue is how to pass the directory information of previous VASP jobs from Firework to Firework. It is possible to do this manually (as was done in the MPWorks codebase), or using the ``pass_job_info`` keyword built into Fireworks, but the standard way to do this in MatMethods is *VaspLocs*. Procedurally, all you need to do is add the ```PassVaspLocs`` FireTask to every Firework that contains a VASP job (see ``matmethods.vasp.examples.vasp_workflows.get_wf_bandstructure_Vasp`` for an example). Downstream jobs like ``CopyVaspOutput`` will have a ``vasp_loc`` variable that can be set to True, and will automatically get the previous VASP dir parsed from before. Similar with ``VaspToDbTask``. Note that a couple of advantages of this system are:

* It is a general way of passing VASP dirs that works with any Firework, and doesn't require you to code the logic of passing VASP directories inside of other functions (e.g., database insertion tasks as was done previously in MPWorks). Thus, the task of reporting and passing the VASP job location is well-separated from the other functions and can just be added in very easily. The only downside is that you have to remember to add in this FireTask.
* The VASPLocs maintains a running dictionary of job type to job location. If you need to grab outputs from multiple jobs (or say, from two jobs back), it is all supported within the framework. Just read the docs, e.g., of ``CopyVaspOutput``.
* In the future, if the job directories are located across different machines and require ``scp`` or some other complex transfer mechanism, that can be built into this infrastracture and the user does not need to worry about it.

Workflow "Powerups"
===================

Workflow powerups are intended to be like function decorators, but for Workflows. For example, let's say you've built a multi-step workflow that computes a band structure. Now, you to make sure that once a workflow starts running, it is prioritized to finish that workflow versus starting other workflows. By passing your workflow through a "powerup", you can get back a decorated workflow that sets the priorities of the Fireworks inside your workflow to endow this behavior (e.g., give all children Fireworks 2X the priority of the root parent). Another planned "powerup" is to endow Workflows with duplicate checking, i.e., to make sure the same structure is not run twice. In the past, such duplicate checking logic would be mixed in with the rest of the Workflow (about setting up VASP parameters, running VASP, etc.), and the end result was a very messy workflow code. It was also difficult to turn duplicate checking off and on as desired since all the logic was intermixed. By moving the duplicate checking to a "powerup", one can simply enable duplicate checking by passing the Workflow through the appropriate powerup.

See the "vasp_powerups.py" file for examples.


===========================
Detailed Installation Notes
===========================

Here are some notes on how to get MatMethods up and running in a production system at your supercomputing center. These notes are geared towards the NERSC supercomputing center. You'll need to fill in details and adapt accordingly for other centers. Hopefully, you are not a complete beginner.

Things you need to do once
==========================

There are some things.

Blah
----


