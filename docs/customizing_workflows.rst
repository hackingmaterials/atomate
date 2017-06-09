.. title:: Customizing workflows
.. _customizing workflows:

=====================
Customizing Workflows
=====================

Introduction
============

The :ref:`creating workflows` guide gave details on constructing workflows. This group of tutorials will give specific examples for customizations to workflows as you create them.

For some of these customizations, preset workflows *cannot* be used. Preset workflows are designed to give generically reasonable options. More full access for customizing the workflows can be achieved by using the :py:mod:`atomate.vasp.workflows.base` workflows instead of the presets.

Objectives
==========

* Provide examples for customizating workflows

Prerequisites
=============

It's best if you are able to create workflows on your own in Python. See the :ref:`creating workflows guide <creating workflows>`

.. _powerups:

Powerups
========

Powerups are all designed to be used as functions where you pass in your original workflow and other keyword arguments and get back the modified workflow. An example is shown below, but does not show all of the powerups. To see more powerups go to the powerups documentation for the package you are using, e.g. VASP is :py:mod`atomate.vasp.powerups`.

An example for adding an INCAR setting to use a different force convergence criteria for the only the structure optimization in the elastic workflow is

.. code-block:: python


    from atomate.vasp.workflows.presets.core import wf_elastic_constant
    from atomate.vasp.powerups import add_modify_incar
    from pymatgen import Structure

    # load your structure, e.g. from a POSCAR
    struct = Structure.from_file('POSCAR')

    # create the workflow
    orig_wf = wf_elastic_constant(struct)


    # use the powerup to change any Fireworks with 'optimization' in the name to use EDIFFG=-0.05
    # note: the 'incar_update' is *required* if you want to update
    # note: not passing the ``modify_incar_params`` keyword will result in the Firework getting
    #       the 'incar_update' key and values from your FireWorker's env
    modified_wf = add_modify_incar(orig_wf, modify_incar_params={'incar_update': {'EDIFFG': -0.05}},
                                   fw_name_constraint='optimization')

    # print if you want to check the tasks. Warning: several lines long.
    # print(orig_wf.fws[0].tasks)
    # print(modified_wf.fws[0].tasks)

VASP Calculation Settings
=========================

Most VASP calculation-specific settings (e.g. those from INCAR, KPOINTS, and POTCAR files) are controlled by `pymatgen's vasp sets`_ . VASP workflows take ``vasp_input_set`` options and you can directly import and use them or customize them before using them in atomate workflows.

.. note:: Using the ``vasp_input_set`` or ``vis`` keywords in workflow constructors usually only controls the first Firework that uses that set. If you want to have multiple Fireworks use custom input sets (or just not the first one, e.g. in a bandstructure workflow) then you have to make a custom workflow yourself.

Using a different functional
----------------------------

To use a different functional, for instance in a optimization calculation, you can do the following:

.. code-block:: python


    from fireworks import Workflow
    from atomate.vasp.fireworks.core import OptimizeFW
    from pymatgen.io.vasp.sets import MPRelaxSet
    from pymatgen import Structure

    def get_optimize_wf(structure, name="optimization wf", vasp_input_set=None,
                        vasp_cmd="vasp", db_file=None, user_kpoints_settings=None,
                        tag="", metadata=None):
        """
        Returns a structure optimization workflow.

        Args:
            structure (Structure): input structure to be optimized and run
            name (str): some appropriate name for the transmuter fireworks.
            vasp_input_set (DictSet): vasp input set.
            vasp_cmd (str): command to run
            db_file (str): path to file containing the database credentials.
            user_kpoints_settings (dict): example: {"grid_density": 7000}
            tag (str): some unique string that will be appended to the names of the fireworks so that
                the data from those tagged fireworks can be queried later during the analysis.
            metadata (dict): meta data

        Returns:
            Workflow
        """
        # input set for relaxation
        vis_relax = vasp_input_set or MPRelaxSet(structure)
        if user_kpoints_settings:
            v = vis_relax.as_dict()
            v.update({"user_kpoints_settings": user_kpoints_settings})
            vis_relax = vis_relax.__class__.from_dict(v)

        # Structure optimization firework
        fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                          db_file=db_file, name="{} structure optimization".format(tag))]

        wfname = "{}:{}".format(structure.composition.reduced_formula, name)

        return Workflow(fws, name=wfname, metadata=metadata)

    # load your structure, e.g. from a POSCAR
    struct = Structure.from_file('POSCAR')

    # create a custom input set
    my_custom_input_set = MPRelaxSet(struct, potcar_functional='LDA')

    # create the workflow
    my_wf = get_optimize_wf(struct, vasp_input_set=my_custom_input_set)

For the supported options, see the VASP documentation and `pymatgen's vasp sets`_ documentation. PBE (default), LDA, PW91, LDA_US were supported at the time of writing.


Custom KPOINTS settings
-----------------------

KPOINTS settings can also be similarly customized using the above example. You can control them with the following keywords (from `pymatgen's vasp sets`_):

* ``force_gamma``: always use gamma centered kpoint generation. Default (False) is to use Automatic Density kpoint scheme, which will use the Gamma centered generation scheme for hexagonal cells, and Monkhorst-Pack otherwise.
* ``user_kpoints_settings``: Override kpoints setting by supplying a dict. E.g., ``{"reciprocal_density": 1000}``. Other options are ``grid_density`` or ``length``.

.. code-block:: python

    from pymatgen.io.vasp.sets import MPRelaxSet
    from pymatgen import Structure

    # load your structure, e.g. from a POSCAR
    struct = Structure.from_file('POSCAR')

    # create a custom input set
    my_custom_input_set = MPRelaxSet(struct, force_gamma=True, {"grid_density": 10} )

    # create the workflow
    my_wf = get_optimize_wf(struct, vasp_input_set=my_custom_input_set)

If you need more control, create the ``Kpoints`` object directly with pymatgen. It is flexible and only a brief example will be shown. See the `full Kpoints documentation`_ for more

.. code-block:: python

    from pymatgen.io.vasp.sets import MPRelaxSet
    from pymatgen.io.vasp.inputs import Kpoints
    from pymatgen import Structure

    # load your structure, e.g. from a POSCAR
    struct = Structure.from_file('POSCAR')

    # the simples way to do this is to create a subclass of the input set you want
    # and override the kpoints property to return what you want.
    class MyInputSet(MPRelaxSet):
        def __init__(self, structure, points=(5,5,5), shift=(0,0,0), **kwargs):
            super(MPRelaxSet, self).__init__(structure, MPRelaxSet.CONFIG, **kwargs)
            self.points = points
            self.shift = shift

        @property
        def kpoints(self):
            # choose either of these
            # use Monkhorst-Pack scheme
            return Kpoints.monkhorst_automatic(kpts=self.points, shift=self.shift)
            # use a Gamma centered scheme
            return Kpoints.gamma_automatic(kpts=self.points, shift=self.shift)

    # create an instance of the custom input set
    my_custom_input_set = MyInputSet(struct, points=(5,5,5), shift=(1,1,1))
    # show that the set applied
    print(my_custom_input_set.kpoints)

    # create the workflow
    my_wf = get_optimize_wf(struct, vasp_input_set=my_custom_input_set)


.. _full Kpoints documentation: http://pymatgen.org/pymatgen.io.vasp.inputs.html#pymatgen.io.vasp.inputs.Kpoints



Custom INCAR settings
---------------------

Custom INCAR settings can also be accomplished using ``VaspInputSet`` objects, but it is often more efficient to use a `add_modify_incar Powerup <powerups>`_


Use a different POTCAR
----------------------

Which POTCAR file you want to use is controlled by the input set as well. The easist way to control it is by updating the ``config_dict`` dictionary of your input set.

.. code-block:: python

    from pymatgen.io.vasp.sets import MPRelaxSet
    from pymatgen import Structure

    # load your structure, e.g. from a POSCAR
    struct = Structure.from_file('POSCAR')

    # create a custom input set
    my_custom_input_set = MPRelaxSet(struct)
    print('Config dict example: {}\n'.format(my_custom_input_set.config_dict))
    print('Before change: {}'.format(my_custom_input_set.config_dict['POTCAR']['Mg']))
    my_custom_input_set.config_dict['POTCAR']['Mg'] = 'Mg'
    print('After change: {}'.format(my_custom_input_set.config_dict['POTCAR']['Mg']))

    # create the workflow
    my_wf = get_optimize_wf(struct, vasp_input_set=my_custom_input_set)

.. warning:: Make sure not to try a nested dictionary update (e.g. ``my_custom_input_set.config_dict.update({'POTCAR': {'Mg': 'Mg'}})`` )! It will wipe out all of the other ``POTCAR`` entries in the dict.


.. _pymatgen's vasp sets: http://pymatgen.org/pymatgen.io.vasp.sets.html
.. _pymatgen.io.vasp.sets.MPHSERelaxSet: http://pymatgen.org/pymatgen.io.vasp.sets.html#pymatgen.io.vasp.sets.MPHSERelaxSet

