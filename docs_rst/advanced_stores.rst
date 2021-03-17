
.. title:: Advanced Storage Stratagies
.. _advanced_storage:

==================
Advanced Storage Stratagies
==================

Storing data greater than 16 Mb 
============

Introduction
------------

As the analysis capabilities of the high-throughput computational material science grows, more and more data is required to perform the analysis.
While the atomic structure and metadata about the calculations can be stored within the 16 Mb limit of MongoDB's data framework the storage requirement for field data like the charge density is often orders of magnitude larger than that limit.
As such, alternative methods are required to handle the storage of specific fields in the parsed task document.
For these fields, with larger storage requirements, the data itself will be removed from the task document and replaced with and ``fs_id`` to reference where the data is stored.
This results in a task document with a much smaller size, which will still be uploaded onto the MongoDB specified in your ``DB_FILE``

Two approaches are currently available in atomate to handle the storage of large data chunk.
The user can implement a GridFS-based chunking and storage procedure as we have done in ``VaspCalcDb``.
But the recommended method for large object storage is to use ``maggma`` stores implemented in the ``CalcDb`` class.
Currently, only the Amazon S3 store is implemented.

Please read the documentation for ``maggma`` for more details _maggma: https://materialsproject.github.io/maggma

Configuration
------------

To storing the larger items to an AWS S3 bucket via the ``maggma`` API, the user needs to have the ``maggma_store`` keyword in present in the ``db.json`` file used by atomate.

.. code-block:: json

    {
        "host": "<<HOSTNAME>>",
        "port": <<PORT>>,
        "database": "<<DB_NAME>>",
        "collection": "tasks",
        "admin_user": "<<ADMIN_USERNAME>>",
        "admin_password": "<<ADMIN_PASSWORD>>",
        "readonly_user": "<<READ_ONLY_PASSWORD>>",
        "readonly_password": "<<READ_ONLY_PASSWORD>>",
        "aliases": {}
        "maggma_store": {
                "bucket" : "<<BUCKET_NAME>>",
                "s3_profile" : "<<S3_PROFILE_NAME>>",
                "compress" : true,
                "endpoint_url" : "<<S3_URL>>"
        }
    }

Where ``<<BUCKET_NAME>>`` is S3 bucket where the data will be stored, ``<<S3_PROFILE_NAME>>`` is the name of the S3 profile from the ``$HOME/.aws`` folder.
Note, this AWS profile needs to be available anywhere the ``VaspCalcDb.insert_task`` is called (i.e. on the computing resource where the database upload of the tasks takes place).

Usage
-----------

Example: store the charge density 

To parse a completed calculation directory.  We need to instantiate the ``drone`` with the ``parse_aeccar`` or ``parse_chgcar`` flag.

.. code-block:: python

    calc_dir = "<<DIR>>/launcher_2019-09-03-11-46-20-683785"
    drone = VaspDrone(parse_chgcar=True, parse_aeccar=True)
    doc = drone.assimilate(calc_dir)
    task_id = mmdb.insert_task(doc)

Some workflows like the ``StaticWF`` will pass the parsing flags like ``parse_chgcar`` to the drone directly.

To access the data using the task_id we can call

.. code-block:: python
    
    chgcar = mmdb.get_chgcar(task_id)

Similar functionalities exist for the band structure and DOS.
Please refer to the documentation of ``VaspCalcDb`` for more details.
