==========
MatMethods
==========

MatMethods is a library for Materials Science workflows. It is currently in development and **unsupported**.

To cite MatMethods, you can cite the following two papers::

    (1) (1) Jain, A.; Ong, S. P.; Chen, W.; Medasani, B.; Qu, X.; Kocher, M.; Brafman, M.; Petretto, G.; Rignanese, G.-M.; Hautier, G.; Gunter, D.; Persson, K. A. FireWorks: a dynamic workflow system designed for high-throughput applications, Concurr. Comput. Pract. Exp., 2015, 22, doi:10.1002/cpe.3505.

    (2) Ong, S. P.; Richards, W. D.; Jain, A.; Hautier, G.; Kocher, M.; Cholia, S.; Gunter, D.; Chevrier, V. L.; Persson, K. a.; Ceder, G. Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis, Comput. Mater. Sci., 2013, 68, 314–319, doi:10.1016/j.commatsci.2012.10.028.


Testing the VASP functionality
==============================

Some of the things you need to do are:

1. Set up VASP_PSP_DIR variable

To test the VASP functionality, run the unit tests in ``matmethods.vasp.firetasks.tests``. Note that some of the tests will require the VASP executable be installed (just read all the messages, they are self-explanatory).

