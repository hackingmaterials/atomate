====================
atomate Changelog
====================

**v0.4.2**

.. caution:: The ``tags_fws`` powerup now has different default parameters!

* updates to piezo workflow (S. Dwaraknath)
* formation energy to Ehull builder (A. Faghaninia)
* tag_fws is more general (A. Faghaninia)
* updates for PMG naming schemes for vars (A. Jain)
* boltztrap runs can add tags (A. Faghaninia)
* can filter which tasks are used in materials builder (A. Faghaninia, A. Jain)

**v0.4.1**
* more fixes for elastic workflow (J. Montoya)
* more validation for VASP runs (A. Faghaninia)
* more flexible ObjectId insertion (A. Faghaninia)
* misc doc updates (A. Jain)

**v0.4**
* rename of "MatMethods" to atomate(!) (A. Jain)
* bulk modulus workflow and equation of state (K. Matthew)
* add features to Gibbs workflows (K. Matthew)
* elastic workflow updates (J. Montoya, K. Matthew)
* Spin orbit coupling (A. Faghaninia)
* HSE line-mode band structure workflow (A. Faghaninia)
* Feff workflows (K. Matthew)
* bug fixes (K. Matthew)
* much code refactoring, cleanup, and many minor improvements (K. Matthew, A. Jain, J. Montoya, S.P. Ong, B. Bocklund, A. Faghaninia)

**v0.3**

* Raman workflow (K. Mathew)
* Gibbs workflow (K. Mathew)
* More efficient task builder (S. Ong)
* tag workflows and add_trackers powerups (A. Jain, A. Faghaninia)
* refactor elastic workflow (K. Mathew)
* bugfixes and tools package (K. Mathew)

**v0.21**

* Lammps workflows and packmol support (K. Mathew)
* Rework some of the RunVaspFake code (K. Mathew)
* Fixes to elastic workflow (J. Montoya)
* Minor refactoring (K. Mathew)
* Minor MD workflow updates (M. Aykol)
* Fix builder for HSE gap and add chemsys (A. Jain)
* WF metadata powerup (A. Jain)
* Minor bug fixes and misc. improvements (K. Mathew, J. Montoya, A. Faghaninia)

**v0.2**

* BoltzTraP transport workflows (A. Jain)
* major builder improvements (merge multiple collections, progressbar, config, more...)
* use FrozenJobErrorHandler by default (A. Jain)
* add basic configuration overrides for preset workflows (A. Jain)
* misc improvements and bugfixes (A. Jain, K. Mathew)
* py3 compatibility fixes (K. Mathew)

**v0.1**

* add some builders
* elastic + piezo workflows (J. Montoya + S. Dwaraknath)
* minor doc improvements (A. Faghaninia)
* misc code improvements and bug fixes, plus upgrades for new pymatgen (A. Jain)

**v0.0.3**

* initial release (A. Jain, S.P. Ong, K. Mathew, M. Aykol)