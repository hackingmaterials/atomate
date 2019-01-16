=================
atomate Changelog
=================

**v0.8.5**

* add NMR workflow (S. Dwaraknath)
* surface workflow refactor (J. Montoya)
* update elasticity package for latest pymatgen (J. Montoya)

**v0.8.4**

* add some QChem functionality! (S. Blau, E. Spotte-Smith)
* better credential parsing in CalcDB (A. Rutt)
* update paramiko for better security (A. Jain)
* add tasks_settings.yaml to site-packages (A. Jain)

**v0.8.3**

* add CHGCAR and AECAR storage options (J. Shen)

**v0.8.2**

* various requirement updates, code refactorings, and bug fixes especially in Qchem and LAMMPs packages (S. Blau)

**v0.8.1**

* add Bader charge parsing (S. Dwaraknath)
* parse_outputs can push fields to new FWS (S. Dwaraknath)
* minor bug fix (S. Torrisi)
* some testing updates (A. Faghaninia, S. Blau)

**v0.7.9**

* Overhaul and update of QChem in atomate (B. Wood, S. Blau)
* Fix integration tests (A. Dunn)
* Stop officially supporting Py2

**v0.7.8**

* some QChem tasks and workflows (B. Wood)
* add full HSE band structure workflow (A. Jain)

**v0.7.7**

* allow list of ErrorHandler in RunVaspCustodian (A. Jain)
* fix adsorption unit test (J. Montoya)

**v0.7.6**

* remove force checking in drone (S. Dwaraknath)
* cleanups to DOS + BS parsing and insertion in GridFS (S. Dwaraknath)
* code cleanup and bugfixes (S. Dwaraknath, specter119)
* some more docs on offline mode (A. Jain)

**v0.7.5**

* standardize drones to datetime.utcnow() (J. Montoya)
* fixed additional field serialization issue (J. Montoya)
* fix defuse_unsuccessful logic, hat tip to @specter119 for pointing it out (A. Jain)
* some doc updates (A. Jain)

**v0.7.4**

* fix pymatgen dep (A. Jain)

**v0.7.3**

* minor drones updates (S. Dwaraknath)
* installation fix (P. Huck)

**v0.7.2**

* Bugfix database getter in builders (S. Dwaraknath)

**v0.7.1**

* update ``atwf`` to find the path to workflows better (M. Dias Costa)
* better surface workflow naming (A. Jain)

**v0.7.0**

* change default behavior when a run looks OK but is unconverged (A. Jain)
* Some test and code cleanups (S. Dwaraknath, J. Montoya)
* update to FW names when no structure provided (S. Dwaraknath)
* remove boltons dependency (A. Faghaninia)
* fix max_force check for selective dynamics (J. Montoya)

**v0.6.9**

* update requirements to include boltons
* bugfix for atwf (S. Dwaraknath)

**v0.6.8**

* New SCAN functional workflow (S. Dwaraknath)
* remove dependence on pymatgen-db (S. Dwaraknath)
* more bandgap properties parsed by drone (transition and is-direct) (S. Dwaraknath)
* option to clean up large output files like WAVECAR (S. Dwaraknath)
* option to recursively copy file tree in CopyFilesFromCalcLoc (A. Faghaninia)
* bugfix: apply vasp_input_set_params when StaticFW have parents (specter119)
* misc bugfixes (S. Dwaraknath, A. Jain)

**v0.6.7**

* New ferroelectrics workflow! (T. Smidt)
* option to parse LOCPOT in VaspDrone (S. Dwaraknath)
* rename set_fworker -> set_execution_option
* more options for BoltztrapFW (A. Faghaninia)
* misc. bugfixes (D. Broberg, K. Mathew, P. Huck)

**v0.6.6**

* powerup to preserve the same FWorker for all jobs in workflow (S. Dwaraknath)
* DriftErrorHandler in VASP custodian jobs (S. Dwaraknath)
* some FireTasks in anticipation of ferroelectrics workflow (T. Schmidt, A. Jain)

**v0.6.5**

* fix delta_volume_percent, set as new key and update FixTasksBuilder (B. Bocklund, A. Jain)
* drone schema version reflects atomate version (M. Horton)
* unit test fix (J. Montoya)

**v0.6.4**

* add config option for half_kpts_first and max force (A. Jain, S. Dwaraknath)
* better logic for band structure parsing (S. Dwaraknath)
* misc bugfix (P. Huck)

**v0.6.3**

* fix Gibbs wf db insertion (A. Dunn, K. Mathew)
* minor doc updates & fixes (A. Jain)

**v0.6.2**

* Fix LepsFW after prev refactor (A. Jain)
* Doc improvements (A. Jain, B. Bocklund)

**v0.6.1**

* many improvements to documentation (A. Jain, B. Bocklund)
* add DFPTFW (K. Mathew)
* simplify LepsFW - move Raman into RamanFW (K. Mathew)
* copy piezo tensor to output (S. Dwaraknath)

**v0.6.0**

* Gibbs preset workflow and anharmonic contributions (B. Bocklund)
* improvements to packmol workflow (K. Mathew)
* modify_potcar powerup (J. Montoya)
* more metadata in some analysis collections (B. Bocklund)
* ability to specify common params in atwf (A. Jain)
* allow powerups in atwf (J. Montoya)
* many improvements to builders performance (A. Jain)
* updates and fixes to installation tutorial (A. Jain, B. Bocklund)
* unit testing updates (J. Montoya)
* misc fixes ...

**v0.5.8**

* major improvements to LAMMPS workflow (B. Wood, K. Mathew)
* doc updates (B. Bocklund)
* minor cleanups (K. Mathew)

**v0.5.7**

* VASP drone stores original inputs (S. Dwaraknath)
* updates to EELS workflow (K. Mathew)
* misc cleanups (A. Jain, S.P. Ong, K. Mathew)


**v0.5.6**

* major improvements to elastic tensor calculations and compatibility with latest pymatgen (J. Montoya, K. Mathew)

**v0.5.5**

* remove PyPI download size by an order of magnitude

**v0.5.4**

* re-attempt to fix packaging of YAML workflow library in pip

**v0.5.3**

* attempt to fix packaging of YAML workflow library in pip
* update doc links

**v0.5.2**

* band gap estimation builder based on dielectric constants
* clean up pypi packaging (S.P. Ong)
* link to new doc links
* misc bugfixes and workflow settings update/fixes (K. Mathew, A. Jain)

**v0.5.1**

* use ruamel instead of pyyaml (S.P. Ong)
* add magnetic moment parsing of output (M.K. Horton)
* misc cleanups, bug fixes, doc improvements (K. Matthew, S. Dwaraknath, A. Jain)

**v0.5.0**

.. caution:: pymatgen has updated its default kpoint scheme! Kpoint settings will change.

* migration to new pymatgen and new kpoint settings
* much improved docs (B. Bocklund, A. Jain)
* *major* code cleanup (J. Montoya, K. Mathew, A. Jain)
* many unit test updates (A. Faghaninia, H. Tang, S.P. Ong, A. Jain)
* fix automated testing on pull requests (K. Mathew)
* misc fixes


**v0.4.5**

* *extensive* code review, code cleanup, and improved code docs - with some minor name refactoring
* new builders: dielectric, structureanalysis (currently gives dimensionality of structure)
* rewrite powerups as in-place with cleaner syntax
* improved installation tutorial (B. Bocklund)
* improve/fix/reorganize some unit tests
* bug fixes (A. Jain, H. Tang, K. Mathew, B. Bocklund)

**v0.4.4**

* NEB workflow (H. Tang)
* adsorption workflow (J. Montoya)
* improvements to Gibbs workflow (K. Mathew)
* misc bugfixes, improvements (A. Faghaninia, A. Jain)

**v0.4.3**

* Add Gibbs energy w/volume (K. Mathew)
* Draft EXAFS workflow (K. Matthew)
* Add slater-gamma formulation to compute the Gruneisen parameter (K. Matthew)
* gamma vasp powerup (S. Dwaraknath)
* More options for elasticity WF (J. Dagdalen)
* Add StdErrorHandler to handlers (A. Jain)
* Auto-detect and remove line_mode parameter in MMVaspDB (A. Jain)
* added unit tests
* misc cleanup, refactoring, and doc udpates
* misc bugfixes


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