# What is this directory?

This directory contains materials science workflows that chains together
individual Fireworks defined within atomate. The workflows are defined
in a special YAML format (different than the typical FireWorks format) that
is intended to be simple and easy to modify.

# What is the format of the workflows?

The overall format is YAML.

At the top there is usually a comment describing the workflow. The rest of
the workflow looks like this:

fireworks:
- fw: atomate.vasp.fireworks.core.OptimizeFW
- fw: atomate.vasp.fireworks.core.StaticFW
  params:
    parents: 0
- fw: atomate.vasp.fireworks.core.NonSCFUniformFW
  params:
    parents: 1
- fw: atomate.vasp.fireworks.core.NonSCFLineFW
  params:
    parents: 1
common_params:
  db_file: db.json
  $vasp_cmd: $HOME/opt/vasp
name: bandstructure

The `fireworks` key is a list of Fireworks; it is expected that all such
Fireworks (when you dig into the Python code) have "structure" as the
first argument and other optional arguments following that. Each Firework
is specified via "fw": <explicit path>.

You can pass arguments into the constructor using the special
keyword `params`, which is a dict. Any param starting with a $ will
be expanded using environment variables.If multiple fireworks share
the same `params`, you can use `common_params` to specify a common
set of arguments that are passed to all fireworks. Local params
take precedent over global params.

Another special keyword is `parents`, which provides
the *indices* of the parents of that particular Firework in the
list. e.g. a parent of 0 means the first Firework in the list is the parent.
Multiple parents are allowed. This allows you to link the Fireworks into a
logical workflow.

Finally, `name` is used to set the Workflow name
(structure formula + name) which can be helpful in record keeping.