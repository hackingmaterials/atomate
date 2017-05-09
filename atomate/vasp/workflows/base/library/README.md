# What is this directory?

This directory contains materials science workflows that chains together
individual Fireworks defined within atomate. The workflows are defined
in a special YAML format (different than the typical FireWorks format) that
is intended to be simple and easy to modify.

# What is the format of the workflows?

The overall format is YAML.

At the top there is usually a comment (hashtag) describing the workflow.
The rest of the workflow looks like this:

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
metadata:
  tag: testing_workflow

The `fireworks` key is a list of Fireworks; it is expected that all such
Fireworks (when you dig into the Python code) have "structure" as the
first argument and other optional arguments following that. Each Firework
is specified as "fw: <path.to.firework.in.atomate>"".

You can pass arguments into the Firework constructor using the special
keyword `params`, which is a dict. Any param starting with a $ will
be expanded using environment variables. If multiple fireworks share
the same `params`, you can use `common_params` to specify a common
set of arguments that are passed to all fireworks. Local params
take precedent over global params.

Another special keyword is `parents`, which provides
the *indices* of the parents of that particular Firework in the
list. The indices start at zero, i.e, the first Firework in your list
has zero. Thus, if you want the second Firework in the list to be a child
of the first Firework, you should specify a parent of 0 for the Firework.
Multiple parents are allowed. This allows you to link the Fireworks into a
logical workflow.

In the above example, we have:
* the first Firework (OptimizeFW) will run before anything else
* the second Firework (StaticFW) will run after the OptimizeFW is complete
* the third and fourth Fireworks (NonSCFUniformFW and NonSCFLineFW) will
run after the StaticFW is complete. Note these two Fireworks can run in parallel.

Next, `name` is used to set the Workflow name (structure formula +
name) which can be helpful in record keeping.

Finally, one can specify a `metadata` key as a YAML dict/hash that will
initialize workflow metadata - this is purely optional and for bookkeeping.