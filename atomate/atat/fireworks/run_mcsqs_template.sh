#!/bin/bash

# this shell script is a template: it's not meant to be run directly!
# please see McsqsFW

# write version string to file for later database ingestion
mcsqs -v 2>&1 | head -n 1 > mcsqs_version.txt

# tracker command: every minute aggregate results and append to tracker file
tracker_cmd="
while sleep 60
do
    mcsqs -best | tail -n 1 | awk '{print strftime(\"%s\") \" \" \$$2}' >> mcsqs_tracker.txt
done
"

# run tracker in background, set same time limit as sqs processes
${timeout_cmd} -s 0 ${time}s bash -c "$${tracker_cmd}" &

# run mcsqs in parallel
for (( id=0 ; id<${ncores} ; id ++ ))
do
    ${timeout_cmd} -s 0 ${time}s mcsqs -n ${size} ${settings} -ip=$$id &
done

# keep Firework running until mcsqs times out
sleep ${time}

# final aggregation of best SQS results
mcsqs -best