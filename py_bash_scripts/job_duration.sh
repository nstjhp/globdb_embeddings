#!/bin/bash

# Get sacct output with Start and End times
# If the cluster has had a problem for example, then elapsed might not = duration.
# Elapsed does not include time spent waiting in the queue, time spent on job setup/teardown, or any delays caused by system overhead.
sacct -S $(date -d "last Friday 12:00" +%Y-%m-%dT%H:%M) -o JobID,Partition,NodeList,Elapsed,State,Start,End -n -P | grep COMPLETED | while IFS='|' read -r jobid partition nodelist elapsed state start end; do
    # Skip the header line
    if [[ "$jobid" == "JobID" ]]; then
        echo "$jobid|$partition|$nodelist|$elapsed|$state|$start|$end|StartToEnd"
        continue
    fi

    # Calculate Start to End Time
    start_epoch=$(date -d "$start" +%s)
    end_epoch=$(date -d "$end" +%s)
    duration_seconds=$((end_epoch - start_epoch))
    duration=$(date -u -d @$duration_seconds +%T)

    # Print the result
    echo "$jobid|$partition|$nodelist|$elapsed|$duration|$start|$end"
done | column -t -s '|' | less
