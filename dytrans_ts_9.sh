#!/bin/bash

SHOTNO=$1

eval "matlab -nodesktop -nojvm -r 'dytrans_ts_9($1); exit'"

exit