#!/bin/bash

if test $# -lt 2 ; then
  echo "Usage: $0 <procs> <command> [args]"
  exit 1
fi

procs=$1
command="$2"
shift 2
mpiexec -np $procs konsole --workdir $PWD -e gdb --eval-command=run --args "$command" "$@"
