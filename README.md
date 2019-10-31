[![Build Status](https://travis-ci.org/EIT-team/PEITS.svg?branch=master)](https://travis-ci.org/EIT-team/PEITS)

# Install guide
The install guide can be found in the Wiki (correct as of 16.1.2017) [here](https://github.com/EIT-team/PEITS/wiki/Installing-PEITS)

# Research Paper
A detailed description of the forward model can be found in https://dx.doi.org/10.1109/TBME.2014.2342280

# Using the Docker container to deploy and run on a HPC facility

Running Docker containers is generally not supported on HPC systems, mostly for security reasons. 
Instead, another containerisation technology called Singularity is preferred.

## Building a container from an image on DockerHub

Singularity allows you to build a container from remote repositories, including
Docker images on DockerHub. The following process has been successfully tested on UCL's
central HPC system Myriad but it might vary when running on another facility (mostly due to Singularity versions, etc.)

With Singularity v2, convert a Docker image to a Singularity image via

```
singularity build <IMAGE_NAME>.sqsh docker://<USERNAME>/<IMAGE>:<TAG>
```

In our case, the Docker image can be found on [DockerHub](https://hub.docker.com/r/uclrits/peits), the tag being 16.04. 

With Singularity v3, you'll probably want to use the new SIF container format:

```
singularity build <IMAGE_NAME>.sif docker://<USERNAME>/<IMAGE>:<TAG>
```

On shared computing environments, the image file is best placed on a shared
network/cluster file system (not in your home directory).

## Running a container

You can log in to an interactive shell within the container with the command

```
singularity shell <IMAGE_NAME>
```

_Note_: remember that Singularity containers are by default read-only, unless
you have created a sandbox one (Singularity v3 only).  In this case, you can use
the option `--writable` to run the container in read-write mode (this will
effectively modify the container), or `--writable-tmps` if you want to mount the
file system as read-write but with non-persistent data.

If you are in an HPC system and you want to submit a job to the queue you need
to run a non-interactive command.  You can execute a specific `<COMMAND>` within
the container with:

```
singularity exec <IMAGE_NAME> <COMMAND>
```

If your program is going to write some data and you cannot use the `--writable`
or `--writable-tmpfs` options you have to bind to the container some local
directories that you can read/write.  You can do this with the `-B`/`--bind`
option:

```
singularity exec \
    -B /path/to/local/dir1:/mount/point1:ro \
    -B /path/to/local/dir2:/mount/point2:rw \
    <IMAGE_NAME> <COMMAND>
```

In this example, `/path/to/local/dir1` will be mounted on `/mount/point1` in
read-only mode, while `/path/to/local/dir2` will be mounted on `/mount/point2`
in read-write mode.

### Real-world example: container with MPI support on HPC cluster

This section describes how we ran a Singularity container with MPI on Myriad for
the Stroke project.  This may not reflect what you will need to do but shows
some details that you may need to take into account when running a Singularity
container on an HPC system.

This is a simplified version of the job script that we used for the Stroke
project:

```sh
#!/bin/bash -l
#$ -l h_rt=0:30:00
#$ -pe mpi 16
#$ -cwd
THIS_DIR="${HOME}/Scratch/stroke" # This must be hard-coded
# Use a unique ID in order to launch multiple jobs using the
# container at the same time
ID="${JOB_ID:-${RANDOM}}"
OUTPUT_DIR="${THIS_DIR}/PEITS-tmp/${ID}/output"
echo "Creating \"${OUTPUT_DIR}\"..."
mkdir -p "${OUTPUT_DIR}"
# Flags that will be used only when some variables are defined
FLAGS=()
if [[ -n "${TMPDIR}" ]]; then
    # MPI wants to write to `${TMPDIR}` if available.
    # If so, bind it in read-write mode to `${SINGULARITY_TMPDIR}`
    FLAGS+=(-B "${SINGULARITY_TMPDIR}":"${TMPDIR}":rw)
fi
if [[ -n "${PE_HOSTFILE}" ]]; then
    PE_DIR="$(dirname "${PE_HOSTFILE}")"
This conversation was marked as resolved by dpshelio
    # MPI needs to read `${PE_HOSTFILE}` if available, mount its directory
    FLAGS+=(-B "${PE_DIR}":"${PE_DIR}":ro)
fi
# Run singularity
singularity exec --no-home -H /home/peitsier \
	    -B "${HOME}/software/PEITS/data":/build/PEITS/data:ro \
	    -B "${HOME}/software/PEITS/partitions":/build/PEITS/partitions:ro \
	    -B "${THIS_DIR}/PEITS-tmp/run":/build/PEITS/run:ro \
	    -B "${OUTPUT_DIR}":/build/PEITS/output:rw \
	    "${FLAGS[@]}" \
	    "${THIS_DIR}/peits.sqsh" \
	    /build/PEITS/run/run_test.sh
```

This is a bit convoluted but the code is commented.  We had to mount a few
directories in read-only mode that contained some data and were not within the
container.  We then mounted in read-write mode the `output` directory, where the
output files will be placed.  The directory mounted on `/build/PEITS/run`
contained the shell script that we wanted to run.  This contained the following
call to `mpirun`:

```sh
mpirun -np ${NSLOTS:-1} "/path/to/the/program"
```

The environment variable `NSLOTS` is set by Myriad's scheduler to the number of
cores requested for the job.  Thus, `mpirun` would spawn `NSLOTS` processes,
defaulting to 1 if the variable is not set.

The `FLAG` variable is set dynamically depending on some variables available in
the environment.  This allowed us to run the same script on the login node for
some quick and not resource-demanding test.
