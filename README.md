[![Build Status](https://travis-ci.org/EIT-team/PEITS.svg?branch=master)](https://travis-ci.org/EIT-team/PEITS)

# Install guide

Follow [the install guide](./INSTALL.md) or use the provided docker container as
detailed below.

You can build the container using:
```bash
docker build -t uclrits/peits:verion-1.1 .
```

or fetch the one for this release as:
```bash
docker pull uclrits/peits:version-1.1
```

Run then the container as:
```bash
docker run -it --name mypeits -v /full/localpath:/mydata  uclrits/peits:version-1.1
```
Where `/full/localpath/` is the path to a directory on your machine (the host)
that you want to see it from within the container (under `/mydata`).
This command will login you into the container. If you exit you will stop it.
If you want to restart it and get inside again:

```bash
docker start mypeits
docker attach mypeits
```

To delete the container when you are done with use
```bash
docker rm mypeits
```
which will delete any content you have written out of the mounted location (_i.e._, `/mydata`).

Alternatively, you can use [docker compose](https://docs.docker.com/compose/) to
build and run the container:

```bash
docker-compose up # to build it
docker-compose run peits # to get inside
docker-compose down # to stop it
```

Change the field `volumes` on `docker-compose.yml` to specify which directory
you want to mount on the docker container.

## Running PEITS on the container

Browse to the PEITS directory and run `dune_peits` as:

```bash
cd /build/PEITS/src
mpirun -np 2 ./dune_peits
```

You need first to set the configuration accordingly to your needs. There are
two example files: `data/parameter_example.cfg` and `data/standardparams_example.cfg`.
Edit them as desired and rename them to `parameter.cfg` and `standardparams.cfg`
respectively.

For more details check the scripts available under the [`tests`](./tests/) directory.


# Research Paper

A detailed description of the forward model can be found in the paper:
[A Fast Parallel Solver for the Forward Problem in Electrical Impedance Tomography](https://dx.doi.org/10.1109/TBME.2014.2342280).
