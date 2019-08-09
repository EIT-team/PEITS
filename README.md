[![Build Status](https://travis-ci.org/EIT-team/PEITS.svg?branch=master)](https://travis-ci.org/EIT-team/PEITS)

# Install guide

The install guide can be found in the [`INSTALL.md`](./INSTALL.md) file or you
can use the provided docker container.

Build the container available within the repository using:
```bash
docker build -t uclrits/peits:16.04 .
```

or fetch the latest one available as:
```bash
docker fetch uclrits/peits:16.04
```

Then you can run the container as:
```bash
docker run -it --name mypeits -v /full/localpath:/mydata  uclrits/peits:16.04
```
that will login you into the container. If you exit you will stop it.
If you want to restart it and get inside again:

```bash
docker start mypeits
docker attach mypeits
```

To delete the container when you are done with it you
```bash
docker rm mypeits
```
That will delete any content you have written out of the mounted location (_i.e._, `/mydata`).

Alternatively, you can use [docker compose](https://docs.docker.com/compose/) to
build and run the container:

```bash
docker-compose up # to build it
docker-compose run peits # to get inside
docker-compose down # to stop it
```

Change the field `volumes` on `docker-compose.yml` to specify which directory
you want to mount on the docker container.

# Research Paper

A detailed description of the forward model can be found in the paper:
[A Fast Parallel Solver for the Forward Problem in Electrical Impedance Tomography](https://dx.doi.org/10.1109/TBME.2014.2342280).
