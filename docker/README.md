# Docker container

To download and mount the Docker image, first download and install the relevant version of [Docker](https://docs.docker.com/get-docker/) for your OS, and then run the following in a command prompt:

```
docker pull qdercon/pstpipeline:v0.2.0
```

The image includes everything required to run the Jupyter notebooks both locally or on a cluster/cloud server in a containerised environment (i.e., local package installs etc. will not be affected). More specifically, it is a Linux environment containing:

* R v4.1.2 plus all package dependencies (see "DESCRIPTION" file for full details)
* Python v3.9.5 plus all dependencies, including rpy2 for running R code in Jupyter notebooks
* JupyterLab
* CmdStan v2.28.1
* Jupyter notebooks and required raw data

The image can also be built locally from the included Dockerfile. For example, on Windows, to clone the repo, extract the relevant data & notebooks, and build the Docker image, you could run the following:

```
git clone https://github.com/qdercon/pstpipeline
cd pstpipeline
cp -a data-raw docker/data-raw
cp -a notebooks docker/notebooks
cd docker
docker build -t pstpipeline-docker .
```

Once downloaded or built, to mount the image, run the following in a command prompt:

```
docker run -it --rm -p 8888:8888 -v <:/path/to/folder>:/root/<mount_folder_name>/ pstpipeline-docker
```

The -v flag and the path that follows is optional; this allows you to "mount" a folder on the disk to enable notebooks and model outputs to be saved locally. The command will output a link beginning with ```http//:127.0.0.1:8888/lab?token=``` which can be copied and pasted into a browser to open JupyterLab.
