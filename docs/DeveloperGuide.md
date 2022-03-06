# Developer Guide

The project team inends to maintain a containerized workspace so as to better facilitate contributions from contributors who might otherwise favor disparate development environments.

## Standing up a local workspace
The [workspace directory](../workspace) contains a Dockerfile which will need to be built locally.

```
docker build -t nv_workspace [project root]/workspace/
```
Once the build is complete you should be able to see the container image
```
$ docker image ls | grep nv_workspace
nv_workspace     latest    bbcac7be2aba   7 minutes ago    1.01GB
```

To start the workplace simply run the container
```
$ docker run -it nv_workspace
virtualenvwrapper.user_scripts creating /root/.virtualenvs/premkproject
virtualenvwrapper.user_scripts creating /root/.virtualenvs/postmkproject
virtualenvwrapper.user_scripts creating /root/.virtualenvs/initialize
...
(nv) root@c5a9028723ef:/#
```
This will start the container and activate the "nv" python virtual environment.