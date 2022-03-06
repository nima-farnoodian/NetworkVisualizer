# Developer Guide

The project team inends to maintain a containerized workspace so as to better facilitate contributions from contributors who might otherwise favor disparate development environments.

## Standing up a local workspace
**Prerequisites:**
[docker](https://docs.docker.com/get-docker/)

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
$ docker run -it -p 8888:8888 -v $(pwd):/nvworkspace nv_workspace
virtualenvwrapper.user_scripts creating /root/.virtualenvs/premkproject
virtualenvwrapper.user_scripts creating /root/.virtualenvs/postmkproject
virtualenvwrapper.user_scripts creating /root/.virtualenvs/initialize
...
(nv) root@708b6fefd888:/nvworkspace#
```
This will start the container and activate the "nv" python virtual environment.

To prep the workspace run the `prep-workspace.sh` script.
```
(nv) root@708b6fefd888:/nvworkspace# ./workspace/prep-workspace.sh
```

Once the workspace is prepped you may run jupyter lab
```
jupyter lab --allow-root --ip=0.0.0.0 &
...
    To access the server, open this file in a browser:
        file:///root/.local/share/jupyter/runtime/jpserver-192-open.html
    Or copy and paste one of these URLs:
        http://b98bf8b21709:8889/lab?token=a00d051b095d8f2656dea4c79bb222c575f592eb1351afee
     or http://127.0.0.1:8889/lab?token=a00d051b095d8f2656dea4c79bb222c575f592eb1351afee
```