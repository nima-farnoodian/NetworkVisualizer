#!/bin/bash
cat >> ~/.bashrc <<EOF
export WORKON_HOME=$HOME/.virtualenvs
export PROJECT_HOME=$HOME/dev
source /usr/local/bin/virtualenvwrapper.sh
mkvirtualenv -p python3.9 nv
EOF