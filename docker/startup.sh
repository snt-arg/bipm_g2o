#!/bin/bash

# Pull or clone bipm_g2o from the master branch
if [ ! -d /workspace/bipm_g2o ]; then
    git clone -b master https://github.com/snt-arg/bipm_g2o.git
else
    echo "Repo already exists, updating..."
    cd /workspace/bipm_g2o
    git checkout master
    git pull origin master
fi

# Enter the project directory
cd /workspace/bipm_g2o

# Start an interactive shell
exec /bin/bash

