#!/bin/bash

this_path=$(pwd)
echo "Adding the path: "$this_path" to the PYTHONPATH"
export PYTHONPATH+=:$this_path
