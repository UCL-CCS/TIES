#!/bin/sh
cwd=$(pwd)

script_path=$(realpath $0)
scipt_dir=$(dirname $script_path)
echo mat $scipt_dir
# set the current working directory

pushd
cd $scipt_dir

vmd -e vis.vmd

popd # for some reason doesn't kick in
cd $cwd