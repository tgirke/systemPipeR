#!/bin/bash
set -e
# git_link is the http git link of a repo
git_link=${1}
# update_dir is what is the folder path in website you want to copy to
update_dir=${2}
down_folder="/tmp/${update_dir}"
folder_base=$(basename ${update_dir})

# clean temp
echo "Clean dwonload folder ${down_folder}"
rm -rf ${down_folder}
# clone repo 
echo "Downloading ${update_dir} to \"${down_folder}\""
git clone ${git_link} ${down_folder}

# update
echo "Clean update folder ${update_dir}"
rm -rf ${update_dir}
mkdir -p ${update_dir}
# copy 
echo "Update files"
cp -r ${down_folder}/docopt.R/* ${update_dir}
echo "done"


