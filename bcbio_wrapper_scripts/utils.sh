#!/bin/bash

# TODO 1    CHECK AND SET UP ENVIRONMENT VARIABLES
## should check env variables

# TODO 2    PASRSE THE INPUT YAML CONFIG FILE

## Function for reading from a simple YAML file, adapted from: https://stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script/21189044#21189044
## Modified the printf line in the awk script code below to add the word "export ", such that it creates environment variables accessible from the rest of the analysis scripts
function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("export %s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

echo "using configuration file: " $1


## read the config file and create the variables using the "bcbio_" prefix
eval $(parse_yaml $1 "bcbio_")

echo "" 
echo "## printing the whole contents of the file:"

parse_yaml $1 "bcbio_"

echo "## set the \$PATH variable like this:"
echo "export PATH=${bcbio_install_path:?}/tools/bin:${bcbio_install_path:?}/anaconda/bin:${bcbio_install_path:?}/anaconda/envs/bcbiovm/bin:${bcbio_install_path:?}/extra3/bin:${bcbio_install_path:?}/extra2/bin:${PATH}"
echo "" 
# TODO 4
# TODO 5
