==== python2 =====
   ## create export for extra 2
   # echo "export PATH=${bcbio_install_path%?}/extra_conda2/bin:\$PATH" >> ~/.bashrc
   # source ~/.bashrc
   export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${bcbio_install_path%?}/extra_conda2/bin:${PATH}

else 
   echo " --- [$(date +"%F %R")] Check if the PATH is set correctly for extra3 environment."
   if [[ ":$PATH:" == *"${bcbio_install_path%?}/extra_conda2/bin"* ]]; then
      echo " --- [$(date +"%F %R")] The path is set correctly."
   else
      echo " --- [$(date +"%F %R")] PATH does not contain then path to the bin directory of the Python2 environment."
      echo " --- [$(date +"%F %R")] Create exports for the Python3 environment path."
      ## create export for extra 2
      # echo "export PATH=${bcbio_install_path%?}/extra_conda2/bin:\$PATH" >> ~/.bashrc
      # source ~/.bashrc

   fi

==== python3 =====
   # echo "export PATH=${bcbio_install_path%?}/extra3/bin:\$PATH" >> ~/.bashrc
   # source ~/.bashrc
   export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${PATH}


else 
   echo " --- [$(date +"%F %R")] Check if the PATH is set correctly for extra3 environment."
   if [[ ":$PATH:" == *"${bcbio_install_path%?}/extra3/bin"* ]]; then
      echo " --- [$(date +"%F %R")] The path is set correctly."
   else
      echo " --- [$(date +"%F %R")] PATH does not contain then path to the bin directory of the Python3 environment."
      echo " --- [$(date +"%F %R")] Create exports for the Python3 environment path."
      ## create export for extra 3
      # echo "export PATH=${bcbio_install_path%?}/extra3/bin:\$PATH" >> ~/.bashrc
      # source ~/.bashrc
      export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${bcbio_install_path%?}/extra3/bin:${PATH}

-======bcbio already installed =========
   ## check if exports were made for that installation if not make them
   echo " --- [$(date +"%F %R")] Check if the PATH is set correctly for bcbio_nexgetn installation."
   if [[ ":$PATH:" == *"${bcbio_install_path%?}"* ]]; then
      echo " --- [$(date +"%F %R")] The path is set correctly."
   else
      echo " --- [$(date +"%F %R")] PATH does not contain then path to bcbio_nextgen installation."
      echo " --- [$(date +"%F %R")] Create exports for the bcbio_nextgen path."
      ## create export for bcbio path
      # echo "export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:\$PATH" >> ~/.bashrc
      # source ~/.bashrc
      export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${PATH}
   fi


=== bcbio install =====
## create symlink for bcbio
# echo "export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:\$PATH" >> ~/.bashrc
# source ~/.bashrc
export PATH=${bcbio_install_path%?}/anaconda/bin:${bcbio_install_path%?}/tools/bin:${PATH}

