# Make.sh = update Makefile.lib, Makefile.shlib, Makefile.list
#           or style_*.h files
# Syntax: sh Make.sh style
#         sh Make.sh Makefile.lib
#         sh Make.sh Makefile.shlib
#         sh Make.sh Makefile.list

# function to create one style_*.h file
# must whack *.d files that depend on style_*.h file,
# else Make will not recreate them

style () {
  # modified C.K. create version_liggghts.h
  builddate=`date +%Y-%m-%d-%H:%M:%S`
  wai=`whoami`
  vers=`cat version_liggghts.txt`
  bra=`cat version_liggghts_branch.txt`

  if [ -d .git ]; then
    githash=`git log -1 --format="%H"`
    echo "#define LIGGGHTS_VERSION \"$bra $vers, compiled $builddate by $wai, git commit $githash\"" > version_liggghts.h
  elif [ -d ../.git ]; then
    cd ..    
    githash=`git log -1 --format="%H"`
    cd src
    echo "#define LIGGGHTS_VERSION \"$bra $vers, compiled $builddate by $wai, git commit $githash\"" > version_liggghts.h
  else
    echo "#define LIGGGHTS_VERSION \"$bra $vers, compiled $builddate by $wai, git commit unknown\"" > version_liggghts.h
  fi;

  list=`grep -sl $1 $2*.h`
  if (test -e style_$3.tmp) then
    rm -f style_$3.tmp
  fi
  for file in $list; do
    qfile="\"$file\""
    echo "#include $qfile" >> style_$3.tmp
  done
  if (test ! -e style_$3.tmp) then
    if (test ! -e style_$3.h) then
      touch style_$3.h
    elif (test "`cat style_$3.h`" != "") then
    rm -f style_$3.h
    touch style_$3.h
      rm -f Obj_*/$4.d
      if (test $5) then
        rm -f Obj_*/$5.d
      fi
      rm -f Obj_*/lammps.d
    fi
  elif (test ! -e style_$3.h) then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    if (test $5) then
      rm -f Obj_*/$5.d
    fi
    rm -f Obj_*/lammps.d
  elif (test "`diff --brief style_$3.h style_$3.tmp`" != "") then
    mv style_$3.tmp style_$3.h
    rm -f Obj_*/$4.d
    if (test $5) then
      rm -f Obj_*/$5.d
    fi
    rm -f Obj_*/lammps.d
  else
    rm -f style_$3.tmp
  fi
}

# create individual style files
# called by "make machine"
# col 1 = string to search for
# col 2 = search in *.h files starting with this name
# col 3 = prefix of style file
# col 4

if (test "$1" = "style") then
#      | search string         | .h file name      | style file prefix |                   |
  style ANGLE_CLASS             angle_              angle               force
  style ATOM_CLASS              atom_vec_           atom                atom                atom_vec_hybrid
  style BODY_CLASS              body_               body                atom_vec_body
  style BOND_CLASS              bond_               bond                force
  style COMMAND_CLASS           ""                  command             input
  style COMPUTE_CLASS           compute_            compute             modify              modify_cuda
  style DIHEDRAL_CLASS          dihedral_           dihedral            force
  style DUMP_CLASS              dump_               dump                output
  style FIX_CLASS               fix_                fix                 modify
  style IMPROPER_CLASS          improper_           improper            force
  style INTEGRATE_CLASS         ""                  integrate           update
  style KSPACE_CLASS            ""                  kspace              force
  style MINIMIZE_CLASS          min_                minimize            update
  style PAIR_CLASS              pair_               pair                force
  style SURFACE_MODEL           surface_model_      surface_model       force
  style NORMAL_MODEL            normal_model_       normal_model        force
  style TANGENTIAL_MODEL        tangential_model_   tangential_model    force
  style COHESION_MODEL          cohesion_model_     cohesion_model      force
  style ROLLING_MODEL           rolling_model_      rolling_model       force
  style READER_CLASS            reader_             reader              read_dump
  style REGION_CLASS            region_             region              domain
  style CFD_DATACOUPLING_CLASS  cfd_datacoupling_   cfd_datacoupling    fix_cfd_coupling
  style CFD_REGIONMODEL_CLASS   cfd_regionmodel_    cfd_regionmodel     fix_cfd_coupling
  style LB_CLASS                ""                  lb
  style SPH_KERNEL_CLASS        sph_kernel_         sph_kernel          pair_sph-fix_sph
  style MESHMODULE_CLASS        mesh_module_        mesh_module         fix_mesh_surface
  style MESHMOVER_CLASS         mesh_mover_         mesh_mover          mesh_mover
elif (test "$1" = "models" -o "$1" = "models_full") then
  sed_ex="sed -E" # BSD sed
  sed --version | grep -i gnu > /dev/null 2>&1
  [ $? -eq 0 ] && sed_ex="sed -r" # GNU sed

  surface_models=`grep -s -E '^SURFACE_MODEL\(' surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  surface_model_ids=`grep -s -E '^SURFACE_MODEL\(' surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/'`
  normal_models=`grep -s -E '^NORMAL_MODEL\(' normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  normal_model_ids=`grep -s -E '^NORMAL_MODEL\(' normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/'`
  tangential_models=`grep -s -E '^TANGENTIAL_MODEL\(' tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  tangential_model_ids=`grep -s -E '^TANGENTIAL_MODEL\(' tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/'`
  cohesion_models=`grep -s -E '^COHESION_MODEL\(' cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  cohesion_model_ids=`grep -s -E '^COHESION_MODEL\(' cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/'`
  rolling_models=`grep -s -E '^ROLLING_MODEL\(' rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/'`
  rolling_model_ids=`grep -s -E '^ROLLING_MODEL\(' rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/'`

  #echo $surface_model_ids
  #echo $normal_model_ids
  #echo $tangential_model_ids
  #echo $cohesion_model_ids
  #echo $rolling_model_ids

  # check for duplicate constants
  sm_duplicates=`grep -s -E '^SURFACE_MODEL\(' surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  nm_duplicates=`grep -s -E '^NORMAL_MODEL\(' normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  tm_duplicates=`grep -s -E '^TANGENTIAL_MODEL\(' tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  cm_duplicates=`grep -s -E '^COHESION_MODEL\(' cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`
  rm_duplicates=`grep -s -E '^ROLLING_MODEL\(' rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\3/' | sort | uniq -d`

  if [ -n "$sm_duplicates" ]; then echo "ERROR: duplicate surface model identifiers:"; echo $sm_duplicates; exit -1; fi
  if [ -n "$nm_duplicates" ]; then echo "ERROR: duplicate normal model identifiers:"; echo $nm_duplicates; exit -1; fi
  if [ -n "$tm_duplicates" ]; then echo "ERROR: duplicate tangential model identifiers:"; echo $tm_duplicates; exit -1; fi
  if [ -n "$cm_duplicates" ]; then echo "ERROR: duplicate cohesion model identifiers:"; echo $cm_duplicates; exit -1; fi
  if [ -n "$rm_duplicates" ]; then echo "ERROR: duplicate rolling model identifiers:"; echo $rm_duplicates; exit -1; fi

  stylefile=style_contact_model.h
  filteredfile=style_contact_model_filtered.tmp
  minfile=style_contact_model_minimal.whitelist
  whiteLfile=style_contact_model.whitelist
  whiteLuserfile=style_contact_model_user.whitelist
  whiteLautofile=style_contact_model_autoExamples.whitelist

  tangential_models="TANGENTIAL_OFF $tangential_models"
  cohesion_models="COHESION_OFF $cohesion_models"
  rolling_models="ROLLING_OFF $rolling_models"

  # if whitelist exists, copy it
  # otherwise use the default minimal one
  if (test -e $whiteLfile) then
    cp $whiteLfile $filteredfile
  else
    cp $minfile $filteredfile
  fi

## add merging of user, autoExamples and filteredFile

  if (test -e $whiteLuserfile) then
    while read -r line; do
      if [ -z "$(grep "${line}" "${filteredfile}" )" ] ; then
        echo "${line}" >> "${filteredfile}"
#        echo "${line}"
      fi
    done < $whiteLuserfile
  fi
  if (test -e $whiteLautofile) then
    while read -r line; do
      if [ -z "$(grep "${line}" "${filteredfile}" )" ] ; then
        echo "${line}" >> "${filteredfile}"
#        echo "${line}"
      fi
    done < $whiteLautofile
  fi

  if (test ! -e $filteredfile) then
    rm -f $stylefile
    touch $stylefile
  elif (test ! -e $stylefile) then
    mv $filteredfile $stylefile
    rm -f Obj_*/force.d
    rm -f Obj_*/modify.d
    rm -f Obj_*/lammps.d
  elif (test "`diff --brief $stylefile $filteredfile`" != "") then
    mv $filteredfile $stylefile
    rm -f Obj_*/force.d
    rm -f Obj_*/modify.d
    rm -f Obj_*/lammps.d
  else
    rm -f $filteredfile
  fi

  echo "Creating list of contact models completed."

fi
