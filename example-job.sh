#!/bin/bash

#  Job script to execute TLEEDM
#  - Runs the bookkeeper utility
#  - copies input and tleedm application to a work folder
#  - executes tleedm
#  - copies output (as specified by manifest file) back here
#  - deletes work directory, if requested

# define paths for scripts:
BOOKIE=/home/tuwien/viperleed/current_packed/utilities/bookkeeper
TLEEDMSOURCE=/home/tuwien/viperleed/current_packed/viperleed

# define location of work directory
WORK=./work

# define script behaviour
#ALLTENSORS=true    # uncomment to copy all tensors and deltas, instead of only the newest
#RMWORK=true        # uncomment to delete work directory after finishing

#####################################################
#              SCRIPT BELOW IS STATIC               #
#      only touch if you know what you're doing     #
#####################################################

# in case bookkeeper has not been run, run it now (will stop if it has nothing to do)
cp "$BOOKIE" ./bookkeeper
chmod +x bookkeeper
echo "Running bookkeeper..."
./bookkeeper

# create work directory, copy tleedm executable and TensErLEED source code
OUTDIR=`pwd`
mkdir "$WORK" 2> /dev/null

# copy all or only newest Tensors and Deltas to work directory.
if [ "$ALLTENSORS" = "true" ]; then
  cp Tensors "$WORK"/Tensors
  cp Deltas "$WORK"/Deltas
else
  TNUM=""
  mapfile -d $'\0' TENSORS < <(find Tensors -regex 'Tensors/Tensors_[0-9][0-9][0-9].zip' -print0 2> /dev/null)
  if [ "$TENSORS" ]; then
    TNAME=${TENSORS[-1]:(-15)}
    TNUM=${TNAME:(-7):3}
    mkdir "$WORK"/Tensors 2> /dev/null
    cp Tensors/$TNAME "$WORK"/Tensors/$TNAME
  fi
  if [ "$TNUM" ]; then
    DNAME=""
    mapfile -d $'\0' DELTAS < <(find Deltas -regex 'Deltas/Deltas_[0-9][0-9][0-9].zip' -print0 2> /dev/null)
    for f in "${DELTAS[@]}"; do
      if [ ${f:(-7):3} == $TNUM ]; then
        DNAME=${f:(-14)}
      fi
    done
    if [ "$DNAME" ]; then
      mkdir "$WORK"/Deltas 2> /dev/null
      cp Deltas/$DNAME "$WORK"/Deltas/$DNAME
    fi
  fi
fi

# copy input files, go there and execute
cp PARAMETERS VIBROCC IVBEAMS DISPLACEMENTS POSCAR _PHASESHIFTS EXPBEAMS.csv "$WORK" 2> /dev/null
cp -r "$TLEEDMSOURCE"/* "$WORK"
cd "$WORK"
chmod +x tleedm
chmod +x tensorleed/EEASiSSS.x
chmod +x tensorleed/beamgen3.out
./tleedm

# copy the files listed in manifest back to home
xargs -a manifest cp -r -t "$OUTDIR"

# clean up
cd "$OUTDIR"
if [ "$RMWORK" = "true" ]; then
  rm -r $WORK
fi
