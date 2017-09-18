#!/usr/bin/env bash

echo "Calling ..."
echo java -cp $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) edu.ucsd.msjava.ui.MzIDToTsv "$@"
echo

java -cp $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) edu.ucsd.msjava.ui.MzIDToTsv "$@"
