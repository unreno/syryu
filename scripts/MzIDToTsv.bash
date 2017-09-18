#!/usr/bin/env bash

java -cp $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) edu.ucsd.msjava.ui.MzIDToTsv "$@"
