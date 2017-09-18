#!/usr/bin/env bash

java -cp $(cygpath /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) edu.ucsd.msjava.ui.MzIDToTsv "$@"
