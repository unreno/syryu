#!/usr/bin/env bash

echo "Calling ..."
echo java -Xmx3500M -jar $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) "$@"
echo

java -Xmx3500M -jar $(cygpath -w /cygdrive/c/ryulab/MSGFPlus/MSGFPlus.jar) "$@"
