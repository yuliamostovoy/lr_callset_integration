#!/bin/bash
#
set -euxo pipefail

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi

TAG="$1"

cp ../scripts/*.java .
docker build --progress=plain -t us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap .
docker tag us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:${TAG}
docker push us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_ultralong_remap:${TAG}
rm -f *.java

