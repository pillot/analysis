#!/bin/sh

if [ $# = 0 ]; then
echo "provide the file containing datasets"
exit 3
fi

rm -f ds.out

root -l -b > cache.log 2>&1 << EOF
gEnv->SetValue("XSec.GSI.DelegProxy","2");
TProof::Mgr("ppillot@nansafmaster2.in2p3.fr")->SetROOTVersion("la-root");
.q
EOF

for query in `cat $1`; do

  root -l -b -q /Users/pillot/Work/Alice/Macros/Facilities/cacheDatasets.C\(\""$query"\",\""ds.out"\",0\) >> cache.log 2>&1 &

done

