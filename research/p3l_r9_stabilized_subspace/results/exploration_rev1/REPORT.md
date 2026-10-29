# E3_REV1 weak-ray repair (exposed DEV only)

This separate preregistered development revision changes the failed fixed-coordinate radius projection to original R7 radial retraction, preserving the fixed weak ray in B while allowing its amplitude to change. Original E3 remains immutable.

| arm | r | available | qualified | pair<=1 | median pair | op RHS | audit RHS | op wall(s) | one-state-only | max ray defect |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|B_RAY|0.01|8|8|8|0.1143444874392381|696|32|1.6985929997172207|0|5.855417636698539e-16|
|C_RAY|0.01|8|8|8|0.11379567898117043|696|32|1.503564499784261|0|0.1409866147552032|
|B_RAY|0.02|8|8|7|0.23866129185363272|704|64|1.4055763757787645|0|2.43586500355971e-16|
|C_RAY|0.02|8|8|7|0.2493527424462399|704|64|1.5177324172109365|0|0.21953175880790562|

Additional actual physics: 2992 RHS, 94 states/LUs, 0 solver failures. Aggregate exploration RHS=16968.

All operational and audit costs include failed/rejected attempts. Current selected endpoint J diagnostics occur after operational selection timing. A and Frozen-J remain original-run references; their wall times are from a separate timing window. No fresh confirmation, posterior, or continuous-fiber claim.

Reference endpoints (original timing window for A/Frozen/U):

|arm|r|qualified|pair<=1|median pair|operational RHS|timing window|
|---|---:|---:|---:|---:|---:|---|
|A_NATIVE_FULLJ|0.01|8|8|0.13276662081042567|696|ORIGINAL_E3|
|F_FROZENJ|0.01|8|8|0.10844494267497981|192|ORIGINAL_E3|
|U_UNCORRECTED|0.01|8|8|0.1686984577640922|64|ORIGINAL_E3|
|B_RAY|0.01|8|8|0.1143444874392381|696|REV1|
|C_RAY|0.01|8|8|0.11379567898117043|696|REV1|
|A_NATIVE_FULLJ|0.02|8|7|0.28691756857431383|696|ORIGINAL_E3|
|F_FROZENJ|0.02|8|6|0.21667018781149971|192|ORIGINAL_E3|
|U_UNCORRECTED|0.02|6|6|0.5706431367373985|64|ORIGINAL_E3|
|B_RAY|0.02|8|7|0.23866129185363272|704|REV1|
|C_RAY|0.02|8|7|0.2493527424462399|704|REV1|

Within REV1 paired timing:
- r=0.01: fixed B relative wall saving=-12.971%; RHS saving=0.000%; 20% efficiency gate=False; endpoints differ in 8/8.
- r=0.02: fixed B relative wall saving=7.390%; RHS saving=0.000%; 20% efficiency gate=False; endpoints differ in 8/8.

This is a development mechanism repair. Original negative results remain part of the evidence; qualifying endpoints do not establish posterior or continuous-fiber validity.
