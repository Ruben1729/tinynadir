# tinyomm

Heap-free bounded CCSDS/GP decoded-field adapter and KVN subset for conventional SGP4. Required orbit/epoch/B* fields are checked; strict mode requires Earth/TEME/UTC/SGP4 metadata. An explicit CelesTrak-GP-defaults mode accepts their documented omitted redundant fields, while rejecting supplied incompatible metadata. Supports9digit catalog ids, scientific notation, optional canonical numeric unit tags, ISO epochs with optional UTC Z, versions2.0/3.0, duplicate detection, and ephemeris type0. Raw JSON/XML/CSV decoders and full CCSDS structural/schema validation are outside this profile. Max64 fields,8192text bytes and256bytes per KVN line. Stack field descriptors occupy2048bytes on64bit ABI; no heap/globals/I/O. Canonical decoded element type is shared with tinytle; SGP4 initialization does not round-trip through TLE precision. The uniform-day time profile rejects actual UTC leap-second instants. All failed outputs preserved.

## Verification

Independent python-sgp4 OMM initialization/propagation cases are frozen for numeric and9digit ids, multiple orbit regimes, epoch fractions, KVN unit tags and signed time offsets. SGP4 numerical limits1mm/1e-5m/s. Parser tests reject missing/duplicate/incompatible fields and verify explicit defaults and atomicity. Metadata ids do not affect dynamics; reference uses satnum0 for ids its internal Alpha-5 encoding cannot represent. References: https://celestrak.org/NORAD/documentation/gp-data-formats.php and https://github.com/brandon-rhodes/python-sgp4 .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyomm/`. Regression fixtures are checked in and run offline without Python or reference libraries.
