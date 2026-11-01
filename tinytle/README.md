# tinytle

Bounded conventional TLE parser supporting numeric and Alpha-5 identifiers (A=10 through Z=33; I/O excluded), checksums, matching catalog ids, the1957/2056 year pivot, Gregorian epoch ranges and SGP4 ephemeris type0. Text buffers need not be NUL terminated; exactly69 characters plus CR/LF are supported. Decoded input fields keep source degrees/revolutions-per-day/B* conventions; these are distinct from SI orbit states. A separate SGP4 adapter initializes from decoded fields and preserves the Vallado conversion order, enabling Alpha-5 or OMM without TLE round trips. Other ephemeris types are explicitly unsupported. The existing tinysgp4 direct parser is unchanged; use this module for expanded inputs. Catalog id metadata up to9digits is retained in decoded fields, not the old internal SGP4 satid. Parser uses no heap; SGP4 initialization has its existing record/core stack costs. Failed outputs/satellites preserved.

## Verification

Frozen independent python-sgp4/Vallado C++ outputs cover numeric/Alpha-5 ids, epoch pivots, near-Earth/deep-space states and positive/negative times. Regression position/velocity tolerances1mm/1e-5m/s. Additional short buffers, checksum, invalid identifiers/epochs and XP-rejection cases test atomicity. Reference: https://github.com/brandon-rhodes/python-sgp4 and https://www.space-track.org/documentation#/tle-alpha5 .

Reference scripts, packages and raw outputs are in ignored `.validation/tinytle/`. Regression fixtures are checked in and run offline without Python or reference libraries.

Epoch-fraction comparisons allow5e-14day (4.32ns): this parser retains the established whole day-of-year conversion order, while python-sgp4 rounds the fraction independently. State comparisons retain their1mm/1e-5m/s limits.

Text parsing, decoded-element initialization and text-to-SGP4 convenience initialization are separate archive objects, so an element-only consumer need not link the unused text parser.
