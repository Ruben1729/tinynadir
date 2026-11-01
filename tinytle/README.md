# tinytle

A heap-free bounded conventional TLE parser and SGP4 initialization adapter. Text buffers need not be NUL terminated: exactly 69 columns plus up to two CR/LF bytes are accepted. Dynamical fields retain their source degrees, revolutions/day and inverse-Earth-radius BSTAR units; propagated orbit states use the existing `tinysgp4` TEME km/km/s interface.

The text profile supports numeric catalog IDs 0–99999 padded with leading zeros or spaces, and Alpha-5 IDs A0000–Z9999 (A=10 through Z=33; I/O excluded). Checksums, matching decoded catalog IDs, dynamical field ranges and ephemeris type 0 are checked. A completely blank BSTAR field means zero. Epochs use the 1957/2056 year pivot; day 000, including its fractional day, denotes the preceding December 31. Classification, launch designator, element/revolution numbers and mean-motion derivative fields are ignored by this adapter. The mean-motion derivatives do not enter conventional SGP4 dynamics.

Decoded elements can also initialize SGP4 directly, using WGS72 improved mode and Vallado's unit-conversion evaluation order. This supports Alpha-5 and OMM without a TLE precision round trip. Decoded initialization accepts uniform-day UTC epochs from 1600 through 9999 and preserves catalog metadata up to nine digits in the decoded fields. The older internal SGP4 satid is not used for these metadata IDs. An initialization is published only after propagation at the epoch produces finite position and velocity; failed parsing or initialization preserves the caller's existing output/satellite.

Text parsing, decoded-element initialization and text-to-SGP4 convenience initialization are separate archive objects, so an element-only consumer need not link the unused text parser. No heap, mutable globals or runtime I/O are used. SGP4 initialization retains its record/core stack costs and performs an additional epoch-state validation.

## Verification

Frozen independent python-sgp4 2.27 / bundled Vallado C++ outputs cover numeric/Alpha-5 IDs, epoch pivots, near-Earth/deep-space states and positive/negative times. Four additional independent legacy-format fixtures cover blank BSTAR, leading-space numeric IDs, epoch day000 and day000.25, checking decoded epochs and propagation at −60/0/+60 minutes. Position/velocity tolerances are 1 mm / 1e-5 m/s. Epoch-fraction comparisons allow 5e-14 day (4.32 ns): this parser retains the established whole day-of-year conversion order, while python-sgp4 rounds the fraction independently.

Additional tests check short buffers, checksums, invalid IDs/epochs, unsupported ephemeris types, atomic failures and rejection of extreme finite mean motion that yields a nonfinite epoch state. Reference scripts, packages and raw outputs are in ignored `.validation/tinytle/`; checked-in numeric fixtures run offline without Python or reference libraries.

References: [python-sgp4](https://github.com/brandon-rhodes/python-sgp4), [CelesTrak TLE conventions](https://celestrak.org/columns/v04n03/), [Space-Track Alpha-5](https://www.space-track.org/documentation#/tle-alpha5).
