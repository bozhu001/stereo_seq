# L3 7/7 BCR completion preflight

- Timestamp: 2026-08-28T00:49:35.297779+02:00
- Manifest patients: 7/7.
- Completed and protected from rerun: IPF/FO23-1-06168; SSC/15491/14.
- Remaining: HC/NL-66; HC/NL-72; IPF/FO23-1-06170; SSC/05957/17B; SSC/07998/15A.
- Existing internal candidate R1/R2 and sidecar QNAME manifests: complete for 7/7.
- Original streaming: complete across chunks 9, 10, 11 and 12 (3,939,764,077 pairs); candidate gzip checkpoint present.
- Decision: Case A. Continue from candidate FASTQ; do not stream original FASTQ.
- SSC_15491_14_TRUST4_0003 read-level/VDJ/spatial validation complete: True.
- Expected remaining TRUST4 wall time: approximately 90–110 minutes sequentially at 16 threads, based on the two completed patients and 31.4 GiB compressed R1+R2 input.
- Active scanner/TRUST4/finalizer at preflight: none (checked separately in host/WSL process table).
