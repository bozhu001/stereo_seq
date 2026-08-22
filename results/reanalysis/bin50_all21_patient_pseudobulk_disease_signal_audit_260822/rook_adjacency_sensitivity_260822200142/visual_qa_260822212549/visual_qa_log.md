# Visual QA — rook adjacency sensitivity report

- Report: `../rook_adjacency_sensitivity_report_260822212549.html`
- Report size: 20,871,077 bytes
- SHA-256: `BEAA483BA274A1316EEFC81B8BC0AE0EA7EB05886EFBF756B6EE9B6DBCDE984E`
- HTML structure: opening and closing HTML tags present
- Dynamic content checks: title, all three QC domains, `B_core`, and R 4.5.2
  session information present
- Image packaging: 6 data-URI images embedded; 0 external image links
- Browser: Microsoft Edge headless mode, local `file:` URL, file access enabled
- Primary long-page screenshot: `01_top.png`

## Visual inspection

The primary screenshot was inspected directly. The title, QC table, sensitivity
definition, queen-versus-rook metric table, paired comparison figure, formal
program-QC table, patient disease-test summary, coordinate-validation heading,
and all-21 coordinate contact sheet rendered without missing-image icons or
horizontal clipping.

The report is unusually long because the all-21 contact sheet is embedded at
native detail. Anchor-based screenshots below that figure were not reliable in
Edge's oversized headless viewport: `02_coordinates.png` compressed the
remaining content to the bottom edge and `03_k8_to_tail.png` was blank. These
two files are retained as QA diagnostics and are not counted as successful
visual checks. The remaining embedded images were instead verified by their
presence in the HTML data URIs and by the QMD's 14/14 input-file existence
check.

## Outcome

PASS for report structure, standalone packaging, primary-page layout, key
tables, paired figure, and all-21 contact sheet. No analysis was rerun during
this QA.
