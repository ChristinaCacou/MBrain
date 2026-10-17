# Amendment 2: execution lock supplement

**Date:** 16 October 2026  
**Original public preregistration:** commit `07949d983493186d63d3232a546247476adae4bf`  
**Prospective protocol version:** 1.0, unchanged

## Purpose

This amendment locks the operational materials required to execute the already-registered *Stentor coeruleus* confirmatory study. It does **not** report biological observations and does not change the registered sample target, treatment arms, acquisition sequence, exclusions, endpoints, decision thresholds, or three-gate confirmatory rule.

No confirmatory biological data had been collected by ChatGPT when this execution supplement was generated.

## Randomization and blinding

A blocked randomization schedule has been generated for 192 planned cells: 48 cells on each of four recording days, with 24 blinded treatment-code A and 24 treatment-code B assignments per day. Baseline MECH-first versus ELEC-first order is independently balanced 24:24 within each day.

The treatment mapping itself is stored separately in a private key. This amendment publishes only the SHA-256 hash of that key, which locks its contents without revealing the mapping. The key must not be opened during acquisition, segmentation, exclusions, or blinded QC.

## Execution materials locked by SHA-256

```text
a9ece0f2853d53ebb3dbecb0eb03f7194e6365b9bbefe34a4d52934c9703c039  Stentor_Wet_Lab_Execution_SOP.docx
467a21ff88d1db7db3feff3d1b98c95b461de1fcf4f1308e5a6dac70d2b8b3b7  Stentor_Confirmatory_Execution_Workbook.xlsx
32974690f8ff083eab7ddb6f76ca9c5775121e5df1595aa8ba8f188c84e0d135  Blinded_Randomization_Schedule.csv
b4bb98f109f75f7db0105d84f5a7e60408520b8b933c2c83463ac180d7a60d7a  Calibration_Log_Template.csv
8aeab402e01535bd63aac0f2d0bccc0e7e02d5bfc310ddee68c585c290275d90  Protocol_Deviation_Log.csv
25f6b7451bace7d5c1893b8808c709a19961009d2e54f98f1bb01004f91cf29f  BENCH_RUN_CARD.txt
d64b29d51257bb45024240337e905f5fc691429fd9d59ea93e01fda90a94c770  Run_Locked_Analysis_After_Unblinding.py
7b319877eb749736c51fc302637cd0dc49cdc1b460310a3354b5955d96701151  Prospective_Stentor_Validation_Analysis.py
a6b08df58cc94973659099396cf48db7c01c6f16adb86eabfb2ad4b1a7b57ab2  Stentor_Validation_Data_Template.csv
8637c944ebf82ab68e6012a055b22aae20588974dd52287512f1595fb07d0e8d  PRIVATE_TREATMENT_KEY_DO_NOT_OPEN_UNTIL_DATA_LOCK.txt
```

## Calibration starting points

The execution SOP records published stepper-motor starting points of approximately 1 mm / 0.288 N for Level 4 and approximately 2 mm / 0.557 N for Level 3. These are not replacements for the registered non-confirmatory calibration cohort. The low-force training setting must still be fixed at the laboratory setting giving 70-90% naive contraction; the high-force setting must be the lowest safe setting giving at least 90% contraction.

## Mechanistic caution recorded before acquisition

Historical electrophysiology indicates that electrically evoked action potentials can themselves affect subsequent mechanical habituation. Accordingly, the registered electrical probe should be interpreted as a strong downstream mechanistic control, not as proof that a completely independent sensory channel carries the same state. This interpretive caution does not change the locked confirmatory analysis.

Any future addition of a photic control or other new confirmatory endpoint must be separately time-stamped before the first confirmatory biological observation and cannot be used post hoc to rescue a failed registered gate.

## Final unblinding rule

The private A/B treatment key is revealed only after the full blinded dataset is frozen and QC is complete. `Run_Locked_Analysis_After_Unblinding.py` verifies the registered analysis SHA-256, maps A/B to SHAM/RNAI in a temporary file, invokes the unchanged locked analysis, and writes a final execution receipt.

The inferential ceiling remains unchanged: even if all three registered gates pass, the experiment rejects one explicit mechanistic null and does not establish phenomenal consciousness.
