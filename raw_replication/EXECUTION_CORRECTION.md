# alm-3 execution correction

The frozen analysis specification committed before alm-3 numerical outcome extraction explicitly excluded callosotomy conditions. The first session-processing/aggregation implementation did not yet apply an explicit animal-level callosotomy filter. Public CRCNS alm-3 metadata were then inspected without changing any predictor, outcome, time window, inferential unit, null, permutation procedure, or success criterion. Seven animals contain the metadata field `surgicalManipulation = Callosotomy`: ANM270340, ANM278733, ANM291063, ANM291064, ANM297634, ANM298738, and ANM298740.

The definitive aggregate therefore excludes those seven animals exactly as required by the pre-outcome lock. Any unfiltered aggregate from the initial execution is an implementation artifact and is not a reported analysis. This correction implements a frozen eligibility rule; it does not revise the hypothesis or endpoint after seeing the replication result.

Frozen lock commit: 698515d08cceaf82f959f53df3c2a713adc81fe4
Corrected aggregator commit: 68b12ae23da7b38c895ac49cef9c224f5816c5e8
Metadata-inspection workflow run: 38272786434
