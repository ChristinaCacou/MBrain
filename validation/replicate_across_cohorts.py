#!/usr/bin/env python3
"""Post hoc two-cohort directional replication screen. Entire primary DLPFC discovery set.
Primary source: GSE174409 edgeR QL (25,122 genes).
Independent BA9 source: GSE182321 edgeR QL (19,327 genes).
Gene names: public GENCODE v48 GTF (mapping ambiguities excluded).
No feature filtering based on the external cohort's p-values or effect sizes.
"""
import csv,gzip,re,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

out=Path("replication"); out.mkdir(exist_ok=True)
p=pd.read_csv("primary/results/DLPFC_additive_OUD_edgeR_QL.csv").rename(columns={"gene_id":"ensembl","PValue":"p_primary","FDR":"q_primary","logFC":"fc_primary"})
e=pd.read_csv("external/independent_outputs/GSE182321_edgeR_full.csv").rename(columns={"gene":"symbol","PValue":"p_external","FDR":"q_external","logFC":"fc_external"})
assert len(p)==25122 and len(e)==19327
assert p.ensembl.is_unique and e.symbol.is_unique
discover=p.loc[p.q_primary<.05].copy()
assert len(discover)==1288
m={}
by_sym=defaultdict(set)
with gzip.open("gencode.gtf.gz","rt") as f:
    for l in f:
        if l.startswith("#"):continue
        cols=l.split("\t")
        if len(cols)<9 or cols[2]!="gene":continue
        a=re.search(r'gene_id "([^"]+)"',cols[8])
        b=re.search(r'gene_name "([^"]+)"',cols[8])
        if not a or not b:continue
        id=a.group(1).split(".")[0]
        sym=b.group(1)
        m.setdefault(id,set()).add(sym)
        by_sym[sym].add(id)
gene_map={k:next(iter(v)) for k,v in m.items() if len(v)==1 and len(by_sym[next(iter(v))])==1}
discover["symbol"]=discover.ensembl.map(gene_map)
discover["mapping_status"]=np.where(discover.symbol.isna(),"not_unique_or_missing_annotation","mapped")
c=discover.merge(e[["symbol","p_external","q_external","fc_external"]],on="symbol",how="left",validate="many_to_one")
c["testable"]=c.p_external.notna()
test=c[c.testable].copy()
# BH across the entire available, primary-discovered confirmatory family,
# NOT across genes chosen using follow-up external significance.
test["q_followup_BH"]=multipletests(test.p_external.to_numpy(),method="fdr_bh")[1]
test["sign_concordant"]=np.sign(test.fc_primary)==np.sign(test.fc_external)
test["replicated"]=test.sign_concordant & (test.q_followup_BH<.05)
c=c.merge(test[["ensembl","q_followup_BH","sign_concordant","replicated"]],on="ensembl",how="left",validate="one_to_one")
test=test.sort_values(["q_followup_BH","p_external"])
c.to_csv(out/"primary_1288_all_candidates.csv",index=False)
test.to_csv(out/"eligible_independent_replication_tests.csv",index=False)
rep=test.loc[test.replicated].copy()
rep.to_csv(out/"replicated_direction_concordant_BH05.csv",index=False)
egr=test.loc[test.symbol=="EGR1"]
negative_controls = test.loc[(test.q_followup_BH<.05)&~test.sign_concordant]
summary={
"primary_fdr05_genes":int(len(discover)),
"primary_unmapped_or_nonunique":int(discover.symbol.isna().sum()),
"not_testable_external":int((~c.testable).sum()),
"external_testable_candidates":int(len(test)),
"followup_fdr05_genes_any_direction":int((test.q_followup_BH<.05).sum()),
"followup_fdr05_direction_concordant":int(len(rep)),
"followup_fdr05_opposite_direction":int(len(negative_controls)),
"EGR1_in_primary_discoveries":bool((discover.ensembl=="ENSG00000120738").any()),
"EGR1_mapped_external":bool(len(egr)==1),
"EGR1_external_p":float(egr.p_external.iloc[0]) if len(egr) else None,
"EGR1_followup_BH":float(egr.q_followup_BH.iloc[0]) if len(egr) else None,
"EGR1_direction_concordant":bool(egr.sign_concordant.iloc[0]) if len(egr) else None,
"EGR1_replicated":bool(egr.replicated.iloc[0]) if len(egr) else None,
"primary_discovery_source":"GSE174409 official edgeR 4.0.16 log2 FC and BH FDR",
"external_test_source":"GSE182321 official edgeR 4.0.16 log2 FC and nominal P",
"mapping_source":"GENCODE human v48",
"limitation":"Procedure devised after viewing EGR1 results, so post hoc. Cohort confounders and source-published nomination of EGR1 remain.",
}
import json
(out/"REPLICATION_SUMMARY.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print("Top replicated:",rep[["ensembl","symbol","fc_primary","q_primary","fc_external","p_external","q_followup_BH"]].head(20).to_string(index=False))
assert test["q_followup_BH"].between(0,1).all()
assert len(test)==test["ensembl"].nunique()
