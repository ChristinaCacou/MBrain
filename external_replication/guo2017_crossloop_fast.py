"""Parallel implementation of the already-frozen GUO2017_CROSSLOOP_LOCK.md.

This changes no scientific definition. The eligible asset IDs were established by the
outcome-blind structural inventories before numerical outcome extraction.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import numpy as np
import pandas as pd
from external_replication.guo2017_crossloop_analysis import (
    OUT, process_asset, list_assets, pooled_stat, perm_primary, boot_ci, exact_dir, loo
)

ELIGIBLE = set(['52011d93-3d4e-40d1-bbd5-a86e42d211e1', 'fd53079b-e639-4436-829b-cb3aac34fe71', '8dc414b5-f735-43ef-847e-eef252f6cc55', '61d441af-7bf1-4603-9d42-491e255da8f9', '8ce1a50f-11bd-4a75-a510-64c3f32bb529', 'bc8f7c92-fd60-43b8-bc5b-fcc796f59a5c', 'b24c782d-1734-4b28-a536-59c08feeee7b', '5c6dc876-22f3-4773-8fa3-716f673ecc59', 'e3eafe27-c7e9-45de-aa58-4e95eaf7f2cf', 'de54dfcd-18bb-43a0-b54e-0e429f7b972b', '29cf87e2-e0d1-41e4-846e-9c6f05b5ed67', '29bf0bb3-b18e-45c9-97b4-adc481b53dcb', 'eeae3679-f20c-487d-9306-561936850e6e', 'e7eda98a-c08c-4943-b221-f032433b0a03', 'b69b65f9-a02e-4745-9dac-b5a0b4459871', '2b56b171-5bfe-49f4-a756-7ee77123af52', '4061eece-6fc4-46f1-a0cf-1ef260104c7e', '29d62c6f-93c2-4311-b6cb-e5d06b188844', '3c85c4bc-d52e-4249-8ff1-227eb1864669', '944547c0-9fe3-4dd1-92c4-e2152b507c7e', '02ed0a10-2d6d-4c38-882c-b0b54d6dbe22', 'ba0f2368-7948-4e94-8260-5cc538174b41', '62ec587b-99af-4001-a265-1097089a73ff', 'e2c8bbeb-7244-4bd4-9020-613b071be8a7', 'f58ecdd7-04aa-43a1-b82d-599bd7b9aa7b', '754d3364-4944-477c-b626-c94ab5d119d0', '93d145fd-8d10-439e-b588-3a941be2a793', '43ca97bd-b382-49b5-9aaa-55a466591794', 'e4aa2177-04e5-44ff-af38-0a58b098851c', '26fa5c75-6891-45ea-8019-ea413da61438', '5592d3cb-205e-4d6b-90db-6535340f8f36', '0ed2b558-2008-428d-809d-e4c735314389', '305a64c7-1eb2-436f-9167-8af74c780974', 'dd648e62-d70c-4d32-b77f-f88ac16b0cf6', 'd3748ac3-ec4b-4380-8d07-dcf452616250', '4135ca35-fcdd-40e3-9b24-d912dff9f2fe', '6da90d19-c15b-4411-b56f-c29b003806c7', '26e9ac40-ed73-4efc-9a6d-ea0d25af4a59', '380e2a16-a8d8-4595-b7d2-e0b42595216d', 'add32df8-7f38-476c-acd0-3b1abf510e7f', 'd5c2f043-4571-42da-ab51-54577759179e', '57313ac8-eeb3-4dcc-9ae5-315880069570', '86a9a2c9-429f-4bf3-bd63-d925c29dbf45', '092550fe-91d3-4585-a097-8517a8d7657e', '2f420877-b560-4677-8a1b-78b24440269b', '95c54b63-55ac-48bd-baba-362d9af5178b', '4d946ae0-6f76-42cd-9cff-ea211c5db1da', '10b5c3c8-05cc-4dda-8df2-8186be74b764', '51eb2318-4855-4115-85d5-26038cda6d0e', 'ec181b81-c2c2-4e2a-ac6d-f907c4909e73', 'b665fdf6-dc89-472c-9316-aa86d7ce2325', '7e2af79b-0290-4cb5-9d14-457fe8da3759', 'fd13975a-4fbd-4da0-85f2-6d1ff7222877', 'ee8800ff-5f80-4af1-8c4c-e856167dd35c', 'f3583340-bfbf-4ac6-acfb-72bf42588699', '81f3f22d-35c3-44b1-bda5-446e4cb6b401', '971b4c2a-84cc-493a-87e5-63cf6839f38c', 'acd40a40-d63f-4bd4-ae05-cb7154f9de99', '475e61ba-ca60-4c7a-80cc-1d954ce5d6af', '4fb553c6-1151-4a57-85a5-d07ff574c59f', '85a9980c-6c21-4499-84cd-73bc1277fd1b', '35b2e771-1e75-46cc-9a52-abc4100e3c88', '78307e23-96ac-4167-aa0d-b64d7628c846', '386482f1-959c-4442-a01e-6b1dab9bd6eb', 'd4816d4b-10f8-410c-837a-4d47ea68e08a', 'aea47d9d-5331-4ac9-a6a9-87435ac02c0b', '623bfb46-404b-44d7-a0df-5aff76c1f3b3', '19a1fc05-60ee-4ea7-90c6-f830f4acfd33', '0590607a-e6e5-4264-a765-df6f60d46510', 'fd52ae02-ba0e-463a-9d5f-ab04486ad127', '73aecd59-8d8a-491d-923e-6d3c4654deab', '3fd30419-474a-4e98-9e14-e483d0ba9b77', 'eedfe77c-357b-469a-92b1-ac10ed6c5ed5', '694e0c05-bb88-47b7-8626-76070dc2aee3', '59bdaa7c-8b7d-425d-abe9-0adaa5db2110', 'acd95096-4ef7-4d54-947b-e4e7289fb9d2', '306de10c-65ab-4dd9-a2c4-00c9ed8646f2', '87340cdb-9869-4ec5-aa2e-e5f4912b6c4e', 'bdb53aa7-03ba-45fb-a212-739aeea56ccf'])

def main():
    assets=[a for a in list_assets() if a['asset_id'] in ELIGIBLE]
    rows=[]
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs={ex.submit(process_asset,a):a['asset_id'] for a in assets}
        for f in as_completed(futs):
            r=f.result(); rows.append(r)
            print(r.get('subject'),r.get('recording_region'),r.get('late_auc'),r.get('impairment'),r.get('error'))
    sess=pd.DataFrame(rows);sess.to_csv(OUT/'session_results.csv',index=False)
    eligible=sess[np.isfinite(pd.to_numeric(sess.late_auc,errors='coerce'))&np.isfinite(pd.to_numeric(sess.impairment,errors='coerce'))].copy()
    mouse=eligible.groupby(['subject','recording_region','target_site'],as_index=False).agg(late_auc=('late_auc','mean'),early_auc=('early_auc','mean'),impairment=('impairment','mean'),n_sessions=('asset_id','size'),mean_control_accuracy=('control_accuracy','mean'),mean_stim_accuracy=('stim_accuracy','mean'))
    mouse.to_csv(OUT/'mouse_results.csv',index=False)
    obs=pooled_stat(mouse);p=perm_primary(mouse,obs) if len(mouse)>=3 else np.nan;ci=boot_ci(mouse) if len(mouse)>=3 else [np.nan,np.nan];pred=loo(mouse) if len(mouse)>=3 else {}
    sec={}
    for g in ['ALM','Thalamus']:
        q=mouse[mouse.recording_region==g]
        sec[g]=exact_dir(q) if 2<len(q)<=8 else {'n_mice':int(len(q))}
    early=pooled_stat(mouse,'early_auc') if len(mouse) and mouse.early_auc.notna().all() else np.nan
    counts=mouse.groupby('recording_region').size().to_dict();success=bool(len(mouse)>=12 and all(counts.get(g,0)>=5 for g in ['ALM','Thalamus']) and np.isfinite(obs) and obs>0 and np.isfinite(p) and p<.05 and pred.get('r2_vs_fold_stratum_mean',-np.inf)>0)
    result={'analysis_status':'FROZEN_NUMERICAL_GUO2017_CROSSLOOP_TEST','lock_commit':'561215ed975cff693af2e93795d81dbbeae7aafc','implementation':'parallel exact-eligibility implementation; no scientific definition changed after lock','n_structurally_eligible_assets':len(assets),'n_eligible_sessions':int(len(eligible)),'n_eligible_mice':int(len(mouse)),'mice_by_recording_region':counts,'primary':{'pooled_within_stratum_rank_association':obs,'one_sided_permutation_p':p,'stratified_bootstrap95':ci,'success':success,**pred},'secondary_by_recording_region':sec,'temporal_control':{'early_window_pooled_association':early,'late_minus_early':float(obs-early) if np.isfinite(early) else None},'interpretation':'Frozen numerical source-data test; qualitative source conclusions were already public.'}
    (OUT/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__': main()
