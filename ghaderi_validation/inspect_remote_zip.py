import json
from remotezip import RemoteZip

URL='https://zenodo.org/records/17973874/files/2026_Ghaderi_NCOMMS_Zenodo_data_code.zip?download=1'
need=['processed_data/Decoding_NbNeurons5_75.mat','processed_data/Decoding_NbNeurons5_50.mat','processed_data/Optoinhibition_mat.mat','processed_data/psth_10ms.mat','data_helpers/Area_list.mat']
with RemoteZip(URL) as rz:
    infos=rz.infolist()
    out=[]
    for z in infos:
        if any(z.filename.endswith(x) for x in need):
            out.append({'name':z.filename,'size':z.file_size,'compressed_size':z.compress_size})
    print(json.dumps(out,indent=2))
    open('ghaderi_manifest.json','w').write(json.dumps(out,indent=2))
