"""Pack only public-score predictions for independent sealing; no truth API."""
from pathlib import Path
import argparse,json

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--results',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    public=json.loads(Path(a.manifest).read_text());source=Path(a.results);out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    counts={'OK':0,'FAIL':0,'NOT_RUN':0}
    for record in public['records']:
        obs=record['observation_id'];methods={}
        for method,K in public['planned_methods_and_starts'].items():
            starts=[]
            for i in range(K):
                task_id=obs+'__'+method.replace(':','_')+'__s'+str(i)
                path=source/(task_id+'.json')
                if not path.exists():counts['NOT_RUN']+=1;continue
                raw=json.loads(path.read_text(encoding='utf-8-sig'));status=raw['status'];counts[status]=counts.get(status,0)+1
                result=raw.get('result') or {'counts':raw.get('partial_costs')}
                fields=('m','ranking_objective','relative_residual','counts','task_wall_s','wrapper_status','wrapper_stop_reason','stage','stages_reached')
                starts.append({'start_id':'s'+str(i),'status':status,'task':raw['task'],'wall_s':raw['wall_s'],
                               'result':{k:result[k] for k in fields if k in result}})
            methods[method]={'starts':starts}
        (out/(obs+'.json')).write_text(json.dumps({'observation_id':obs,'methods':methods,'truth_access':'NONE'},indent=2))
    print(json.dumps({'packed_observations':len(public['records']),'prediction_dispositions':counts,'truth_access':'NONE'}))
if __name__=='__main__':main()
