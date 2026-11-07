"""Post-run stdlib distribution metadata, zero numerical imports/FW/GPU context."""
import argparse,json,sys,time,subprocess
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();plan=json.loads(a.plan.read_text());command=plan['tasks'][0]['command'];python=command[0];vendor=command[command.index('--vendor-path')+1];c=time.process_time();wall=time.perf_counter()
 code="import sys,json,importlib.metadata as m; sys.path.insert(0,"+repr(vendor)+"); packages={name:dict(version=m.distribution(name).version,distribution_root=str(m.distribution(name).locate_file('')),metadata_directory=str(m.distribution(name)._path)) for name in ('numpy','scipy')}; print(json.dumps(dict(sys_executable=sys.executable,sys_version=sys.version,sys_path=sys.path,packages=packages)))"
 process=subprocess.Popen([python,'-c',code],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);stdout,stderr=process.communicate()
 from twod_queue import external_cpu
 external=external_cpu(process)
 result=dict(status='POST_RUN_INVENTORY_COMPLETE' if not process.returncode else 'POST_RUN_INVENTORY_FAILED',created_utc=datetime.now(timezone.utc).isoformat(),scope='Post-run inventory at recorded time; not an original startup receipt or proof of historical unchanged package state',production_plan=str(a.plan),production_launcher_command=command,metadata_child_command=[python,'-c',code],vendor_path_from_real_plan=vendor,metadata=json.loads(stdout) if not process.returncode else None,stderr=stderr,returncode=process.returncode,external_child_cost=external,controller_CPU_seconds_before_write=time.process_time()-c,wall_seconds=time.perf_counter()-wall,forbidden_imports_not_requested=['numpy','scipy','torch','cuda'],numerical_probe=0,new_FW=0,GPU_context_created=False,GT_access='NONE')
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return process.returncode
if __name__=='__main__':raise SystemExit(main())
