"""Separate measurement-only warm-start Full-GN supplement; primary B0 is read-only."""
from pathlib import Path
import argparse,json,os,sys,subprocess,csv
sys.dont_write_bytecode=True
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[key]='1'
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from gsda.fem_model import FEMModel
from gsda.fwi import FWIOptions,complete_outputs,dump
from gsda.gn import refine,GNOptions
from gsda.accounting import BudgetBook
from gsda.guard import OnlineGuard


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--config',type=Path);ap.add_argument('--measurement',type=Path)
    ap.add_argument('--baseline',type=Path);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--options',type=Path,required=True);ap.add_argument('--gn-options',type=Path)
    ap.add_argument('--job',default='FULL_GN_SUPPLEMENT')
    ap.add_argument('--cases',nargs='+');ap.add_argument('--online-root',type=Path);ap.add_argument('--b0-root',type=Path)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if args.cases:
        if args.online_root is None or args.b0_root is None:ap.error('batch requires --online-root and --b0-root')
        records=[]
        for case in args.cases:
            book=BudgetBook(ROOT,job=f'{args.job}_{case}');record=dict(case=case,status='NOT_RUN')
            try:
                book.check(new_job=True);folder=args.output/case;folder.mkdir(parents=True,exist_ok=True)
                command=[sys.executable,str(Path(__file__).resolve()),'--config',str(args.online_root/case/'config.json'),'--measurement',str(args.online_root/case/'observation.npz'),'--baseline',str(args.b0_root/case),'--output',str(folder),'--options',str(args.options),'--job',f'{args.job}_{case}']
                if args.gn_options:command+=['--gn-options',str(args.gn_options)]
                with (folder/'stdout.log').open('w') as log:process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
                record.update(returncode=process.returncode,status='PROCESS_COMPLETED' if process.returncode==0 else 'FAILED')
                if (folder/'ONLINE_SEALED.json').exists():record['online_status']=json.loads((folder/'ONLINE_SEALED.json').read_text())['status']
            except Exception as error:record['reason']=str(error)
            records.append(record);dump(args.output/'BATCH_STATUS.json',dict(cases=args.cases,attempts=records,remaining=args.cases[len(records):],cohort='DATA_EXPOSED_MEASUREMENT_ONLY_REFINEMENT'))
            print(json.dumps(record),flush=True)
        return
    if any(x is None for x in [args.config,args.measurement,args.baseline]):ap.error('single case requires --config --measurement --baseline')
    guard=OnlineGuard(ROOT).install();book=BudgetBook(ROOT,job=args.job);status='FAILED';model=None;result=None
    try:
        book.check(new_job=True);config=json.loads(args.config.read_text());config['assembly_strategy']='cached_combined'
        options=json.loads(args.options.read_text());prior=FWIOptions(**options.get('fwi_options',options))
        gn_source=json.loads(args.gn_options.read_text()) if args.gn_options else {}
        gn_options=GNOptions(**gn_source.get('gn_options',gn_source))
        with np.load(args.measurement,allow_pickle=False) as observation:
            if set(observation.files)!={'y'}:raise ValueError('online NPZ must contain only y')
            y=observation['y'].copy()
        model=FEMModel(config,book,phase='R4_FULL_GN_SUPPLEMENT')
        result=refine(model,y,args.baseline,args.output,prior,gn_options);status=result['status']
    except BaseException as error:
        book.failure('R4_FULL_GN_SUPPLEMENT','refine_gn',error);dump(args.output/'FAILURE.json',dict(type=type(error).__name__,message=str(error)));raise
    finally:
        artifact_error=None
        if model is not None and (args.output/'ONLINE_SEALED.json').exists():
            try:
                complete_outputs(model,args.output)
                fields=['time_unix','kind','outer','alpha','objective','data_objective','prior','scaled_projected_gradient_inf','damping','accepted','actual_reduction','predicted_reduction','ratio','optimizer_success','optimizer_stop','iterations']
                with (args.output/'optimization_trace.csv').open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader()
                    for line in (args.output/'TRACE.jsonl').read_text().splitlines():
                        row=json.loads(line);writer.writerow({key:row.get(key) for key in fields})
            except Exception as error:artifact_error=str(error);book.failure('R4_FULL_GN_SUPPLEMENT','output aliases',error)
        b0cost={}
        if args.baseline is not None and (args.baseline/'RECEIPT.json').exists():
            old=json.loads((args.baseline/'RECEIPT.json').read_text());b0cost={key:old[key] for key in ['full_fe','rhs','wrapper_process_cpu_seconds','wall_seconds']}
        receipt=book.receipt(status=status,output=str(args.output),artifact_error=artifact_error,b0_cost=b0cost,not_cold_start=True,cohort='DATA_EXPOSED_MEASUREMENT_ONLY_REFINEMENT')
        if b0cost:
            receipt['cumulative_full_fe']=receipt['full_fe']+b0cost['full_fe'];receipt['cumulative_rhs']=receipt['rhs']+b0cost['rhs']
            receipt['cumulative_wrapper_process_cpu_seconds']=receipt['wrapper_process_cpu_seconds']+b0cost['wrapper_process_cpu_seconds']
            receipt['cumulative_wall_seconds']=receipt['wall_seconds']+b0cost['wall_seconds']
        dump(args.output/'RECEIPT.json',receipt);dump(args.output/'runtime_cost.json',receipt);guard.save(args.output/'FILE_ACCESS_SUMMARY.json')
    print(json.dumps(dict(status=status,output=str(args.output),cost=book.snapshot(pid=os.getpid()))),flush=True)

if __name__=='__main__':main()
