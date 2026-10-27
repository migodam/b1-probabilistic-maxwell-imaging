from pathlib import Path
import json,re,collections,subprocess
out=Path('research/p3l_r6_task_optimal_maxwell/audit/r6_evidence'); inp=Path('research/p3l_r6_task_optimal_maxwell/inputs/Research_P3LR5T')
a=json.loads((out/'SOURCE_AND_ASSET_AUDIT.json').read_text())
notes={1:'Evidence boundary; B0 chronology is superseded by completed R5B. No upgrade of old A/B failures.',2:'Measure/pushforward, disintegration, source partition; standard conditional mathematics, no new experiment.',3:'Decoder atom and 3D gauge; Full66 existence differs from mass/posterior importance and 3D overlap scope.',4:'Cone theorem, projected nondegeneracy, threshold and eight counterexamples; physical assumptions unverified here.',5:'Resolvent, task adjoint, residual/inverse, missing mass and task bounds; formulas are not cost/certificate success.',6:'Bayes action, representation loss, false teacher example, action-aware distillation and finite-task sufficiency.',7:'Gaussian derivative identity, within/between decomposition, TV/task separation and likelihood-path identity.',8:'Five proposed algorithms; NOT_RUN and mature baselines preserved; plug-in allocation is not exact finite-sample optimum.',9:'Physics features are deterministic functions of observations; only computational/statistical advantage is open.',10:'Prior-art relations; no global novelty claim; singular learning and lifting assumptions distinguished.',11:'Status ledger; conditional theorem versus numerical diagnostic versus unvalidated algorithm.',12:'Prospective stages; R5T NOT_RUN B0 text kept historical and reconciled to R5B.',13:'Paper split and resource judgments belong to parent Codex; no final novelty decision in this audit.'}
lines=['# Complete Source section coverage checklist','', 'All 1,625 Source lines and 455 NEXT_EXPERIMENTS lines were read in this audit. READ means coverage, never proof approval. Every subsection below inherits the stated audit focus.','', '|Source line|Section|Audit focus|','|---:|---|---|']
section=0
for h in a['theory_section_coverage']:
 s=h['heading']
 if s.startswith('# '):
  cn={'一':1,'二':2,'三':3,'四':4,'五':5,'六':6,'七':7,'八':8,'九':9,'十':10,'十一':11,'十二':12,'十三':13}
  for k in sorted(cn,key=len,reverse=True):
   if s.startswith('# '+k+'、'):section=cn[k];break
 text=notes.get(section,'Title/front matter or appendix; source, reproducibility and reference qualification.')
 if '附录' in s:section=0;text='Analytical/decoder receipts and supplied citations; no physical or global novelty conclusion.'
 lines.append(f"|{h['line']}|{s.lstrip('# ')}|READ — {text}|")
(out/'THEORY_SECTION_COVERAGE.md').write_text('\n'.join(lines)+'\n')
# Additional original references exist now, not assumed from old receipt.
ref=json.loads(Path('research/p3l_r5_parallel_ab/audit/EVIDENCE_INPUT_REFERENCES.json').read_text()); refs=[]
for group in ['indispensable_historical_references','new_R5_required_artifacts']:
 for x in ref[group]:
  p=Path(x['path']); refs.append(dict(path=str(p),role=x['role'],exists=p.exists(),bytes=p.stat().st_size if p.is_file() else None))
a['additional_R5_reference_paths_current']=refs
# R5B canonical freeze references and checkpoint manifest paths, no rehash large models.
b=Path('../b1-p3l-r5-b-confirm/research/p3l_r5_b_independent_confirmation')
freeze=json.loads((b/'configs/FREEZE.json').read_text());a['R5B_freeze_references_current']=[dict(path=x['path'],exists=Path(x['path']).exists(),bytes=Path(x['path']).stat().st_size if Path(x['path']).is_file() else None) for x in freeze['files']]
a['missing_theory_input_scope']={'declared_name':'Research_P3LR5_Complete_Research_and_Codex_Handoff.zip','declared_pdf':'Research_P3LR5/Research_P3LR5_ZH.pdf','declared_tex':'Research_P3LR5/source/Research_P3LR5_ZH.tex','status':'NOT_FOUND_IN_R5T_EXTRACTED_INPUT_OR_R5_RESEARCH_SUBTREE','qualification':'R5 SOURCE_INVENTORY lists size/name only; not a substitute for full report. No global filesystem absence claim.'}
(out/'SOURCE_AND_ASSET_AUDIT.json').write_text(json.dumps(a,indent=2))
r=json.loads((inp/'REFERENCES_AND_EVIDENCE.json').read_text());sources=r['sources']
ledger={'search_date':'2026-10-27','skill':'scholarqa-research','transport':'scholarqa-cli','initial_exit_code':1,'retry_exit_code':1,'search_status':'BLOCKED_429','new_verified_citations':[],'claims':[
 {'claim':'Goal-oriented inverse inference / task-weighted posterior quadrature relationship','paper_id':None,'evidence':None,'evidence_tier':None,'relationship':'UNRESOLVED_SEARCH_GAP','confidence':'none; 429'},
 {'claim':'Action-aware distillation prior art','paper_id':None,'evidence':None,'evidence_tier':None,'relationship':'UNRESOLVED_SEARCH_GAP','confidence':'none; 429'},
 {'claim':'Risk-based inverse optimization stopping prior art','paper_id':None,'evidence':None,'evidence_tier':None,'relationship':'UNRESOLVED_SEARCH_GAP','confidence':'none; 429'}],
 'source_supplied_references':[dict(id=s['id'],description=s.get('description'),urls=s.get('urls'),verification_in_this_audit='NOT_VERIFIED; supplied R5T bibliography only; cannot support a new substantive literature conclusion') for s in sources if s['id'].startswith('L')],
 'no_generic_web_substitutes':True,'no_duplicate_MCP_transport':True,'novelty_judgment':'NOT_MADE; no evidence retrieved is not evidence of novelty'}
(out/'CITATION_LEDGER.json').write_text(json.dumps(ledger,indent=2,ensure_ascii=False))
summary={'status':'EVIDENCE_AUDIT_COMPLETE_WITH_LITERATURE_BLOCKED_429','new_Maxwell_states':0,'new_Maxwell_RHS':0,'new_training_updates':0,'GT_arrays_opened':False,'old_history_written':False,'R5B512_status':'REUSED_HISTORICAL_EXPOSED_EVIDENCE_ONLY','old_80pct_gate':'HOLD_PRESERVED','analytical_replay':'PASS_ANALYTICAL_UNIT_NOT_MAXWELL','decoder_replay':'PASS_DECODER_ONLY_PRIOR_DIAGNOSTIC','remaining_P0_units':'NOT_RUN_BY_THIS_AUDITOR','final_novelty_judgment':'PARENT_CODEX_REQUIRED','report':'REPORT_ZH.md'}
(out/'AUDIT_RESULT.json').write_text(json.dumps(summary,indent=2))
print('coverage rows',len(a['theory_section_coverage']),'additional references',len(refs),'missing',sum(not x['exists'] for x in refs))
