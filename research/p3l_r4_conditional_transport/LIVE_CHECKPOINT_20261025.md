# R4 live checkpoint — IN PROGRESS, do not finalize yet

Worktree b1-p3l-r4-conditional-transport, branch research/p3l-r4-conditional-transport-20261025. Latest committed code95ac480e, plus fresh physical queue/new audit docs uncommitted. User R4 task is active, originalP1 fixed; no historical ledgers touched, no SHA256, no cloud/3D/Mie. Native6.1sol agents, medium/high, noDeepSeek/ultra inheritance.

## Already closed new work

- Full66 rawSpline/diag/Gaussian/triangle/CFM actual3seed32k training and curves;15anchor models1024/256;45riskactions inclphase8/16 revisions,3momentproviders;6boundedgeometryCFMfits.
- Actual mixed forward/reverseKL refinement:3family pairedseed1,3000updates each,203264+11264RHS/pair. All3 replication gates FAIL properNLLcondition; no seed2/3. Group results/physics_transport_paired_v2/CONTROLLER_RESULT.json. V1MPSdouble reduction0RHSfailure retained; V2 onlyCPUdouble reduction fix/prefix, no new budget.
- FirstrealVJPpreflight64RHS originalFAIL compareddifferentFloat64/Float32 inputs. Saved-query sameinputauditPASS, error2.4e-8,0newRHS, no threshold relaxation.
- Postfit9modelDEV generation+score+paid45504RHS andcorrectedMH49536RHS complete: results/POSTFIT_DEV_V2_EXECUTION_COMPLETE.json. Wrong parentrelease key guard failed beforephysics, originalrelease/logretained, explicit BINDING_FIXED metadata repair; no generation replay. Maintoolfinish_r4_postfit_release_binding.py.
- Same-cacheFull/Richardson/GMRES/fixedaccuracyDEV4arms105010RHS, allpaths/decisions/proofsagree. Boundary3beta×3seed×4parents×256 all4arms113687RHS; no misses, largermovesmostlyreject, asymmetricCFMCPUoverlapdisclosed. Audits under AUDIT/global_gmres and interval_boundary.
- QUIET FRESH LOCKED8 fourarms COMPLETE: results/IV/interval_chain_locked_v1/COMPARISON.json; Full65536RHS71.397s,Richardson4608RHS24.153s(66.17%wall saving),GMRES14691RHS37.605s,fixedaccuracy20759RHS46.832s;0misses,105594totalRHS. Float64validated,notdirectedintervalcertificate/notposteriormixing.
- GeometryEulerdefect→Heunrevision retained; strict2%distance auditHeun9/288,Euler4/288;directallmatched. NounmatchedGTriskranking. Archive strictRMS>=.02+noiseball+DeltaPhi<=2+priorDelta<=2+progress:worldcounts INCLUDINGanchor DEV0..7=15,12,7,16,28,16,15,10.119entries over8observations,111newcandidates,candidate_only/nullweights; no posterior/modeclaim. Sourcepoolincludesordinary/strongcontrols, notbracketonly. results/III/geometric_candidate_archives_v2.
- DINO6fits accuracygateALLFAIL; ranker3fits worseconstantfield-update. Local16/8/8parentcorpusnotfullLazyDINOreplication. Actualfacts AUDIT/PHYSICS_MECHANISM_FACTS.json.
- BoundedSPDgeometryCFM matched baseline initialization/normalization/losscoordinates: Jr=Jx*T,Wr=I+3Vdiag(s²/(1+s²))Vt,randomsameeigenspectrum. All6fits selected3000. DEV128NF32seedlossmean geometryenergy.132520718/MSE.037114368 vsnone.158894336/.053325340,random.164919606/.057273991. Still MUCHworseCNN; noFullchecks yet. ROI90empcoverage.9784,rankcoverage.9551 overnominal.9. results/geometry_weighted_cfm_dev_v1.

## Fresh768 campaign ACTIVE — NO heldoutGT has been scored

- FROZEN_LOCKED.json +FROZEN_CONFIRMATION.json and FROZEN_DENSITY_BOTH_COHORTS.json created BEFORE anynewobs. Riskcomponent FINALv3 path AUDIT/risk/fresh_freeze_v3/RISK_COMPONENT_FROZEN.json,45riskmethods.26densitymodels,allCFMfixedNFE32,K128. Bothcohorts identicalfrozenchoices, no postL1adapt.
- PrimaryriskSplinefusionselectedbyR2VAL256mean.0071267302458 vsDiag.0071301649635; DEVDiagbetterbutNOTswitching. ControlsstrongCNN,Diagfusion,permuted,decoder/random,oldreoptimizedR2class. NewlocalCNNsamefunctionclassR2 (headrescale), don'tcallnovelarchitecture.
- Closed6newfusionBEST.pt have moment_K=6418 headerbug; weightsunchanged. fusion_provider_compat.py readonlyrecoversK32/0 fromsealedproviders/TRAIN+VALgates. Closedbytespreserved. CanonicalfreshCPUbatch64savedarraysdefinescores; old2e-6replayFAILpreserved, v2 1e-5technicalheuristiconly.
- Main physical queue session79894: tools/queue_r4_fresh_physics.py; log logs/fresh_physics_queue.log; steps/results under results/fresh_physics_queue_v1. Generates all512+256truthobservations (done), then64-parent native96+probit64contexts chunks, quietkernelafterfirst64(done), exportsprovenanchorconditionsafterallgroups. Perchunk~86s/RSS~260MB. Envstoredlosslesscodec. Watch currentprogress, doNOTreplayexistingchunks.
- Risk background watcher session20557 waits global results/FRESH_PHYSICS_AND_CONTEXTS_COMPLETE.json60s cadence; then all45methods features/infer/export BOTHgroups via BookIV, NO SCORE. log logs/fresh_risk_online_wait_v1.log; resultmarker results/fresh/risk_online_queue_v1/{COMPLETE,ERROR}.json. Outputs results/fresh/{locked,confirmation}/risk45.
- Density watcher session34881 waits SAME globalmarker; then26models infer CPU16 using tools/evaluate_r4_density_frozen_nfe.py, bothgroups,K128, BookIV7200s/2GiB, NO SCORE. log logs/fresh_density_infer_both.log; marker results/FRESH_DENSITY_ONLINE_BOTH_COMPLETE.json; outputs results/fresh/{group}/density_online/ONLINE_COMPLETE.json.
- Need ALL risk+density BOTHcohorts onlineclosed BEFORE firstheldoutGTscore. Mainfrozenmodel/risk/sourceidentities immutable; only addnewtools/reports.

## Next needed actual work

1. Finishphysicalqueue/watchers. Checkfailures/partialtruth-free; doNOTreplaceparents/resetbudget.
2. Mainrelease FullPPC4groups (rawSpline_s1explicit;anchorCFMnone/geometry/random_s1NF32densitynone) fixedsample0all768 +first16L1first64. Inertrequest configs/INERT_FRESH_DENSITY_FULL_PPC_REQUEST.json.56832RHS StageIV. No candidatefilter. Existing run_r4_probability_audit freeze/run, separately offlineSBCscore afterallonlineclosed.
3. Rawriskdiagnostics all768×6actions (R2/strong_s1/Spline_s1 raw/clipped),24576+12288RHS StageIV. risk45/exports/RAW_FORWARD_FROZEN.json +tools/r4_raw_material_forward.py run; explicitstageIV release/CLI matching.
4. Representative first64L1 primarySpline_s1+R2raw/clipped +ONE sharedordinaryFull-Phi perparent;16states/256RHS eacharm,total81920hardRHS StageIV. Riskwatcher willprepare request afteronline, no physical. tools/r4_fresh_representatives.py --stageIV matchingrequest/release/actualreceipt; originalPhi282.873792 unchanged. SeparateGTscoreafterclosed. Notposterior means.
5. Offline score ALL risk45/density26 bothgroups separatePIDs/Books; usecanonicalarrays, allfailures. CI98.3333 for3primary tests; othercomparisons exploratory notglobalFWER. Frozenfamilyrisk3seedlossmean, notensembleimages. Couldthenaggregate768descriptively; no retuningafterL1.
6. AdaptiveCFM numericaldensity audit NEW agent/r4_density authorized0RHS max900s: DEV0/1 only3modelsS1rawmatched/anchor-none/geometry, CPUdoublefull66divDOP853 twotols1e-5/1e-4 and1e-6/1e-5. NeverchangesfrozenNF32/lockedmethod. ExistingfixedNFEtracefailed3/4; auditdistinguishesintegrationerror fromdistributionquality. Agentwillreport.
7. FinalrequiredChineseR4reports/figures/alloriginalfailures, resourceclosure, independentstatsaudit. GO-A/B/C/D/E separate; don'tletq/samplerfailsnegatevalidriskbenefit orviceversa.
8. Publish existingB1Git/Pages +evidenceRelease, anonymousrawMarkdown/readback. R4 NOTYETPUBLIC. Preserve~100GB+physicalarrays locally, packagewithstreaming/splitcompressedassets andtruthfulmanifest; noSHA256. Existingtools/package_p3r2_evidence.py usefulpattern, nofollowoldsymlinks/noOOM.
9. FinalconciseChinese links andmainmetrics, no prematureCOMPLETE. Memoryused MEMORY.md698-700, finalappend exactlyoneoai-mem-citationblock.

## Scientific boundaries to retain

Full66properNLLgain !=correctposterior; fixedsample0Gamma256PPCnecessary, notsufficient. CFMnumericdensityunqualifiedcan'texactMH. RiskmeanimagecanoutsideP1decoder; physicalrepresentative separate. Observednoiseballnotpairdistance<=1/exactfiber or95%credible. Archive111candidatesarenotposterior mass. StaticMaxwellweightsmayhelpCFMbutneitherdynamicbracketnorA17incrementestablished. PlainGaussianproviderasgoodFlowforrisk; directCNNretraininggainnotnewarchitecture. ClosestpriorartverifiedLazyDINO2026,amortizedVI2024,preconditionedNF2021,retrospectiveacceptancebounds2014; sourceURLs inAUDITdocs. No blanketTAPnoveltyclaim.
