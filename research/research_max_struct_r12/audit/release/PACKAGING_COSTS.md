# R12原始证据准备费用（与科学费用分开）

此汇总没有运行科学正演、反演或评分。按不重叠执行范围列出可测CPU；压缩phase/外层child、运输/解码、历史与resume分别去重。硬件USB断开不是科学失败，原件和中断产物均保留。

|范围|可测CPU秒|wall秒|边界|
|---|---:|---:|---|
|local_primary_ZIP|5.793095|18.35275920899585|Original sixZIP creation; no nested perassetCPU addition|
|local_supplement_ZIP|3.438815|3.78425983292982|124paidcandidate/history/transport files|
|remote_ZIP_archive_phases|64.921875|207.1626125003677|Original perarchive clocks retained; original wholeprocesssummary overwritten by firstresume, wholeprocessCPU remainsnull|
|remote_Zstd|86.40625|45.442448399728164|Existingzstandard6,longwindow27,2threads; all originalZIPs retained|
|local_Zstd|12.407518|20.204833831638098|Original01/02 clocks reused; recovered03-06/supplement1thread; rejectedsupplementwrapper cost still charged|
|local_decoded_byte_comparison|4.277785999999998|6.9246090007945895|Direct decoded byte equality, no hashes|
|copy_from_internal_staging|2.34221|7.485627207439393|External copy plus complete bytecomparison; internal originals preserved|
|existing_wrapper_decode_checks|1.042272|4.155628249980509|Includes interruptedlocal03 frame rejection, not a scientificfailure|
|network_and_remotecontrol_localprocesses|33.932323000000004|None|Includes abortedtransferattempts/resumes; remotephysicalCPU separatelyabove; concurrentwallnot summed|
|internal_completed_transport_before_userpause|18.442564|None|Prior completedmetadata/twoassettransfers only; interruptedremote03 terminalCPUunknown|
|remoteZIP_decode_on_external|3.5200280000000013|13.292924832552671|Sourcecopy/remote01 internaldecode not added again|
|prior_internal_remote01_decode|0.673812|1.457419749815017|Performed once before userpausedinternalwrites; source movedbycopy only|

已知不重叠CPU下界237.198548秒；精确全程CPU=`null`，全程wall=`null`，货币费用=`null`。并行运输wall不能相加，缺失不当0。

最初remoteZIP组合回执在首次resume时被短resume时钟覆盖，但每包原CPU/wall完整保留，故使用逐包phase和而不伪造wholeprocess值。原local03截断Zstd保留并重新压缩，解码对照原ZIP字节，没有hash。

ROOT_STORAGE_INCIDENT与EXTERNAL_DISCONNECT_INTERRUPTION记录硬件/操作中断；全部source科学数组、ZIP原件、内盘旧副本和partial保留。root实际split/upload费用另见PUBLICATION_PARTS_MANIFEST.json及发布回执，不在本准备下界里重复计。Git最新文稿优先于原始ZIP中的历史文稿快照。
