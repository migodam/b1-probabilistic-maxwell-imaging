# 远端原始归档直传draft release：未执行的发布helper

源文件：`src/upload_draft_remote.py`（放到现有Windows主机），`src/launch_remote_upload.py`（本机运行）。本次只检查旧代码、官方API说明和新源语法，**没有运行gh auth token、没有读凭据、没有创建release、没有SSH或上传**。

## 源审查与最小改动范围

R12的`audit/release/upload_ready.py`与`publish_parts.py`已有有限asset清单、实际费用/返回码回执、gh release upload的模式，但依赖本机gh和硬编码历史tag，后者还另生成128MiB传输片。R13本次检索未发现可直接复用的Windows stdlib原始流式上传器。新helper只复用有限清单/逐asset回执思路，保留既有600MB归档，不重新打包、不搬GB到本机，也不跑数值。

GitHub官方API将附件二进制上传到uploads.github.com，使用已存在release ID；成功为201。release metadata含draft属性，GET先核为draft。参见[asset上传说明](https://docs.github.com/en/rest/releases/assets#upload-a-release-asset)、[release接口](https://docs.github.com/en/rest/releases/releases#get-a-release)。本helper固定现有repo `migodam/b1-probabilistic-maxwell-imaging`，使用官方文档当前API版本2026-03-10。

## 主线程实际执行时

先由主线程另行创建并确认**draft** release，取得真实release ID，将远端helper复制至明确指定的Windows路径。所有host/python/script/archive/manifest/receipt路径均由命令行给定，无公开代码硬编码私有SSHkey或host。

远端manifest使用package_evidence实际执行后生成的PACKAGE_RECEIPT.json（parts字段，每项archive/bytes）；也支持显式assets字段name/bytes。archive-dir下每个实际文件的bytes必须与该回执一致；不能拿仅有计划上限的PACKAGE_PLAN冒充实际包。若要上传manifest/restore说明本身，主线程另造同结构、实际stat过的assets清单；不自动造URL或改已有asset。

本机示例参数全为占位，执行前替换真实值：

```text
python src/launch_remote_upload.py \
  --host WINDOWS_SSH_ALIAS \
  --remote-python 'D:\Python\python.exe' \
  --remote-script 'D:\R14\src\upload_draft_remote.py' \
  --archive-dir 'D:\R14\release_shards' \
  --manifest 'D:\R14\release_shards\PACKAGE_RECEIPT.json' \
  --remote-receipt 'D:\R14\private_operational\upload_attempt01.jsonl' \
  --release-id EXISTING_DRAFT_RELEASE_ID \
  --local-receipt '/external/R14/private_operational/upload_launcher01.json'
```

可用`--identity-file`提供用户自行指定的SSHkey路径；代码无固定key路径。PowerShell命令通过UTF16LE EncodedCommand发送，只有非敏感路径与release参数；每个PowerShell字符串按单引号规则转义，token不参加命令拼接。

launcher以参数数组运行`gh auth token --hostname github.com`，只捕获内存，不输出该命令stdout/stderr；再通过SSH的stdin传一行token。token不写文件、不放环境变量、不进进程参数。远端仅从stdin读token到内存，逐1MiB读本地归档经HTTPS发送，明确Content-Length，不将整个asset载入内存。本机只回收小JSON日志，不经过归档字节。local receipt中记录host/远端回执位置，故必须留private_operational，不进Git。

## 可核查行为与边界

远端先GET release且要求draft=true，再分页读取全部已有assets。相同名字不做DELETE/PATCH/clobber：若已有uploaded asset字节数相同，仅记录SKIPPED_EXISTING_SAME_SIZE_NOT_CONTENT_VERIFIED，不声称内容已核验；字节数不同则记录冲突失败。其余资产各尝试一次；失败不自动重试、不删除可能残留starter asset，不影响任何科学运行或原始归档。重试由主线程显式新attempt回执启动。

每asset记录实际源文件bytes、成功send调用字节、HTTP状态、服务器size/asset ID/API返回URL、CPU/wall；失败的send调用可能只发送部分字节，未确认部分明确未知，不伪造精确网络流量。只接受201且名字/size/state与实际文件相符作为UPLOADED_TO_DRAFT；错误body/header/token与原始exception文本不落盘。token在有限进程结束后释放，不声称Python字符串可被严格清零。

这是draft附件上传成功记录，不是匿名公开readback或科学验收。URL从真实API响应取得，draft URL尚不代表公众可访问。发布draft、检查权限、公开匿名读回、更新DOWNLOADS、完整恢复核对均由主线程后续执行。本helper无create/publish/delete接口、无哈希、无新依赖、无服务，不执行任何数值重试。尚未在真实远端/网络环境运行，只有源语法检查，不能称已完成发布验证。
