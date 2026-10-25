"""Anonymous, bounded delivery readback; no scientific execution or hashes."""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / 'research/p3l_r4_conditional_transport'
REPO = 'migodam/b1-probabilistic-maxwell-imaging'
TAG = 'research-p3lr4-20261025'
TARGET = 'c37650dbf5bfac537f5789d7d1bbcb2b875f88c3'
WEB = 'https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r4_conditional_transport/'
API = 'https://api.github.com/repos/' + REPO
RELEASE = 'https://github.com/' + REPO + '/releases/tag/' + TAG


def get(url, limit=1024*1024, extra_headers=None):
    headers = {'User-Agent': 'R4-anonymous-delivery-readback', 'Cache-Control': 'no-cache'}
    headers.update(extra_headers or {})
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=30) as response:
        data = response.read(limit)
        return data, {'url': url, 'status': response.status,
                      'bytes_read': len(data), 'content_length': response.headers.get('Content-Length'),
                      'content_range': response.headers.get('Content-Range'),
                      'redirect_host_only': urllib.parse.urlsplit(response.geturl()).hostname}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='NEW delivery-audit directory')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    rows = []
    try:
        release = json.loads(get(API + '/releases/tags/' + TAG)[0])
        assert release['draft'] is False and release['html_url'] == RELEASE
        assert release['target_commitish'] == TARGET
        remote = {}
        for page in (1, 2, 3):
            batch = json.loads(get(API + f'/releases/{release["id"]}/assets?per_page=100&page={page}')[0])
            for asset in batch:
                assert asset['name'] not in remote
                remote[asset['name']] = {key: asset[key] for key in
                                        ('id', 'name', 'size', 'state', 'browser_download_url')}
            if len(batch) < 100:
                break
        packdir = S / 'packaging/full_evidence_20261025_v2'
        pack = json.loads((packdir / 'RECEIPT.json').read_text())
        expected = {v['file']: v['bytes'] for v in pack['parts']}
        small = json.loads((S / 'packaging/final_delivery_assets/SMALL_UPLOAD_RECEIPT.json').read_text())
        for asset in small['assets']:
            expected[asset['name']] = asset['bytes']
        extra = S / 'packaging/final_delivery_assets/R4_SMALL_ASSETS_UPLOAD_RECEIPT.json'
        expected[extra.name] = extra.stat().st_size
        for name, size in expected.items():
            assert remote[name]['size'] == size and remote[name]['state'] == 'uploaded', name
        assert all(v['state'] == 'uploaded' for v in remote.values())
        checks = [(WEB + 'report.html', '白化残差能量'),
                  (WEB + 'gallery.html', 'P3L-R4'),
                  ('https://raw.githubusercontent.com/' + REPO +
                   '/main/research/p3l_r4_conditional_transport/RESEARCH_REPORT_ZH.md', '13.10%'),
                  (RELEASE, 'P3L-R4')]
        for url, needle in checks:
            data, row = get(url)
            row['expected_text_present'] = needle in data.decode(errors='replace')
            assert row['status'] == 200 and row['expected_text_present'], url
            rows.append(row)
        # Tiny prefix checks, not full downloads or content-integrity certificates.
        locals_ = {
            'compact_core.tar.gz.part0001': packdir / 'compact_core.tar.gz.part0001',
            'raw_physics_and_arrays.tar.part0113': packdir / 'raw_physics_and_arrays.tar.part0113',
            'Research_P3LR4_Complete_Results_and_Handoff.zip': S / 'packaging/pro_handoff_publication/Research_P3LR4_Complete_Results_and_Handoff.zip',
            'Research_P3LR4_Postpack_Supplement.zip': S / 'packaging/postpack_supplement_20261025_v1/Research_P3LR4_Postpack_Supplement.zip',
        }
        for name, local in locals_.items():
            canonical = remote[name]['browser_download_url']
            attempts = []
            for revision in (False, True):
                url = canonical + ('?revision=' + TARGET if revision else '')
                try:
                    data, row = get(url, 128, {'Range': 'bytes=0-127'})
                    row['local_128_byte_prefix_matches'] = data == local.open('rb').read(128)
                    assert row['status'] in (200, 206) and row['local_128_byte_prefix_matches'], name
                    row['asset'] = name
                    row['attempts'] = attempts
                    rows.append(row)
                    break
                except urllib.error.HTTPError as exc:
                    attempts.append({'canonical_url': canonical, 'cache_bust': revision, 'status': exc.code})
                    if revision:
                        raise
        audit = {'status': 'ANONYMOUS_NAME_SIZE_STATE_VERIFIED', 'remote_assets': remote,
                 'expected_before_receipt_assets': expected, 'asset_count': len(remote),
                 'part_count': len(pack['parts']), 'part_bytes': sum(v['bytes'] for v in pack['parts']),
                 'authentication_used': False, 'hashes_used': False}
        (args.out / 'REMOTE_ASSET_AUDIT.json').write_text(json.dumps(audit, indent=2) + '\n')
        result = {'status': 'PUBLICATION_ANONYMOUS_READBACK_PASSED', 'release_url': RELEASE,
                  'release_id': release['id'], 'tag': TAG, 'target_commit': TARGET,
                  'published_at': release['published_at'], 'report_url': WEB,
                  'gallery_url': WEB + 'gallery.html',
                  'anonymous_markdown_url': checks[2][0],
                  'reading_package_url': remote['Research_P3LR4_Complete_Results_and_Handoff.zip']['browser_download_url'],
                  'supplement_url': remote['Research_P3LR4_Postpack_Supplement.zip']['browser_download_url'],
                  'full_part_count': len(pack['parts']), 'full_part_bytes': audit['part_bytes'],
                  'asset_count_before_receipt_attachments': len(remote),
                  'scientific_qualification_changed': False, 'research_ledger_sequence': 1522869,
                  'archive_source_commit': '60066394392fe04fb6a5b30371ffa34c978ade9d',
                  'RHS': 0, 'training': 0, 'inference': 0, 'authentication_used': False,
                  'hashes_used': False, 'checks': rows, 'wall_seconds': time.monotonic()-start}
    except BaseException as exc:
        result = {'status': 'PUBLICATION_READBACK_INCOMPLETE', 'error': repr(exc),
                  'checks': rows, 'wall_seconds': time.monotonic()-start, 'RHS': 0}
        (args.out / 'RECEIPT.json').write_text(json.dumps(result, indent=2)+'\n')
        raise
    shutil.copy2(Path(__file__), args.out / 'as_run_code.py')
    (args.out / 'RECEIPT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'release_url', 'full_part_count', 'full_part_bytes', 'asset_count_before_receipt_attachments', 'wall_seconds')}))


if __name__ == '__main__':
    main()
