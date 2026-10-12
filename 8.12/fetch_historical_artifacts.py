#!/usr/bin/env python3
"""Retrieve specific historical GitHub Actions artifacts (read-only).

Uses the GitHub REST API and repository's GITHUB_TOKEN. Unzips only the two
required CSV reports from each of the 4 independent shock artifacts.
Missing/ambiguous artifacts are a hard failure, not a successful empty audit.
"""
import argparse
import io
import json
import os
from pathlib import Path
import sys
import urllib.request
import urllib.parse
import urllib.error
import zipfile

SHOCKS = ('GPR', 'VIX', 'EPU', 'OIL')
FILENAMES = ('stability.csv', 'irf_stable_only.csv')


def get(url, token, *, binary=False):
    req = urllib.request.Request(url, headers={
        'Authorization': f'Bearer {token}',
        'Accept': 'application/vnd.github+json',
        'X-GitHub-Api-Version': '2022-11-28',
        'User-Agent': 'tvp-gvar-historical-audit',
    })
    # Artifact ZIP endpoint returns an HTTP 302 to a time-limited signed
    # storage URL. Never forward the GitHub bearer token to another host:
    # doing so may produce HTTP 401 and leaks an authorization header.
    class SignedArtifactRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, oldreq, fp, code, msg, headers, newurl):
            nextreq = super().redirect_request(oldreq, fp, code, msg, headers, newurl)
            if nextreq is None:
                return None
            oldhost = urllib.parse.urlsplit(oldreq.full_url).hostname
            newhost = urllib.parse.urlsplit(newurl).hostname
            if oldhost != newhost:
                for key in list(nextreq.headers):
                    if key.lower() in ('authorization', 'x-github-api-version'):
                        del nextreq.headers[key]
                for key in list(nextreq.unredirected_hdrs):
                    if key.lower() in ('authorization', 'x-github-api-version'):
                        del nextreq.unredirected_hdrs[key]
            return nextreq

    opener = urllib.request.build_opener(SignedArtifactRedirect())
    try:
        with opener.open(req, timeout=120) as response:
            data = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f'GitHub artifact download returned HTTP {exc.code}; '
            'verify actions:read permissions and signed URL validity'
        ) from exc
    return data if binary else json.loads(data.decode('utf-8'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default='243e106e-coder/tvp-gvar-replication-private')
    ap.add_argument('--run-id', type=int, default=37872702188)
    ap.add_argument('--output', default='results/audit_artifacts')
    args = ap.parse_args()
    token = os.environ.get('GITHUB_TOKEN') or os.environ.get('GH_TOKEN')
    if not token:
        sys.exit('ERROR: GITHUB_TOKEN is required for private repo artifacts')
    base = f'https://api.github.com/repos/{args.repo}'
    run = get(f'{base}/actions/runs/{args.run_id}', token)
    if run.get('conclusion') != 'success':
        sys.exit('ERROR: historical run was not successful')
    listing = []
    page = 1
    while True:
        data = get(f'{base}/actions/runs/{args.run_id}/artifacts?per_page=100&page={page}', token)
        listing += data.get('artifacts', [])
        if len(data.get('artifacts', [])) < 100:
            break
        page += 1
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {'repository': args.repo, 'run_id': args.run_id,
                'historical_commit': run.get('head_sha'), 'artifacts': []}
    for shock in SHOCKS:
        expected = f'independent-{shock}-{args.run_id}'
        matches = [a for a in listing if a['name'] == expected and not a.get('expired', False)]
        if len(matches) != 1:
            sys.exit(f'ERROR: expected one non-expired artifact {expected}, got {len(matches)}')
        a = matches[0]
        blob = get(f'{base}/actions/artifacts/{a["id"]}/zip', token, binary=True)
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            paths = [n for n in z.namelist() if not n.endswith('/')]
            found = {}
            for filename in FILENAMES:
                candidates = [n for n in paths if Path(n).name == filename]
                if len(candidates) != 1:
                    sys.exit(f'ERROR: {expected} has {len(candidates)} copies of {filename}')
                dest = out / shock / filename
                dest.parent.mkdir(parents=True, exist_ok=True)
                content = z.read(candidates[0])
                if len(content) == 0:
                    sys.exit(f'ERROR: empty {expected}/{filename}')
                dest.write_bytes(content)
                found[filename] = {'zip_member': candidates[0], 'bytes': len(content)}
        manifest['artifacts'].append({'shock': shock, 'artifact_name': expected,
                                      'artifact_id': a['id'], 'files': found})
        print(f'Downloaded validated {shock} reports')
    (out/'source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print('PASS all 8 reports obtained and source manifest recorded')


if __name__ == '__main__':
    main()
