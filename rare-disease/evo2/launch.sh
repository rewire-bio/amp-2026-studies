#!/usr/bin/env bash
# Score the rare-disease windows on one Runpod H100 with the pinned evo2-cloud-operations image.
# The job runs detached on the pod; this script polls for DONE, copies results back and deletes the pod.
# A hard cap (MAX_HOURS, default 5, about USD 20 at USD 3.99/h) also stops the pod from the inside.
#
#   RUNPOD_API_KEY=... bash rare-disease/evo2/launch.sh          # from the repo root, after make_windows.py
#   DATA=rare-disease/evo2/data3 OUT=rare-disease/evo2/run3 bash rare-disease/evo2/launch.sh   # another window set
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OPS=${OPS:-$HERE/../../../evo2-cloud-operations}
source "$OPS/cloud/common.sh"

PLATFORM_ID=RD-H100
IMAGE_REF=${IMAGE_REF:-ghcr.io/rewire-bio/evo2-cloud-operations@sha256:22970720ff8bc05a00460a6066a70ed71951001aa52965141591a0ca6dcc968b}
GPU=${GPU:-NVIDIA H100 80GB HBM3}
MAX_HOURS=${MAX_HOURS:-5}
DATA=${DATA:-$HERE/data}
OUT=${OUT:-$HERE/run}
API=https://api.runpod.io/v2
mkdir -p "$OUT"
KEY="$OUT/.ssh-key"
BUNDLE=$(mktemp -d)

rp() {
  local body status
  body=$(curl -sS -w '\n%{http_code}' -H "Authorization: Bearer $RUNPOD_API_KEY" -H "Content-Type: application/json" "$@")
  status=${body##*$'\n'}; body=${body%$'\n'*}
  if [[ $status != 2* ]]; then log "Runpod API $status: $body"; return 1; fi
  printf '%s' "$body"
}

POD_ID=""
cleanup() {
  set +e
  if [[ -n $POD_ID ]]; then log "Deleting pod $POD_ID"; rp -X DELETE "$API/pods/$POD_ID" >/dev/null; mark terminated; fi
  rm -rf "$KEY" "$KEY.pub" "$BUNDLE"
}
trap cleanup EXIT

# Code bundle: the evo2-cloud-operations model loader (src, followup) plus this directory's scorer.
git -C "$OPS" archive --format=tar -o "$BUNDLE/code.tar" HEAD src followup
tar -rf "$BUNDLE/code.tar" -C "$HERE" score_windows.py
gzip "$BUNDLE/code.tar"
shasum -a 256 "$BUNDLE/code.tar.gz" | cut -d' ' -f1 > "$OUT/code_sha256.txt"
git -C "$OPS" rev-parse HEAD > "$OUT/ops_commit.txt"

new_ssh_key
mark launch_requested
BODY=$(jq -n --arg name "evo2-rd-$(date -u +%Y%m%d%H%M%S)" --arg image "$IMAGE_REF" --arg gpu "$GPU" \
  --arg key "$(cat "$KEY.pub")" --arg max "$((MAX_HOURS * 3600))" '{
    name: $name, image: $image, gpu: {id: $gpu, count: 1}, cloud: "SECURE", disk: 150,
    ports: ["22/tcp"], env: {PUBLIC_KEY: $key, MAX_SECONDS: $max}, cmd: ["/opt/evo2ops/pod-start.sh"]}')
POD=$(rp -X POST "$API/pods" -d "$BODY") || { log "No pod could be created"; exit 3; }
POD_ID=$(echo "$POD" | jq -r .id)
log "Pod $POD_ID created ($GPU)"
mark pod_created
for _ in $(seq 1 120); do
  POD=$(rp "$API/pods/$POD_ID" || true)
  SSH_HOST=$(echo "$POD" | jq -r '.ssh.direct.host // empty' 2>/dev/null)
  SSH_PORT=$(echo "$POD" | jq -r '.ssh.direct.port // empty' 2>/dev/null)
  [[ -n $SSH_HOST && -n $SSH_PORT ]] && break
  sleep 15
done
[[ -n ${SSH_HOST:-} ]] || { log "Pod never exposed SSH"; exit 3; }
mark running
echo "$POD" | jq '{provider: "runpod", gpu: .gpu, region: .dataCenterId, instance_id: .id, image: .image,
  cost_per_hr: .cost, host_cuda: .cudaVersion}' > "$OUT/instance.json"
SSH_USER=root
log "Pod at $SSH_HOST:$SSH_PORT, $(jq -r .cost_per_hr "$OUT/instance.json") USD/h"

wait_for_ssh
ssh_cmd "mkdir -p /work/data /work/code /work/out"
copy_to "$BUNDLE/code.tar.gz" /work/
copy_to "$DATA/refs.tsv.gz" "$DATA/vars.tsv.gz" /work/data/
ssh_cmd "tar -xzf /work/code.tar.gz -C /work/code"
mark workload_start
ssh_cmd "cd /work && setsid nohup bash -c 'nvidia-smi --query-gpu=name,driver_version --format=csv,noheader > /work/out/host.txt; \
  python /work/code/score_windows.py > /work/out/score.log 2>&1; touch /work/out/DONE' > /dev/null 2>&1 < /dev/null & disown" || true

deadline=$((SECONDS + MAX_HOURS * 3600 - 600))
until ssh_cmd "test -f /work/out/DONE" 2>/dev/null; do
  [[ $SECONDS -lt $deadline ]] || { log "Time limit reached before DONE"; break; }
  log "$(ssh_cmd 'tail -1 /work/out/score.log' 2>/dev/null || echo waiting)"
  copy_from "/work/out/." "$OUT/" 2>/dev/null || true
  sleep 300
done
mark workload_end
copy_from "/work/out/." "$OUT/" || log "Copying outputs back failed"
jq -r .status "$OUT/receipt.json" 2>/dev/null || true
