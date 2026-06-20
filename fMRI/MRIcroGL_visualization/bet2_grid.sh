    #!/usr/bin/env bash
    set -euo pipefail

    # === 输入 ===
    IN="${1:?Usage: bet2_grid.sh <input_T1.nii> [OUTDIR]}"
    OUTDIR="${2:-bet2_grid}"
    #
    Fs=(0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55)
    Gs=(-0.2 -0.1 0 0.1 0.2)

    mkdir -p "$OUTDIR"
    echo "[1/4] Reorient to std"
    fslreorient2std "$IN" "$OUTDIR/reorient.nii.gz"

    echo "[2/4] Crop FOV (remove neck)"
    robustfov -i "$OUTDIR/reorient.nii.gz" -r "$OUTDIR/crop.nii.gz"

    echo "[3/4] Bias-field correction (FAST)"
    fast -B -o "$OUTDIR/crop_bias" "$OUTDIR/crop.nii.gz"
    REF="$OUTDIR/crop_bias_restore.nii.gz"

    echo "[4/4] BET2 grid search"
    for f in "${Fs[@]}"; do
    for g in "${Gs[@]}"; do
    tag="f${f//./}g${g//./}"
    out="$OUTDIR/brain_${tag}"
    bet2 "$REF" "$out" -f "$f" -g "$g" -m
    #
    # fslmaths "${out}_mask" -fillh -dilM -ero "${out}_mask"
    #
    slicer "$REF" "${out}_mask.nii.gz" -S 3 120 "$OUTDIR/quick_${tag}.png"
    echo " -> $tag done"
    done
    done

    # 
    pushd "$OUTDIR" >/dev/null
    slicesdir -o "$(basename "$REF")" brain_f*_mask.nii.gz
    popd >/dev/null

    echo "done:$OUTDIR/slicesdir/index.html"


