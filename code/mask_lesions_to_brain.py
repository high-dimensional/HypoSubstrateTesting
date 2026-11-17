#########################################################################
#                                      ###     ###    #######   ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                       #########     #######   #########
#                                                                       #
#                                                                       #
#                      HypoSubstrateTesting                             #
# Testing hypothesized neural substrates with lesions                   #
#                                                                       #
# This script checks whether lesion voxels lie outside the brain and,   #
# if so, multiplies every lesion with the binarised brain mask. Lesions #
# and brain mask are brought onto the same grid by nearest neighbour.   #
# Writes brain_check.csv, brain_mask_report.csv (only if masking was    #
# needed) and brain_mask_status.txt (CLEAN or MASKED) to the log        #
# directory.                                                            #
#                                                                       #
# input:                                                                #
# --lesionDir: binary lesion .nii files                                 #
# --brainMask: brain mask (.nii/.nii.gz)                                #
# --outDir: directory for the masked lesions                            #
# --logDir: directory for the logs                                      #
# --jobs: number of parallel processes                                  #
#                                                                       #
#                                                                       #
# Author: Patrik Bey, patrik.bey@ucl.ac.uk                              #
#                                                                       #
#                                                                       #
# last update: 2026/11/10.                                              #
#                                                                       #
#                                                                       #
#########################################################################


import argparse
import glob
import os
from multiprocessing import Pool

import nibabel as nib
import numpy as np
import pandas as pd

_brain = None


def brain_on_grid(mask_path, ref):
    m = nib.load(mask_path)
    a = np.asarray(m.dataobj) > 0.5
    ii = np.indices(ref.shape).reshape(3, -1)
    w = ref.affine @ np.vstack([ii, np.ones(ii.shape[1])])
    v = np.rint((np.linalg.inv(m.affine) @ w)[:3]).astype(int)
    ok = np.all((v >= 0) & (v < np.array(m.shape)[:, None]), axis=0)
    out = np.zeros(ii.shape[1], bool)
    out[ok] = a[tuple(v[:, ok])]
    assert out.sum() == a.sum(), "brain mask not preserved on lesion grid"
    return out.reshape(ref.shape)


def init(mask_path, ref_path):
    global _brain
    _brain = brain_on_grid(mask_path, nib.load(ref_path))


def check(path):
    data = np.asarray(nib.load(path).dataobj) > 0
    return dict(name=os.path.basename(path), voxels=int(data.sum()), voxels_outside=int((data & ~_brain).sum()))


def mask(args):
    path, out_dir = args
    img = nib.load(path)
    data = np.asarray(img.dataobj) > 0
    out = (data & _brain).astype(np.uint8)
    n_in, n_out = int(data.sum()), int(out.sum())
    assert n_out > 0, f"{path}: lesion empty after brain masking"
    new = nib.Nifti1Image(out, img.affine)
    new.set_sform(img.affine, code=2)
    new.set_qform(img.affine, code=2)
    nib.save(new, os.path.join(out_dir, os.path.basename(path)))
    return dict(name=os.path.basename(path), voxels_in=n_in, voxels_out=n_out, voxels_removed=n_in - n_out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lesionDir", required=True)
    ap.add_argument("--brainMask", required=True)
    ap.add_argument("--outDir", required=True)
    ap.add_argument("--logDir", default=None)
    ap.add_argument("--jobs", type=int, default=32)
    a = ap.parse_args()
    log_dir = a.logDir or os.path.join(os.path.dirname(os.path.abspath(a.outDir)), "logs")
    os.makedirs(log_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.lesionDir, "*.nii")))
    if not files:
        raise SystemExit(f"no .nii lesions in {a.lesionDir}")

    with Pool(a.jobs, initializer=init, initargs=(a.brainMask, files[0])) as pool:
        chk = pd.DataFrame(pool.map(check, files, chunksize=8))
        chk.to_csv(os.path.join(log_dir, "brain_check.csv"), index=False)
        n_aff = int((chk.voxels_outside > 0).sum())
        print(f"CHECK: {n_aff}/{len(chk)} lesions have voxels outside the brain "
              f"({chk.voxels_outside.sum()} of {chk.voxels.sum()} voxels, "
              f"{100 * chk.voxels_outside.sum() / chk.voxels.sum():.3f}%)", flush=True)

        status_file = os.path.join(log_dir, "brain_mask_status.txt")
        if n_aff == 0:
            open(status_file, "w").write("CLEAN\n")
            print("STATUS CLEAN: no lesion voxels outside the brain; no masking performed", flush=True)
            return

        print(f"MASKING: multiplying all {len(files)} lesions with the brain mask -> {a.outDir}", flush=True)
        os.makedirs(a.outDir, exist_ok=True)
        rows = pool.map(mask, [(f, a.outDir) for f in files], chunksize=8)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(log_dir, "brain_mask_report.csv"), index=False)
    open(status_file, "w").write("MASKED\n")
    print(f"STATUS MASKED: removed {df.voxels_removed.sum()} of {df.voxels_in.sum()} voxels "
          f"({100 * df.voxels_removed.sum() / df.voxels_in.sum():.3f}%) from {(df.voxels_removed > 0).sum()} lesions",
          flush=True)


if __name__ == "__main__":
    main()
