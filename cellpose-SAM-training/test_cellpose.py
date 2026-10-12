"""Smoke test for the Cellpose-SAM SMG basement-membrane / epithelium models.

Runs a model in 3D on the small test volume in test_data/ (a crop of a volume held out
from SMG-BM-15V training) and scores it against the annotation, using the same
post-processing and metrics as the notebooks. With --train it first fine-tunes a
throw-away model for a few epochs on 2D slices of the test volume, to check that
training works.

Usage (from this folder):
    python test_cellpose.py --model models/SMG-BM-v0907   # trained weights from Figshare
    python test_cellpose.py --train                       # quick training check
"""
import argparse
import os

import numpy as np
from cellpose import io, models, transforms
from scipy.ndimage import binary_fill_holes, label

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "test_data")


def keep_largest_object(binary_image):
    labeled_array, num_features = label(binary_image)
    if num_features == 0:
        return np.zeros_like(binary_image)
    sizes = np.bincount(labeled_array.ravel())
    sizes[0] = 0
    return (labeled_array == sizes.argmax()).astype(binary_image.dtype)


def segment(model, img):
    _, flows, _ = model.eval(img, normalize=False, do_3D=True, z_axis=0, channel_axis=-1,
                             cellprob_threshold=0.0, batch_size=64)
    mask = flows[2] > 0
    mask = binary_fill_holes(mask, structure=np.ones((3, 3, 3), "int"))
    return keep_largest_object(mask).astype("uint8")


def scores(pred, gt):
    pred, gt = pred > 0, gt > 0
    tp = np.sum(pred & gt)
    fp = np.sum(pred & ~gt)
    fn = np.sum(~pred & gt)
    return (fp + fn) / (tp + fn), tp / (tp + fp + fn)


def training_slices(img, mask, n=5):
    """2D YX, ZX and ZY slices with a dummy third channel, as in the notebooks."""
    slices, labels = [], []
    for p in [(0, 1, 2, 3), (1, 0, 2, 3), (2, 0, 1, 3)]:
        img0 = img.transpose(p)
        mask0 = mask.transpose(p[:3])
        for i in np.linspace(5, img0.shape[0] - 5, n).astype(int):
            slices.append(np.concatenate((img0[i], img0[i][:, :, :1] * 0), axis=-1).transpose(2, 0, 1))
            labels.append(mask0[i].astype("uint16"))
    return slices, labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join("models", "SMG-BM-v0907"),
                    help="path to trained weights (default: models/SMG-BM-v0907)")
    ap.add_argument("--train", action="store_true",
                    help="fine-tune a throw-away model on the test volume first")
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--cpu", action="store_true", help="run on the CPU (slow)")
    a = ap.parse_args()

    c1 = io.imread(os.path.join(DATA, "test-c1-rcan.tif"))
    c2 = io.imread(os.path.join(DATA, "test-c2-rcan.tif"))
    gt = io.imread(os.path.join(DATA, "test-mask-annotated.tif"))
    img = transforms.normalize_img(np.stack([c1, c2], axis=-1), axis=-1)
    print("test volume", img.shape)

    if a.train:
        from cellpose import train
        model = models.CellposeModel(gpu=not a.cpu)
        slices, labels = training_slices(img, gt)
        model_path, _, _ = train.train_seg(model.net, train_data=slices, train_labels=labels,
                                           nimg_per_epoch=len(slices), min_train_masks=0,
                                           n_epochs=a.epochs, batch_size=1, learning_rate=1e-5,
                                           weight_decay=0.1, normalize=False, save_path=HERE,
                                           model_name="test-finetune")
        print("training ran; model saved to", model_path)
        print("loading the fine-tuned model (about 1.2 GB) ...", flush=True)
        model = models.CellposeModel(gpu=not a.cpu, pretrained_model=str(model_path))
    else:
        if not os.path.exists(a.model):
            raise SystemExit("model not found: %s (download it from Figshare, or use --train)" % a.model)
        print("loading", a.model, "...", flush=True)
        model = models.CellposeModel(gpu=not a.cpu, pretrained_model=a.model)

    print("segmenting the test volume in 3D (no progress output; a minute or two on a GPU,"
          " much longer on the CPU) ...", flush=True)
    seg = segment(model, img)
    io.imsave(os.path.join(DATA, "test-seg.tif"), seg)
    err, iou = scores(seg, gt)
    print("voxel error rate %.4f, IoU %.4f (written test_data/test-seg.tif)" % (err, iou))


if __name__ == "__main__":
    main()
