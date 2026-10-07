"""c-jonas (Dryad doi:10.5061/dryad.5f9v7; Jonas et al. 2016 PNAS) -> iEEG-BIDS DERIVATIVE dataset.

Source: 3 zips; per participant P<nn>/ two Letswave 5 pairs (.lw5 header + v7.3 .mat data 'data' float32, MATLAB dims
[sequences x channels x 1 x 1 x 1 x time]) for conditions face_periodic and face_nonperiodic; READ_ME.txt; Stimuli/ (JPEG).
README: data were imported from Micromed TRC, LOW-PASS FILTERED AT 30 Hz (Butterworth order 4) and segmented into
sequences (2 s after sequence onset to ~65 s); Letswave history also shows LW_crop and, for some files, LW_downsample.
-> derivative. Values are float32 in the source and are written unchanged as BrainVision IEEE_FLOAT_32, sequences
back-to-back (one segment per sequence). Letswave events -> events.tsv.
Privacy: some Letswave event codes start with a patient code (letters of the name, e.g. 'XXXXXX1_greyscale_...');
that prefix is removed in events.tsv and in the de-identified .lw5 copies in sourcedata (tokens never printed);
Letswave history dates and MAT text-header dates -> day 01. Stimuli/ images are NOT redistributed (the article says
the actual face images could not be shown for copyright reasons); they remain available from Dryad.
README-mentioned per-participant coordinate .txt files are absent from the deposit: labels only, no coordinates.
Usage: python b3w3_convert_jonas.py <sourcedata_dl> <bids_root>
"""
import hashlib, io, json, os, re, shutil, sys, zipfile
import numpy as np, h5py, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w3_common import wtsv, wjson, write_vhdr, scrub_mat_header, sha256_file

SRC, OUT = sys.argv[1], sys.argv[2]
SD = os.path.join(OUT, "sourcedata", "dryad-5f9v7-deidentified"); os.makedirs(SD, exist_ok=True)
report = {"files": [], "sourcedata": [], "excluded": []}
TOKENS = set()
PREFIX = re.compile(r"^([A-Za-z]{3,}\d*)_(?=greyscale)")
DATE = re.compile(r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})")


def scrub_code(c):
    m = PREFIX.match(c)
    if m:
        TOKENS.add(m.group(1))
        return "seq_" + c[m.end():]
    return c


def chtype(l):
    if re.match(r"^ECG", l): return "ECG"
    if re.match(r"^s[A-Z]", l): return "EEG"
    if re.match(r"^(PULS|BEAT|SpO2)", l): return "MISC"
    if re.match(r"^MKR", l): return "TRIG"
    return "SEEG"


scans = {}
for zn in sorted(x for x in os.listdir(SRC) if x.endswith(".zip")):
    z = zipfile.ZipFile(os.path.join(SRC, zn))
    for n in sorted(x for x in z.namelist() if x.endswith(".lw5")):
        pid, cond = re.match(r"(P\d+)/P\d+_face_(periodic|nonperiodic)_\d+sequences\.lw5$", n).groups()
        h = sio.loadmat(io.BytesIO(z.read(n)), squeeze_me=True, struct_as_record=False)["header"]
        ds = [int(v) for v in np.ravel(h.datasize)]; ne, nch, ns = ds[0], ds[1], ds[5]
        fs = 1.0 / float(h.xstep); x0 = float(h.xstart)
        labels = [str(c.labels) for c in np.atleast_1d(h.chanlocs)]
        with h5py.File(io.BytesIO(z.read(n[:-4] + ".mat")), "r") as hf:
            a = hf["data"][()]  # (time,1,1,1,ch,epochs)
        assert a.shape == (ns, 1, 1, 1, nch, ne), (a.shape, ds)
        assert a.dtype == np.float32
        x = np.concatenate([a[:, 0, 0, 0, :, e] for e in range(ne)], axis=0)  # (ne*ns, nch)
        S = f"sub-{pid}"; D = os.path.join(OUT, S, "ieeg"); os.makedirs(D, exist_ok=True)
        task = "face" + cond
        stem = f"{S}_task-{task}_run-1"
        np.ascontiguousarray(x).astype("<f4").tofile(os.path.join(D, stem + "_ieeg.eeg"))
        ev = []
        for e in np.atleast_1d(getattr(h, "events", [])):
            ep = int(e.epoch); lat = float(e.latency)
            s0 = int(round((lat - x0) * fs))
            onset = (ep - 1) * ns / fs + (lat - x0)
            ev.append([round(onset, 6), 0, scrub_code(str(e.code)), ep, round(lat, 6), (ep - 1) * ns + s0])
        ev.sort(key=lambda r: r[0])
        segs = [("New Segment", f"sequence {k + 1}", k * ns, 1) for k in range(1, ne)]
        write_vhdr(os.path.join(D, stem + "_ieeg"), nch, fs, labels, ["µV"] * nch,
                   comment=f"b3w3_convert_jonas.py: Letswave {os.path.basename(n)} + .mat, {ne} sequences back-to-back, float32 unchanged", markers=segs)
        wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "sequence", "latency_in_sequence", "sample"], ev)
        typ = [chtype(l) for l in labels]
        wtsv(os.path.join(D, stem + "_channels.tsv"), ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "group", "status", "status_description", "description"],
             [[l, t, "µV", "n/a", 30.0, fs, re.sub(r"[ ]?\d+$", "", l) if t == "SEEG" else "n/a", "good", "n/a",
               "label from the Letswave header (chanlocs); no coordinates released" if t == "SEEG" else "non-SEEG channel recorded with the SEEG (label from the Letswave header)"]
              for l, t in zip(labels, typ)])
        hist = [str(getattr(e, "description", "")) for e in np.atleast_1d(h.history)]
        wjson(os.path.join(D, stem + "_ieeg.json"), {
            "TaskName": task,
            "TaskDescription": ("Fast periodic visual stimulation: natural images of objects at 6 Hz (sinusoidal contrast modulation) with a face image as every fifth stimulus (1.2 Hz), sequences of ~70 s (article)."
                                if cond == "periodic" else
                                "Control condition: the same images at 6 Hz with faces presented non-periodically (article)."),
            "SamplingFrequency": fs, "PowerLineFrequency": 50,
            "SoftwareFilters": {"Letswave": {"Description": "Low-pass 30 Hz, Butterworth order 4 (README); Letswave history of this file: " + ", ".join(hist)}},
            "HardwareFilters": "n/a",
            "iEEGReference": "Midline prefrontal scalp electrode (FPz) in 21 participants or an intracerebral white-matter contact in 7 (article; not stated per participant)",
            "Manufacturer": "Micromed",
            "RecordingType": "epoched", "EpochLength": ns / fs,
            "RecordingDuration": ne * ns / fs,
            "SEEGChannelCount": typ.count("SEEG"), "EEGChannelCount": typ.count("EEG"), "ECGChannelCount": typ.count("ECG"),
            "MiscChannelCount": typ.count("MISC"), "TriggerChannelCount": typ.count("TRIG"), "ECOGChannelCount": 0,
            "ElectricalStimulation": False,
        })
        scans.setdefault(S, []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a"])
        import mne
        rr = mne.io.read_raw_brainvision(os.path.join(D, stem + "_ieeg.vhdr"), preload=False, verbose="error")
        k = min(2000, x.shape[0])
        rt = bool(rr.n_times == x.shape[0] and rr.ch_names == labels and np.allclose(rr.get_data(start=0, stop=k) * 1e6, x[:k].T, rtol=1e-6, atol=1e-6)
                  and np.allclose(rr.get_data(start=x.shape[0] - k) * 1e6, x[-k:].T, rtol=1e-6, atol=1e-6))
        report["files"].append(dict(file=n, zip=zn, sub=pid, task=task, sequences=ne, samples_per_sequence=ns, sfreq=fs, nch=nch, n_events=len(ev),
                                    types={t: typ.count(t) for t in set(typ)}, history=hist, roundtrip_ok=rt))
        print(stem, ne, nch, ns, fs, "events", len(ev), "rt", rt, flush=True)
for S, r in scans.items():
    wtsv(os.path.join(OUT, S, f"{S}_scans.tsv"), ["filename", "acq_time"], sorted(r))


def scrub_lw5(b):
    m = sio.loadmat(io.BytesIO(b), squeeze_me=False, struct_as_record=True)
    H = m["header"][0, 0]
    nchg = 0
    if "events" in H.dtype.names and H["events"].size:
        E = H["events"]
        for i in range(E.size):
            c = str(np.ravel(E.flat[i]["code"])[0]) if E.flat[i]["code"].size else ""
            nc = scrub_code(c)
            if nc != c:
                E.flat[i]["code"] = np.array([nc]); nchg += 1
    Hs = H["history"]
    for i in range(Hs.size):
        d = Hs.flat[i]["date"]
        if d.size:
            s = str(np.ravel(d)[0])
            ns_ = DATE.sub(lambda mm: f"01-{mm.group(2)}-{mm.group(3)}", s)
            if ns_ != s:
                Hs.flat[i]["date"] = np.array([ns_]); nchg += 1
    bio = io.BytesIO()
    sio.savemat(bio, {"header": m["header"]}, format="5", do_compression=True)
    nb, _ = scrub_mat_header(bio.getvalue())
    return nb, nchg


for zn in sorted(x for x in os.listdir(SRC) if x.endswith(".zip")):
    z = zipfile.ZipFile(os.path.join(SRC, zn))
    for i in z.infolist():
        if i.is_dir():
            continue
        if i.filename.startswith("Stimuli/"):
            report["excluded"].append(i.filename); continue
        b = z.read(i); o = hashlib.sha256(b).hexdigest()
        rel = os.path.join(zn[:-4], i.filename)
        if i.filename.endswith(".lw5"):
            nb, k = scrub_lw5(b); note = f"Letswave header rewritten (scipy savemat): {k} fields changed (event-code patient prefix removed, history dates day -> 01)"
        elif i.filename.endswith(".mat"):
            nb, ch = scrub_mat_header(b); note = "MAT text header date -> Mmm 01 yyyy" if ch else "unchanged"
        else:
            nb, note = b, "unchanged"
        t = os.path.join(SD, rel); os.makedirs(os.path.dirname(t), exist_ok=True)
        open(t, "wb").write(nb)
        report["sourcedata"].append([rel, i.file_size, o, hashlib.sha256(nb).hexdigest(), note])
for f in sorted(os.listdir(SRC)):
    if f.endswith(".txt"):
        shutil.copy(os.path.join(SRC, f), os.path.join(SD, f))
        report["sourcedata"].append([f, os.path.getsize(os.path.join(SRC, f)), sha256_file(os.path.join(SRC, f)), sha256_file(os.path.join(SD, f)), "unchanged"])
wtsv(os.path.join(SD, "DEIDENTIFICATION_MANIFEST.tsv"), ["path", "bytes_original", "sha256_original", "sha256_here", "change"], report["sourcedata"])
# token check over everything written (BIDS + sourcedata; lw5 decompressed via loadmat repr)
hits = 0
TB = [t.encode() for t in TOKENS]
for root, _, fs_ in os.walk(OUT):
    for fn in fs_:
        p = os.path.join(root, fn)
        if fn.endswith(".lw5"):
            txt = repr(sio.loadmat(p, squeeze_me=True, struct_as_record=False)["header"].__dict__).encode() + repr([str(e.code) for e in np.atleast_1d(getattr(sio.loadmat(p, squeeze_me=True, struct_as_record=False)["header"], "events", []))]).encode()
        elif fn.endswith((".eeg", ".mat")):
            continue
        else:
            txt = open(p, "rb").read()
        hits += sum(1 for t in TB if t in txt)
# data .mat files: raw byte search for the tokens
for root, _, fs_ in os.walk(SD):
    for fn in fs_:
        if fn.endswith(".mat"):
            b = open(os.path.join(root, fn), "rb").read()
            hits += sum(1 for t in TB if t in b)
import collections, glob
cc = collections.Counter()
for fe in glob.glob(os.path.join(OUT, "sub-*", "ieeg", "*_events.tsv")):
    for line in open(fe).read().splitlines()[1:]:
        cc[line.split("\t")[2]] += 1
report["event_codes_after_scrub"] = dict(cc)
print("event codes after scrub:", sorted(cc.items(), key=lambda kv: -kv[1]))
report["privacy"] = {"patient_code_tokens_removed": len(TOKENS), "token_hits_after_scrub": hits}
report["stimuli_excluded"] = len(report["excluded"])
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
for f in (__file__, os.path.join(os.path.dirname(os.path.abspath(__file__)), "b3w3_common.py")):
    shutil.copy(f, os.path.join(OUT, "code", os.path.basename(f)))
rep2 = dict(report); rep2["excluded"] = sorted({os.path.dirname(e) for e in report["excluded"]})
json.dump(rep2, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1)
print(json.dumps({"files": len(report["files"]), "rt_all": all(f["roundtrip_ok"] for f in report["files"]), "privacy": report["privacy"],
                  "stimuli_excluded": len(report["excluded"]), "subjects": len(scans)}, indent=1))
assert hits == 0, "patient-code tokens remain"
print("CONVERT_DONE")
