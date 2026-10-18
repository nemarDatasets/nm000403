"""fix1018 fix_epochsplit.py <id> - epoched BrainVision files: events.tsv = exactly one row per epoch; every other row moves
unchanged to derivatives/stimulus_events/<same path> (same recipe as nm000487 on fix-20261017, accepted 17 Oct 2026).

Why: mne_bids 0.20 read_epochs_bids cuts one epoch of EpochLength at EVERY events.tsv row of an epoched continuous file. Files
whose events.tsv mixes epoch rows with within-epoch events (stimuli, triggers) or holds stimulus rows only fail ('Event time
samples were not unique', 'No matching events found', KeyError on a trial_type whose value is n/a), files with a header-only
events.tsv fail with 'max() iterable argument is empty', and the epochs it would cut at stimulus rows are not the file's epochs.

Number of epochs per file = n_samples (BrainVision header + .eeg size) / round(EpochLength x sfreq); must be an integer and equal
the .vmrk 'New Segment' count when the .vmrk has them. Modes:
  letswave: rows are Letswave events (onset, duration, trial_type, sequence, latency_in_sequence, sample); epoch rows are
            generated as in nm000487: trial_type sequence_start, onset (k-1) x EpochLength, duration EpochLength, sequence k,
            latency_in_sequence = SequenceStartTime (n/a when absent), sample (k-1) x n_per.
  rows:     the file already has one epoch row per epoch (trial_type == EPOCH_TT at sample (k-1) x n_per); they are kept as they
            are (a `value` column that is n/a on every kept row is dropped, else mne_bids builds no event id for them).
Only rows identical in every column to an earlier row are dropped (counted). Signal files unchanged."""
import os, sys, re, json, collections, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fx_common import Repo, headers, F

CONF = {
    "nm000403": dict(mode="letswave"),
    "nm000485": dict(mode="letswave"),
    "nm000498": dict(mode="letswave"),
    "nm000355": dict(mode="rows", tt="epoch"),
    "nm000533": dict(mode="rows", tt="epoch"),
    "nm000557": dict(mode="rows", tt="trial"),
}
LW = ["onset", "duration", "trial_type", "sequence", "latency_in_sequence", "sample"]


def segs_of(R, rec):
    v = rec[:-5] + ".vmrk"
    if not R.exists(v):
        return None
    t = R.read(v)
    return sum(1 for l in t.splitlines() if re.match(r"Mk\d+=New Segment,", l))


def fmt_on(x):
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    return s if s and s != "-0" else "0.0"


def run(nm):
    C = CONF[nm]
    R = Repo(nm)
    H = headers(nm)
    summ = []
    der = "derivatives/stimulus_events"
    for rec, r in sorted(H["recs"].items()):
        b = rec[: rec.rfind("_")]
        js = b + "_ieeg.json"
        J = R.json(js)
        if J.get("RecordingType") != "epoched":
            continue
        ev = b + "_events.tsv"
        h = r["hdr"]
        ns, sf = h["n_samples"], h["sfreq"]
        EL = float(J["EpochLength"])
        n_per = int(round(EL * sf))
        assert ns % n_per == 0, (rec, ns, n_per)
        n_ep = ns // n_per
        nseg = segs_of(R, rec)
        assert nseg in (None, 0, 1, n_ep), (rec, nseg, n_ep)   # 1 = the single New Segment BrainVision writes at the file start
        hdr, rows, nl, tr = R.tsv(ev)
        si = hdr.index("sample") if "sample" in hdr else None
        ti = hdr.index("trial_type") if "trial_type" in hdr else None
        rec_s = dict(file=ev, n_epochs=n_ep, samples_per_epoch=n_per, rows_before=len(rows), vmrk_segments=nseg)
        if C["mode"] == "letswave":
            assert hdr == LW, (ev, hdr)
            epoch_rows = [x for x in rows if x[2] == "sequence_start"]
            other = [x for x in rows if x[2] != "sequence_start"]
            xs = J.get("SequenceStartTime")
            new = []
            for k in range(1, n_ep + 1):
                on = (k - 1) * EL
                new.append([fmt_on(on), repr(EL), "sequence_start", str(k), repr(float(xs)) if xs is not None else "n/a", str((k - 1) * n_per)])
            if epoch_rows and not other and [x[:6] for x in epoch_rows] == new:
                rec_s.update(kind="already", rows_to_derivative=0, exact_duplicates_removed=0, rows_unmapped=0, rows_after=n_ep)
                summ.append(rec_s); continue
            assert not epoch_rows, ev
            new_hdr = hdr
        else:
            tt = C["tt"]
            epoch_rows = [x for x in rows if x[ti] == tt]
            other = [x for x in rows if x[ti] != tt]
            smp = [int(float(x[si])) for x in epoch_rows]
            assert smp == [k * n_per for k in range(n_ep)], (ev, smp[:5], n_ep, n_per, len(smp))
            new_hdr = list(hdr)
            new = [list(x) for x in epoch_rows]
            if "value" in hdr and all(x[hdr.index("value")] == "n/a" for x in new):
                vi = hdr.index("value")
                new_hdr = hdr[:vi] + hdr[vi + 1:]
                new = [x[:vi] + x[vi + 1:] for x in new]
                rec_s["value_column_dropped"] = True
        seen, keep = set(), []
        for x in other:
            k = "\t".join(x)
            if k in seen:
                continue
            seen.add(k); keep.append(x)
        unm = 0
        for x in keep:
            try:
                s = int(float(x[si])) if si is not None else int(round(float(x[hdr.index("onset")]) * sf))
            except ValueError:
                unm += 1; continue
            if not (0 <= s < ns):
                unm += 1
            elif C["mode"] == "letswave" and not ((int(float(x[3])) - 1) * n_per <= s < int(float(x[3])) * n_per):
                unm += 1
        rec_s.update(kind="split", rows_to_derivative=len(keep), exact_duplicates_removed=len(other) - len(keep), rows_unmapped=unm,
                     rows_after=len(new))
        R.write(ev, R.tsv_dump(new_hdr, new, nl, tr), f"events.tsv = {len(new)} epoch rows")
        if keep:
            R.write(f"{der}/{ev}", R.tsv_dump(hdr, keep, nl, tr), "stimulus rows (unchanged)")
        ej = b + "_events.json"
        if R.exists(ej):
            R.write(f"{der}/{ej}", R.read(ej), "events.json copy")
        R.targets.append(rec)
        summ.append(rec_s)
    split = [s for s in summ if s["kind"] == "split"]
    tot = dict(files=len(summ), split=len(split), epochs=sum(s["n_epochs"] for s in summ),
               rows_to_derivative=sum(s["rows_to_derivative"] for s in split),
               duplicates=sum(s["exact_duplicates_removed"] for s in split), unmapped=sum(s["rows_unmapped"] for s in split),
               header_only_before=sum(1 for s in split if s["rows_before"] == 0))
    R.notes.append(json.dumps(tot))
    # root events.json copies + epoch-row description
    for f in sorted(os.listdir(R.root)):
        if re.match(r"task-.*_events\.json$", f):
            R.write(f"{der}/{f}", R.read(f), "root events.json copy")
            if C["mode"] == "letswave":
                E = R.json(f)
                E["onset"] = {"Description": "Start of one sequence (epoch) in this file: (sequence - 1) x EpochLength (sequences back-to-back)", "Units": "s"}
                E["duration"] = {"Description": "EpochLength (length of the sequence)", "Units": "s"}
                E["trial_type"] = {"Description": "Every events.tsv holds exactly one row per sequence (epoch) of the data file (since 18 Oct 2026), "
                                   "so that readers of epoched data (mne_bids.read_epochs_bids) load each file as its sequences. The "
                                   "stimulus-level Letswave events are in derivatives/stimulus_events/ (README, 'Events and epochs').",
                                   "Levels": {"sequence_start": "Start of one sequence (epoch) of the file; not a stimulus event"}}
                E["latency_in_sequence"] = {"Description": "Time of the epoch start on the Letswave sequence time axis (= SequenceStartTime in ieeg.json when given; n/a otherwise)", "Units": "s"}
                E["sample"] = {"Description": "First sample of the sequence in this file (0-based) = (sequence - 1) x round(EpochLength x SamplingFrequency)"}
                R.write_json(f, E, "events.json: sequence_start rows")
    dd = R.json("dataset_description.json")
    R.write_json(f"{der}/dataset_description.json", {
        "Name": f"{nm} stimulus-level events (moved out of the epoch-defining events.tsv, 18 Oct 2026)",
        "BIDSVersion": dd.get("BIDSVersion", "1.10.0"), "DatasetType": "derivative", "License": dd.get("License", "n/a"),
        "Authors": dd.get("Authors", []),
        "GeneratedBy": [{"Name": f"{nm} fix-20261018 events split", "Description":
            f"Within-epoch event rows of the epoched sub-*/*_events.tsv files of {nm}, moved here unchanged (same columns, line text and "
            "order; onset/sample relative to the start of the recording) after dropping rows identical in every column to an earlier row. "
            "The raw events.tsv now holds one row per epoch so that mne_bids.read_epochs_bids loads each file as its epochs. Code: "
            "../../code/fix-20261018_split_stimulus_events.py; per-file counts: stimulus_events_summary.tsv."}],
        "SourceDatasets": [{"URL": f"https://github.com/nemarDatasets/{nm}", "Version": dd.get("Version", "n/a")}]}, "derivative dd")
    cols = ["file", "kind", "n_epochs", "samples_per_epoch", "vmrk_segments", "rows_before", "exact_duplicates_removed", "rows_to_derivative", "rows_unmapped", "rows_after"]
    R.write(f"{der}/stimulus_events_summary.tsv", "\t".join(cols) + "\n" + "".join("\t".join(str(s.get(c, "n/a")) for c in cols) + "\n" for s in summ), "summary")
    epoch_txt = ("one `sequence_start` row per sequence (onset (k-1) x EpochLength, duration EpochLength, sequence k)" if C["mode"] == "letswave"
                 else f"its `{C['tt']}` rows (one per epoch, already present, unchanged" + ("; the `value` column, n/a on every one of them, was dropped)" if any(s.get("value_column_dropped") for s in summ) else ")"))
    R.write(f"{der}/README.md", f"# {nm} stimulus-level events\n\nWithin-epoch event rows of the epoched recordings of {nm}, moved here on "
            f"18 Oct 2026 from the raw `events.tsv`, which now holds {epoch_txt} so that `mne_bids.read_epochs_bids` loads each file as its epochs.\n\n"
            f"- {tot['split']} files, {tot['rows_to_derivative']} rows; same file names, columns, line text and order as the raw files before the move; "
            "`onset`/`sample` are relative to the start of the recording.\n"
            f"- Rows identical in every column to an earlier row were dropped ({tot['duplicates']} rows).\n"
            f"- {tot['unmapped']} rows lie outside the recording" + (" or outside the sequence they name" if C["mode"] == "letswave" else "") + " (kept as delivered).\n"
            "- `stimulus_events_summary.tsv`: per raw file, epochs, rows before, duplicates dropped, rows here.\n", "derivative README")
    os.makedirs(R.p("code"), exist_ok=True)
    R.write("code/fix-20261018_split_stimulus_events.py", open(os.path.abspath(__file__)).read(), "code copy")
    R.summary = summ; R.tot = tot; R.conf = C
    return R


if __name__ == "__main__":
    R = run(sys.argv[1])
    print(R.tot)
    R.done()
