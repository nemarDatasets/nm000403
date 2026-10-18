# nm000403 stimulus-level events

Within-epoch event rows of the epoched recordings of nm000403, moved here on 18 Oct 2026 from the raw `events.tsv`, which now holds one `sequence_start` row per sequence (onset (k-1) x EpochLength, duration EpochLength, sequence k) so that `mne_bids.read_epochs_bids` loads each file as its epochs.

- 56 files, 42199 rows; same file names, columns, line text and order as the raw files before the move; `onset`/`sample` are relative to the start of the recording.
- Rows identical in every column to an earlier row were dropped (33755 rows).
- 0 rows lie outside the recording or outside the sequence they name (kept as delivered).
- `stimulus_events_summary.tsv`: per raw file, epochs, rows before, duplicates dropped, rows here.
