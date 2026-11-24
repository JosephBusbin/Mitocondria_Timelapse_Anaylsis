# Mito Construct Filter: Finding Mitochondria Movies by File Name

**Tool:** `scripts/construct_classifier.py`
**Settings file:** `config/construct_lexicon.json`
**Tests:** `scripts/test_construct_classifier.py`

This guide is for lab members. You only need to be able to copy commands into a
terminal (PowerShell on the lab PC, or bash on the cluster).

---

## 1. What it does and why

Our 4D time-lapses (`.sld` SlideBook files and `.tif` exports) are named in many
different ways. Extracting and analyzing a movie is slow, so we first want to know
**which movies actually show mitochondrial fluorescence**.

The filter reads **only the file and folder names**. It never opens a movie. For
each file, it looks for three kinds of clues:

1. **Mito words.** `mls` (Mitochondrial Localization Signal) is our main tag, as in
   `MAXmlsEGFP_T7LFeb26_45TL_16bit.tif`. `mito`, `cox8` and `tomm20` also count.
2. **Plasmid IDs from the lab plasmid table.** Names often glue IDs together with
   different zero-padding than the table uses:
   - `05feb2025obf01avc01-B.sld` contains `obf01` (= OBF001,
     `psox2bd::mitopaGFP_p2a-tagRFPt`, mito) and `avc01` (not mito).
   - `max_d318dec21jeb09jr19t1L_0-49.tif` contains `jeb09` (= JEB009, not mito) and
     `jr19` (= JR19, `pUAS-mlsEGFPpA`, mito, EGFP).

   `OBF001`, `obf01`, `obf0001` and `OBF01` are all treated as the same plasmid.
3. **Common mito markers** that are often used but not recorded in our table:
   TMRM, MitoTracker, JC-1, Drp1, TFAM, and others (see section 6). These are
   **never trusted automatically**. They send the file to human review.

### The four labels

| Label | Meaning | What to do |
|---|---|---|
| `MITO` | A mito-targeted fluorescent protein that is always bright (EGFP, mCherry, tagRFPT...) or a mito sensor (CEPIA). | Goes into the pipeline. |
| `MITO_PHOTOACTIVATABLE` | The only mito FP is photoactivatable or photoconvertible (paGFP, dendra2, mEos2, kaede...). | **Left out of the default list** (see below). |
| `NOT_MITO` | Plasmid IDs were recognized, none are mito, and there is no mito word. | Skipped. |
| `REVIEW` | Nothing recognized, an unknown ID, conflicting clues, or only a common marker. | A person checks it (see section 4). |

### Why photoactivatable movies are kept separate (the OBF001 example)

OBF001 is **mito-paGFP**. paGFP is almost invisible until it is activated with 405 nm
light. A movie of `05feb2025obf01avc01-B.sld` can therefore look like it has *no*
mitochondrial signal, except in the region that was activated. dendra2 and mEos2
start green and turn red after conversion. These movies need different analysis, so
they get their own label. If a file has both a normal mito FP and a
photoactivatable one, it is labelled `MITO` (the normal FP wins) and the reason says so.

### Red-mito warning

The pipeline assumes **channel 1 is the mitochondria**. In some constructs the mito
FP is red, not green. For example, in MF002 to MF010
(`pCMV_EGFP_p2a-mitoTAGRFPT`), **tagRFPT is the mito FP and EGFP is in the cytosol**.
For these files, the filter sets `mito_color = red` and `red_mito_flag = TRUE`. Check the
channel order before running the pipeline on them.

---

## 2. Running it on the Windows lab PC (PowerShell)

1. Open **PowerShell** and go to the repository folder (where `scripts\` is).
2. Activate the lab Python environment. This is the same environment used for
   extraction (see `windows_lab_extraction_guide.md`):

   ```powershell
   . "$env:LOCALAPPDATA\mito_pipeline_tools\activate_mito4d.ps1"
   ```

3. Run the filter on a data folder. The `` ` `` at the end of a line continues the
   command on the next line:

   ```powershell
   python scripts\construct_classifier.py `
       --plasmid-csv "D:\lab\Plasmid_Primer_List_Bestman_Lab.csv" `
       --scan-dir "D:\SlideBook_data" `
       --out "D:\mito_filter_results"
   ```

   - `--scan-dir` searches all subfolders. It picks up `.sld`, `.sldy`, `.tif`, `.tiff`,
     `.ome.tif`, `.nd2` and `.czi` files. It skips lock files (`~$...`) and hidden files.
   - To classify a typed-up list of names instead, use `--names-file names.txt`
     (one name per line; lines starting with `#` are ignored).
   - Add `--overrides "D:\mito_filter_results\overrides.csv"` to apply your corrections
     (section 5).
   - To avoid typing the plasmid table path every time, set it once per session:
     `$env:MITO_PLASMID_CSV = "D:\lab\Plasmid_Primer_List_Bestman_Lab.csv"`.

4. The screen shows a summary, for example:

   ```
   Classified 3 file(s):
     MITO                       2
     MITO_PHOTOACTIVATABLE      1
     NOT_MITO                   0
     REVIEW                     0   (0 high priority)
     red-mito (check channel order): 0
   ```

On the cluster (bash), the same command is:

```bash
module load python/3.12.7
python3 scripts/construct_classifier.py --plasmid-csv /path/to/plasmids.csv \
    --scan-dir /path/to/movies --out /path/to/results
```

The tool exits with code 0 on success. It exits with code 2 on bad input, for example
a missing file, a plasmid table without the expected columns, or a bad override label.

---

## 3. Reading `mito_classification.csv`

The output has one row per file. Open it in Excel.

| Column | Meaning |
|---|---|
| `rel_path` | File path relative to the scanned folder (forward slashes). |
| `label` | `MITO`, `MITO_PHOTOACTIVATABLE`, `NOT_MITO` or `REVIEW`. |
| `review_priority` | `high` or `normal` for `REVIEW` rows; empty otherwise. |
| `mito_fp` | The FP on the mitochondria (e.g. `EGFP`, `tagRFPT`, `paGFP`). Empty if unknown. |
| `fp_class` | `constitutive`, `photoactivatable`, `photoconvertible`, `sensor` or `unknown`. |
| `mito_color` | Channel the mito FP appears in: `green` (YFP counts as green), `red`, `far_red`, `blue`, `unknown`. |
| `red_mito_flag` | `TRUE` when the mito FP is red. Check the channel order. |
| `matched_ids` | Plasmids found, with their table names, e.g. `JR19=pUAS-mlsEGFPpA`. |
| `evidence` | Short codes, separated by `;` (table below). |
| `reasons` | A plain-English sentence explaining the label. |

Evidence codes:

| Code | Meaning |
|---|---|
| `mls_tag` | `mls` was found in the name. |
| `keyword:mito`, `keyword:cox8`, `keyword:tomm20` | A direct mito word was found. |
| `plasmid:JR19(mito,EGFP)` | Plasmid JR19 was found; it is mito with EGFP. |
| `plasmid:JEB009(no)` | Plasmid found; the table's localization is `no`. |
| `addgene:23348(mito,paGFP)` | A numeric (Addgene) ID was found. |
| `common_marker:tmrm` | A common mito marker was found. `(explained_by_plasmid)` means a recognized plasmid's own name contains it, e.g. `drp1` with jeb026 `pME_drp1k38a`. |
| `unknown_id:jeb99` | Looks like a lab ID, but it is not in the table. |
| `override:<pattern>` | Your override decided the label. |

---

## 4. Reading `mito_review.csv`

This file contains only the `REVIEW` rows, with **high priority first**.

- **High priority:** the only clue is a common mito marker (TMRM, MitoTracker, Drp1...),
  or a mito plasmid that has **no** fluorescent protein (e.g. 45161 `pcDNA3_Drp1K38A`).
  These are likely mito experiments, but a person must confirm what is fluorescent.
- **Normal priority:** nothing was recognized, an ID is not in the table (a typo, or a
  plasmid missing from the table), or the clues conflict. For example, `mito_jeb09.sld`
  has the word `mito` with no FP attached, but JEB009 is not mito.

For each row, open the movie (or ask whoever imaged it), then record the answer as an
override (section 5) so the next run is correct.

---

## 5. Writing overrides

An overrides file is a small CSV with three columns: `pattern,label,note`. Lines
starting with `#` are ignored. Example `overrides.csv`:

```csv
pattern,label,note
# whole experiment folder imaged with MitoTracker Green
2025-03-10_mtg_series/*,MITO,MitoTracker Green in ch1 (confirmed by JB)
# one specific movie
cell3_TMRM_2024.sld,MITO,TMRM red; check channel order
# a dead movie
*_bleached_*,NOT_MITO,photobleached before imaging
*jeb99*,REVIEW,ask who made jeb99
```

Rules:

- `pattern` uses wildcards: `*` = anything, `?` = one character. Upper/lower case does
  not matter.
- The pattern is compared to both the **relative path** (with `/`) and the **bare file
  name**.
- **The first matching row wins**, so put specific patterns above general ones.
- `label` must be exactly one of `MITO`, `MITO_PHOTOACTIVATABLE`, `NOT_MITO`, `REVIEW`.
  Any other value stops the run with an error, so typos are not silently ignored.
- The report keeps the FP columns from the name-based guess. `reasons` also says what the
  heuristic would have chosen.

---

## 6. The common-marker list and why it goes to REVIEW

Labs often image mitochondria with dyes or tagged proteins that never appear in the
plasmid table. A name that mentions them is a strong hint, but not proof. Examples:
"TMRM" might be a different channel from the one we analyze, and "drp1" might be a
construct with no fluorescent protein at all. So these markers always produce
`REVIEW` (high priority) and never `MITO` on their own.

| Group | Tokens matched in names |
|---|---|
| Dyes | `mitotracker`, `mtg`, `mtr`, `mtdr`, `cmxros`, `tmrm`, `tmre`, `jc-1`/`jc1`, `rhodamine123`/`rh123`, `mitosox`, `nao` |
| Reporters | `mitotimer`, `mt-keima`/`mito-keima`/`keima`, `mito-rogfp`, `mito-gcamp`, `mito-dsred`/`dsred2-mito`, `mito-bfp`, `su9`, `tom70`, `omp25`, `mts`, `mt` |
| Proteins | `tfam`, `pink1`, `parkin`, `mfn1`/`mfn2`, `opa1`, `drp1`/`dnm1l`, `miro`/`rhot1`, `mcu`, `vdac`, `hsp60`, `cytc`/`cycs`, `pgc1`/`ppargc1` |

The settings file records a one-line rationale for each token.

**Direct keywords versus common markers.** `mls`, `mito` (also the misspelling
`mit0`), `cox8`/`cox8a` and `tomm20`/`tom20` are **direct** keywords. They are the
targeting sequences our own constructs use, so they can produce `MITO`. Everything in the
table above is a common marker. When a reporter name contains the word "mito" (e.g.
`mito-GCaMP`, `mitoTimer`, `MitoTracker`), the common-marker rule wins and the file goes
to REVIEW.

**Drp1.** Our Drp1K38A plasmids (45161, jeb026, jeb027, bkd001) carry no mito FP.
- `drp1k38a_only.sld` → `REVIEW` (high).
- If the name also contains the matching plasmid ID (e.g. `jeb26_drp1k38a.sld`), the
  marker is "explained" by that plasmid and the file is `NOT_MITO`.

**Short tokens** (`mt`, `mts`, `mtr`, `mtg`, `mtdr`, `nao`, `miro`, `mcu`) only count
when they are not part of a longer word. Digits, `_`, `-` and spaces count as word
edges. So `05feb2025mt-x` matches, but `smtx`, `mtsat` and `smcurve` do not.

---

## 7. Adapting the filter for another lab (`config/construct_lexicon.json`)

All the vocabulary lives in the JSON file. You do not need to change the Python code.

1. **`columns`**: type the exact column headers of your plasmid spreadsheet (saved as
   CSV). You need at least the ID and localization columns.
2. **`localization_synonyms`**: map every value in your localization column to a
   canonical word. Only `mito` counts as mitochondrial.
3. **`direct_mito_keywords`**: the words your lab writes in mito movie names.
4. **`common_markers`**: add or remove dyes, reporters and proteins. Use
   `"boundary": "both"` for short tokens.
5. **`fp_vocabulary`**: every FP spelling you use. Give each one a `class`
   (constitutive / photoactivatable / photoconvertible / sensor) and a `color`.
6. **`numeric_id`**: the digit length of numeric IDs (Addgene numbers are 5–6 digits).
   This keeps dates such as `2025` or `20250205` from matching.
7. **`file_extensions`**: which files `--scan-dir` picks up.

After editing, run the tests:
`python scripts/test_construct_classifier.py -v`. You can also build a fake folder of
empty files with known answers and check the filter against it:

```bash
python3 scripts/construct_classifier.py --plasmid-csv plasmids.csv --make-synthetic scratch/fake_tree --seed 1
python3 scripts/construct_classifier.py --plasmid-csv plasmids.csv --scan-dir scratch/fake_tree --out scratch/fake_out
# compare scratch/fake_out/mito_classification.csv with scratch/fake_tree/expected_labels.csv
```

---

## 8. Limitations

- **The filter only reads file names.** A wrong or incomplete name gives a wrong label.
  A movie labelled `MITO` may still have no usable signal: the dish may not have
  expressed the construct, the movie may be photobleached, or the channels may be
  swapped. **Always confirm with the Stage 3 raw-signal QC** (in-cell vs background
  contrast in the mito channel) before trusting analysis results.
- The FP on the mitochondria comes from the construct name in the plasmid table (the
  FP right next to `mito`/`mls`/`cox8`/`tomm20`). If there is none, it falls back to
  Fluorescent Protein 1. Unusual table entries can be misread. Check `mito_fp` for new
  constructs.
- Plasmids that are missing from the table show up as `unknown_id`, or are not seen at all
  if their prefix is new. Add them to the table.
- Numeric Addgene IDs (5–6 digits) can, in rare cases, collide with other 5–6 digit
  numbers in names, such as times like `142530`. The `evidence` column shows
  `addgene:<id>` so you can spot this.
- It does not read SlideBook series names on its own. The Python function
  `classify_name(..., extra_names=[...])` accepts them for callers that do.

**Planned:** this filter will be wired into `scripts/collect_4d_sld_timelapses.py`
through `--plasmid-csv` / `--mito-filter` options, so collection can skip non-mito
movies automatically. Until then, run it on its own as shown above.
