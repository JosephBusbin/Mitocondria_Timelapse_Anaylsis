# Windows Lab Workstation Guide: 4D SlideBook (.sld) Extraction

This manual provides step-by-step instructions for extracting, validating, and staging 4D mitochondrial time-lapse data from proprietary 3i SlideBook (`.sld` / `.sldy`) microscope containers on a local Windows lab computer.

---

## Overview & Architecture

- **Goal**: Convert multi-gigabyte raw SlideBook containers into losslessly compressed 16-bit OME-TIFF hyperstacks (`.ome.tif`), extract dimensional metadata, compute SHA-256 transfer checksums, and prepare batch manifests for W&M SciClone HPC.
- **Zero Administrator Privileges**: The setup installs exclusively into user-space (`%LOCALAPPDATA%\mito_pipeline_tools`). It **never** modifies the Windows Registry, global system variables, or administrative directories.
- **Native 3i SlideBook 6/7 Support**: Uses 3i's official `SlideBook6Reader.jar` and `SBReadFile.dll` so you can read both monolithic `.sld` files and streamed `.sldy` datasets without requiring commercial SlideBook software or companion folders.

---

## Step 1: Open PowerShell and Navigate to the Repository

1. Press `Windows Key + R`, type `powershell`, and press `Enter`.
2. Navigate to your project folder:
   ```powershell
   cd "$env:USERPROFILE\path\to\Mitocondria_Timelapse_Anaylsis"
   ```
   *(Replace `path\to` with wherever you cloned the repository, e.g. `Documents\GitHub`.)*

---

## Step 2: Perform a Dry Run of the Environment Setup

Before downloading or configuring anything, verify what the setup script will do:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\setup_lab_windows.ps1 -DryRun
```

### Expected Output
The script will report:
- Whether an existing Java runtime (e.g., from FIJI) was found, or if the portable Temurin JRE 17 (~42 MB) will be downloaded.
- The Bio-Formats version (`8.5.0` by default) and CLI tools (`showinf.bat` / `bfconvert.bat`) target location.
- The `SlideBook6Reader.jar` target location (`...\mito_pipeline_tools\readers\`).
- Isolated Python 3.12 virtual environment destination (`mito4d_env`).
- Confirmation that no files were downloaded or modified (`[DRY-RUN]`).

---

## Step 3: Run the Full Environment Setup

Run the setup script to provision the isolated local extraction stack:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\setup_lab_windows.ps1
```

> [!TIP]
> If FIJI is installed on the computer (for example at `C:\Fiji.app`), you can specify `-FijiDir "C:\Fiji.app"` to reuse its bundled Java runtime directly:
> ```powershell
> powershell -ExecutionPolicy Bypass -File .\scripts\setup_lab_windows.ps1 -FijiDir "C:\Fiji.app"
> ```

### What Happens Automatically
1. **TLS 1.2 Enforcement**: Enables secure downloads from Adoptium, Open Microscopy Environment (OME), and GitHub.
2. **Directory Creation**: Creates `%LOCALAPPDATA%\mito_pipeline_tools`.
3. **Java Setup**: Reuses existing Java or downloads portable Eclipse Temurin JRE 17.
4. **Bio-Formats & 3i Reader**: Installs `bftools` 8.5.0 (upgrading an older 7.3.0 install in place) and downloads `SlideBook6Reader.jar` from 3i's official repository into `...\mito_pipeline_tools\readers\`. The jar lives outside `bftools\` so future upgrades never delete it; a jar left in `bftools\` by an older setup is moved automatically. The reader is still required: Bio-Formats' built-in SlideBook reader crashes on SlideBook 6 `.sld` files.
   - Pin a different version with `-BfToolsVersion 8.4.0`, or a newer 3i jar with `-SbReaderUrl <url>`.
5. **Python Virtual Environment**: Uses portable `uv` to create a lightweight Python 3.12 venv.
6. **Verification**: Generates `activate_mito4d.ps1` and runs verification checks.

Look for the final success banner:
```text
======================================================================
[SUCCESS] Local .sld extraction environment verified!
  Python      : PYTHON_OK:3.12.7
  Java Home   : ...\mito_pipeline_tools\jre\jdk-17.0.21+9-jre
  Bio-Formats : bftools 8.5.0 (...\mito_pipeline_tools\bftools)
  3i Reader   : ...\mito_pipeline_tools\readers\SlideBook6Reader.jar (reads SlideBook 6 .sld and 7 .sldy)
======================================================================
```

---

## Step 4: Activate the Extraction Environment

In your PowerShell window, dot-source the activation script:

```powershell
. "$env:LOCALAPPDATA\mito_pipeline_tools\activate_mito4d.ps1"
```

You should see:
```text
[OK] Activated mito4d_env | Java: ... | bftools 8.5.0 | 3i SlideBook6Reader
```

> [!IMPORTANT]
> Keep the leading **dot and space** (`. `). Without them, PowerShell just prints the quoted path back to you and nothing is activated. `$env:LOCALAPPDATA` expands to your own `C:\Users\<you>\AppData\Local`, so the same line works for every user.

> [!NOTE]
> This activation sets `JAVA_HOME`, `BF_JAVA_HOME`, `BF_CP`, and updates `PATH` **only for the current PowerShell window**. It leaves all system-wide settings completely untouched.

---

## Step 5: Perform a Dry-Run Extraction

Test discovery on your raw SlideBook data without writing multi-gigabyte images to disk:

```powershell
python scripts\collect_4d_sld_timelapses.py "data\practice" "data\test_output" --mode extract-ometiff --dry-run
```

### Expected Output Summary
```text
Using 3i SlideBook native reader JAR: ...\readers\SlideBook6Reader.jar
Found 1 valid SlideBook container(s) (7.67 GB total raw container size).
Filtering criteria: SizeZ >= 15, SizeT >= 10, SizeC >= 2
Using Java runtime: ...\jre\jdk-17.0.21+9-jre

--- Scan & Manifest Summary ---
Total .sld files scanned          : 1 (7.67 GB)
Unparseable or skipped files      : 0 of 1 candidate(s)
Total internal series inspected   : 49
Qualifying 4D time-lapse series   : 1 across 1 .sld file(s)
Estimated uncompressed 4D data    : 1.15 GB (~0.63 GB with lossless LZW)
Inventory CSV saved to            : data\test_output\sld_4d_inventory.csv
File health audit saved to        : data\test_output\sld_health_audit.csv
SLURM array list saved to         : data\test_output\movie_list.txt (relative paths; resolved against $MITO_RAW_ROOT)
Checkpoint 2 provenance saved to  : data\test_output\extraction_provenance.json

--- File Health by Cohort Year ---
  2025: total=1 (OK_4D=1)

Extracting 1 qualifying 4D series to LZW-compressed .ome.tif...
  [OK] [DRY-RUN] Would extract series 20 -> data\test_output\05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif

Completed extraction: 0.00 GB written to data\test_output
```

Notice that Series #20 ($Z=49, T=10, C=2$) was identified as qualifying 4D data, while non-qualifying 2D snaps and projections were safely skipped.

---

## Step 6: Perform Full Extraction with SHA-256 Checksums

Execute the real extraction to convert the qualifying 4D series into compressed OME-TIFF files:

```powershell
python scripts\collect_4d_sld_timelapses.py "data\practice" "data\test_output" --mode extract-ometiff --checksum
```

### Expected Output
```text
Extracting 1 qualifying 4D series to LZW-compressed .ome.tif...
  [OK] Extracted: 05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif (0.456 GB)

Completed extraction: 0.46 GB written to data\test_output
```

The script will write the following files to `data\test_output`:
1. `05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif` (the extracted 4D hyperstack)
2. `extracted_checksums.sha256` (cryptographic SHA-256 hash list for HPC transfer auditing)
3. `sld_4d_inventory.csv` (complete metadata inventory for all 49 internal series)
4. `sld_health_audit.csv` (one row per `.sld` file found: status, size, header bytes, error summary)
5. `movie_list.txt` & `movie_list_2025.txt` (SLURM job array manifests for SciClone HPC)
6. `extraction_provenance.json` (machine, OS, tool version, reader jar, and parameter audit record)

> [!NOTE]
> `movie_list.txt` entries are written **relative** to the raw OME-TIFF root (e.g. `2025\...ome.tif` becomes `2025/...ome.tif`) so the same list works wherever the data lands on SciClone. The HPC job resolves them against `$MITO_RAW_ROOT`. To write absolute paths instead, pass `--hpc-raw-root /sciclone/data10/<storage>/raw_ometiff` (or set `MITO_RAW_ROOT` before running).

---

## Step 7: Verify Hyperstack & Checksum Integrity

### A. Inspect Hyperstack Dimensions with Bio-Formats
Run `showinf.bat -nopix` to confirm that all 980 planes and dimensional calibrations are intact:

```powershell
showinf.bat -nopix "data\test_output\05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif"
```

Look for these critical lines in the output:
- `Checking file format [OME-TIFF]`
- `Image count = 980` ($49 \times 10 \times 2$)
- `Width = 952`, `Height = 664`
- `SizeZ = 49`, `SizeT = 10`, `SizeC = 2`
- `Dimension order = XYCZT (certain)`
- `Pixel type = uint16`
- `Metadata complete = true`

### B. Verify the SHA-256 Checksum
Compute the hash with PowerShell and compare it with the entry in `extracted_checksums.sha256`:

```powershell
Get-FileHash -Algorithm SHA256 "data\test_output\05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif" | Format-List Hash
Get-Content "data\test_output\extracted_checksums.sha256"
```

The hash output by `Get-FileHash` must match the corresponding line in `extracted_checksums.sha256`:
```text
Hash : CF0F3644A2B30743D3C7C9CAB9027AC4E406168DBCB690F871D61B19BAD90F6C

extracted_checksums.sha256:
cf0f3644a2b30743d3c7c9cab9027ac4e406168dbcb690f871d61b19bad90f6c  05feb2025obf01avc01-B__S20_05feb2025avc01obf01-B_t2r_-_4.ome.tif
```
*(PowerShell displays hexadecimal letters in uppercase, while sha256 manifests use lowercase; they are identical).*

> [!NOTE]
> **Why do hashes differ between fresh extractions?**
> The OME-TIFF standard requires Bio-Formats to embed a unique UUID (`UUID="urn:uuid:..."`) into the OME-XML header of each newly created container. Because this UUID is randomly generated upon each extraction run, a newly generated TIFF will have a unique whole-file SHA-256 hash. The pixel buffers, z-slices, and dimensions remain 100% bitwise identical.
>
> **The real purpose of `extracted_checksums.sha256`:**
> Once you transfer or `rsync` your extracted `.ome.tif` files from the lab PC to SciClone HPC, run:
> ```bash
> sha256sum -c extracted_checksums.sha256
> ```
> on the HPC to verify that **zero bits were corrupted or dropped during network transit**.

---

## Step 8: Run the Test Suite

Run the automated unit and integration tests to verify your environment:

```powershell
python scripts\test_pipeline_suite.py -v
```

### Expected Result
```text
test_01_sld_metadata_and_safety: ok
test_01b_criteria_thresholds: ok
test_01c_companion_dirs: ok
test_01d_sanitize_and_size_estimation: ok
test_01e_live_sld_inspection_if_available: ok
test_01f_checksum_generation_and_verification: ok
test_01g_extract_ometiff_dry_run: ok
test_01h_skip_existing_files: ok
test_01i_health_audit_prefilter_statuses: ok
test_01j_error_classification_uses_raw_stderr: ok   (the fake-showinf half is skipped on Windows)
test_01k_movie_list_entries_relative_or_absolute: ok
... (tests 02-02e skipped: require numpy + tifffile) ...
... (tests 03-08 skipped: require HPC analysis packages) ...

Ran 22 tests
OK (skipped=12)
```

---

## Advanced Options & Common Scenarios

| Scenario | Command Flag | Example Use Case |
| :--- | :--- | :--- |
| **Short Time-Lapses ($T < 10$)** | `--min-t 5` | Extracts shorter acquisitions (e.g., Series #39 with $T=6$). |
| **Single-Channel Acquisitions** | `--min-c 1` | Extracts movies with only one fluorophore (e.g., mito-GFP alone, Series #30). |
| **Smaller Z-Stacks** | `--min-z 10` | Permits thinner z-stacks. |
| **Alternative Compression** | `--compression Deflate` | Deflate (gzip) compression instead of standard LZW. |
| **Re-extract / Overwrite** | `--overwrite` | Forces re-extraction of previously extracted files. *(By default, existing files are safely skipped in seconds).* |
| **Copy Raw Files Instead** | `--mode copy-sld` | Stages the raw `.sld` container and companion folders instead of converting to TIFF. |

---

## Step 9: Understanding & Troubleshooting Lab Archive Warnings

When running `collect_4d_sld_timelapses.py` on real laboratory hard drives containing years of microscopy data, you may see warnings such as:
```text
[WARN] [SLIDEBOOK_IMAGE_RECORDS_ERR] Skipping aborted_experiment.sld (0.1 MB): Corrupted Image Records / aborted acquisition (likely aborted capture or companion folder missing)
```
or a Bio-Formats exit code 1 with:
```text
Exception in thread "main" java.lang.AssertionError: Failed with exception: The file is not a valid SlideBook document. Problem reading the Image Records.
    at loci.formats.in.SlideBook6Reader.openFile(Native Method)
```

### Is This an Error in the Script?
**No.** The script is working exactly as designed. The warning is the error-handling guard doing its job: trapping unreadable or corrupted files so that a single bad acquisition does **not** crash your multi-hour batch extraction.

### What Causes This on a Lab Microscope PC?
3i's native C++ reader (`SBReadFile.dll`) throws this assertion when the file header exists, but the internal **Image Records table** (the metadata catalog indexing acquisitions, channels, and planes) cannot be read. On physical lab PCs, this typically stems from four real-world situations:

1. **Aborted Acquisitions (Most Common)**:
   When an experiment is canceled mid-run, if SlideBook crashes, if the laser shutter fails, or if the drive fills up, SlideBook writes the initial file header but never finalizes the Image Records table. The file remains on disk as a truncated stub.
2. **Orphaned Streamed Containers (Missing Companion Folders)**:
   When recording in "Streamed to Disk" mode, SlideBook saves metadata in `filename.sld` and raw pixel records in an adjacent companion folder (`filename_data` or `filename.sld_Files`). If a user copied only the `.sld` file and left the companion folder behind on another drive, `SBReadFile.dll` cannot find the image records.
3. **Zero-Byte Stubs & Temporary Lock Files**:
   Microscope workstations accumulate temporary lock files (e.g., `.~filename.sld`, `~$filename.sld`) or 0-byte files from crashed sessions. The pipeline automatically pre-filters files `< 64 KB` and temp prefixes to save processing time.
4. **SlideBook Is Currently Open**:
   If SlideBook 6 is running with the file open, Windows places a write lock on the file, preventing the Bio-Formats reader from accessing the Image Records.

### How to Inspect Unreadable Files
The script records **every** `.sld`/`.sldy` file it finds, readable or not, in one audit spreadsheet. This replaces the older `unparseable_sld_files.csv`.
```text
data\test_output\sld_health_audit.csv
```
Columns:
- `cohort_year`, `rel_path`, `sld_filename`: where the file sits in the archive
- `file_size_mb`, `modified_utc`: size and last-modified time
- `has_companion_dirs`: whether a `_data` / `.sld_Files` / `_Files` / `.dir` folder sits next to it
- `magic_hex`, `is_null_header`: the first 16 bytes, and whether the first 4 KB are all zeros (a pre-allocated file from a crashed capture)
- `status`: one of the values below
- `total_series`, `qualifying_4d_series`: series counts when the file could be read
- `error_summary`: the first 300 characters of the original Java error

| `status` | Meaning | What to do |
| :--- | :--- | :--- |
| `OK_4D` | Read successfully; contains qualifying 4D series | Nothing |
| `OK_NON_4D` | Read successfully; no series meet the Z/T/C criteria | Nothing (or relax `--min-*`) |
| `EMPTY_NO_IMAGES` | Readable container with zero images | Usually an empty session file |
| `SKIPPED_STUB` | Smaller than 64 KB | Aborted or placeholder file |
| `SKIPPED_TEMP_OR_LOCK` | Name starts with `~`, `.~`, `._`, `.#`, or `temp` (never opened) | Ignore; close SlideBook if it is still running |
| `UNREADABLE_FILE` | Windows denied access to the file | Check permissions or locks |
| `EMPTY_OR_NULL_FILE` | First 4 KB are all zeros | Crashed acquisition; the data is unrecoverable |
| `SLIDEBOOK_IMAGE_RECORDS_ERR` | 3i reader could not read the Image Records table | See causes 1, 2, and 4 above |
| `MISSING_SB6_READER` | Bio-Formats fell back to its built-in reader and crashed | `SlideBook6Reader.jar` was not found: re-run setup, then activate |
| `JNI_OR_JAVA_ERR` | Native library, class-path, or out-of-memory failure | Re-run setup; try `--workers 1` |
| `SLIDEBOOK_SESSION_CONCURRENCY_ERR` | 3i reader reported a file open/close conflict | Re-run with `--workers 1` |
| `OME_XML_PARSE_ERR` | showinf output could not be parsed | Report it, with the file name |
| `OTHER_ERR` | Anything else | Read `error_summary` |

The same counts appear per cohort year at the end of the console output and under `health_audit_by_year` in `extraction_provenance.json`.

To spot-check any suspect file:
1. Open PowerShell and check its size:
   ```powershell
   Get-Item "C:\Path\To\file.sld" | Select-Object Name, Length
   ```
2. Try opening the file directly in the official 3i SlideBook software on the lab PC. If SlideBook itself fails to open the document, the file is confirmed damaged or aborted.
3. Check whether an adjacent companion folder (`_data` or `.sld_Files`) exists on an external drive or backup.
