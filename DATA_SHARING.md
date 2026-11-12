# Data-sharing status

Last reviewed: 2026-09-05

The repository contains scientific measurements, derived MATLAB workspaces, legacy scripts,
figures, archives, and a separate Julia implementation. Public visibility does not itself
grant a license to reuse these materials.

## Current permission status

- The Julia implementation has no selected project license yet. Default copyright therefore
  applies until a root `LICENSE` file is added by the rights holder.
- No explicit license or institutional permission record is present for the thesis-era TGA
  measurements and outputs under `code/`. The MATLAB source is archived separately.
- The retained TGAnalysis data and results have their original attribution and license files.
  See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
- Generated Julia reports and figures inherit unresolved project/data rights until the code
  and source-data licensing decisions are recorded.

The repository is already configured with the public GitHub origin
<https://github.com/axeld92/isoconvAnalysis>. Before tagging or advertising a v1.0 release,
the maintainer should confirm that public redistribution of the thesis code, measurements,
instrument exports, and derived outputs is authorized.

## Scientific-data contents

The retained core data are non-human thermogravimetric measurements of waste tire,
*Acrocomia aculeata* endocarp, and dry-mass mixtures. No human-subject or clinical records were
identified. The 21 NETZSCH text exports do retain laboratory metadata (`UNA - FCQ - AI`),
operator initials (`AD`), acquisition dates, sample identifiers, instrument model, calibration
filenames, and crucible/sample masses. Confirm whether those fields may be published or should
be minimized in a separately licensed data deposit.

The retained data and results under `code/` include 32 MAT files, 71 PDFs, two RAR archives
of figures, and two XLSX files.

## Release choices requiring maintainer approval

1. **Code-only release:** license the Julia source, exclude `code/` from the curated release
   bundle, and document how authorized users supply data locally. GitHub's automatically
   generated source archive would still reflect whatever is tracked in the tagged commit.
2. **Code-and-data release:** confirm author/institutional rights, select a separate data
   license, add a machine-readable data citation, and state exactly which raw and derived
   artifacts the data license covers.
3. **Restricted-data release:** publish the Julia source and metadata/configuration, deposit
   only cleared data elsewhere, and provide checksums plus access instructions.

Do not infer a data license from the eventual software license. Software, original data,
third-party code, and generated outputs should each have an explicit scope.

## Required v1.0 clearance record

Before release, record:

- the software copyright holder(s) and chosen root license;
- who owns the measurements and whether the laboratory/institution authorizes redistribution;
- whether operator/laboratory metadata may remain in the public exports;
- the selected data license and exact covered paths, or the exclusion/access policy;
- disposition of the two figure archives and retained TGAnalysis data and results;
- citation metadata for the thesis, software release, and any data deposit.

This document records technical due diligence, not legal advice. Questions about ownership or
institutional policy require confirmation from the relevant rights holder or institution.
