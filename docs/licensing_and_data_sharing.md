# M9 licensing and data-sharing review

Last reviewed: 2026-09-05

The inventory below describes the original M9 review snapshot. MATLAB source was subsequently
removed from the project tree. The retained measurements and outputs are covered by
[`DATA_SHARING.md`](../DATA_SHARING.md).

## Review decision

The locked direct Julia dependencies present no identified software-license conflict: all
five third-party runtime packages, three documentation/development packages, and the
test-only Aqua package carry MIT licenses at the resolved versions. Package sources and
binary artifacts are not vendored into the Julia implementation.

The repository itself is **not ready for an unrestricted v1.0 distribution**. It has no root
software license, no separate data license, no `CITATION.cff`, and no recorded confirmation
that the thesis-era measurements and legacy author code may be redistributed. GitHub's
[licensing guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository)
states that default copyright applies when no license is present; making a repository public
does not grant general permission to reproduce, distribute, or create derivatives.

This is a release blocker requiring a maintainer/rights-holder decision. It does not block
local validation or continued documentation work.

## Direct dependency review

The review inspected the license files installed for the exact versions in `Manifest.toml`,
`test/Manifest.toml`, and `docs/Manifest.toml`. Full rationale remains in
[`dependencies.md`](dependencies.md).

| Environment | Direct third-party packages | Result |
|---|---|---|
| Runtime | MAT 0.12.1, Distributions 0.25.130, Optim 2.2.1, SciMLBase 3.50.2, OrdinaryDiffEqTsit5 2.1.4 | MIT at each locked version |
| Test | Aqua 0.8.16 | MIT |
| Documentation/development | CairoMakie 0.15.13, Documenter 1.17.0, JuliaFormatter 2.12.5 | MIT at each locked version |

An automated review of resolved package-source license files across all three environments
found MIT-like notices for the overwhelming majority, one BSD-like notice, and several Julia
standard-library or plotting-support notices that require manual classification. This scan is
not a substitute for the upstream notices embedded in downloaded JLL binary artifacts. That
remaining obligation is immaterial to a source-only release, but a compiled application or
bundled Julia depot would require a generated third-party-notice inventory.

## Vendored third-party code

`code/github_repo/` and `code/TGAnalysis-master/TGAnalysis-master/` are duplicate archive
copies of [mayankskii/TGAnalysis](https://github.com/mayankskii/TGAnalysis). Both retain the
same MIT license and copyright notice. They are unrelated to the Julia runtime and should be
excluded from a curated v1.0 validation bundle unless retaining the historical archive is an
explicit release goal. If kept, the nested license files must remain. The root summary is in
[`../THIRD_PARTY_NOTICES.md`](../THIRD_PARTY_NOTICES.md).

## Tracked-data and repository-size review

The current release-candidate tree contains 908 non-ignored files totaling about 356.10 MiB;
792 files and about 346.06 MiB are under `code/`. The local compressed Git object store is
about 262.07 MiB.

| Category | Files | Size (MiB) | Release interpretation |
|---|---:|---:|---|
| Julia project, configuration, tests, and docs | 116 | 10.03 | Candidate software release material |
| Instrument text exports | 21 | 13.44 | Core provenance; redistribution permission unresolved |
| MAT data outside TGAnalysis | 30 | 123.44 | Core plus historical/derived workspaces; scope must be selected |
| Legacy result/article outputs | 304 | 136.02 | Historical context; not required by the Julia runtime |
| Other legacy files | 163 | 6.33 | Mixed author scripts and archive context |
| Duplicate TGAnalysis trees | 274 | 66.84 | Third-party MIT archive; exclude or retain with notice |

The tracked binary extensions include 32 MAT files (169.71 MiB), 71 PDFs (12.53 MiB), two RAR
archives (43.46 MiB), and two XLSX files. The largest single file is
`code/resultados.rar` at about 38.60 MiB, so no current object exceeds GitHub's 50 MiB warning
threshold. The aggregate still makes source archives and clones unnecessarily heavy.
GitHub's [large-file guidance](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
recommends keeping repositories small, using Git LFS for large binaries, and distributing
release binaries outside normal Git history.

No `.gitattributes` LFS rules are configured. The existing rule `code/** -text` is valuable
and must remain because it prevents line-ending normalization from changing legacy evidence.

## Privacy and scientific metadata

The inspected TGA exports contain material measurements, not human-subject data. They do
include laboratory identity, operator initials, acquisition timestamps, sample identifiers,
instrument model, calibration filenames, and crucible/sample masses. The privacy risk appears
low, but institutional ownership and publication policy cannot be inferred from file content.
The maintainer should confirm permission and whether operator/laboratory fields should be
retained in a public data deposit.

## Clearance matrix

| Item | Status | Required action before v1.0 |
|---|---|---|
| Direct Julia dependency licenses | Reviewed; compatible | Preserve manifests; retain license evidence |
| Transitive source dependencies | Screening complete | Recheck when manifests change |
| Bundled binary-runtime notices | Not applicable to current source-only layout | Generate notices if distributing a compiled depot/application |
| TGAnalysis third-party archive | MIT notice present | Exclude from curated bundle or preserve both copyright and license |
| Julia project license | Missing | Rights holder selects and adds root `LICENSE` |
| Thesis-era MATLAB author code | Permission not recorded | Confirm ownership and license scope |
| TGA raw/derived data | Permission and data license not recorded | Confirm institutional rights; select exact release paths and data license |
| Instrument metadata | Low privacy risk, not cleared | Confirm laboratory/operator metadata policy |
| RAR archives | Contents not release-audited | Exclude or inventory separately |
| Citation metadata | Missing | Add `CITATION.cff` after author/title/identifier confirmation |

## Recommended release boundary

Until rights are confirmed, prepare a **curated code-only candidate bundle** containing the
Julia package, configuration schema/defaults, tests, user/developer documentation, manifests,
and third-party notices. Exclude `code/` and keep data-derived audit reports, figures, and
scientific summaries in the local validation bundle until their publication scope is cleared.
Do not publish or tag the code-only bundle until the project license is selected. If the data
are later cleared, publish a separately scoped data deposit with its own license and citation,
then point `config/datasets.toml` or a release-specific catalog to that deposit.

The public-facing current status and maintainer decisions are listed in
[`../DATA_SHARING.md`](../DATA_SHARING.md). This review is technical due diligence, not legal
advice.
