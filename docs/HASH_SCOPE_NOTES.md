# Hash scope notes (2026-10-08 audit)

This sidecar note labels what the hash and checksum fields in the files below actually cover. The data files themselves are unchanged. It is a documentation clarification, not a license verdict and not a source re-admission.

Audited commit: `8df8fb86812cfc30f10459fa06bf7fa82945d3ec`

## `data_cache/external/geo/clinical_pull_manifest.json`
- Lines (verified against the audited commit): 4-6, 10-12, 16-18, 22-24, 28-30 (15 lines)
- Scope: Each sha256 here is the hash of the locally generated clinical CSV, not of the source GEO file it was derived from. It is not source-payload integrity.

