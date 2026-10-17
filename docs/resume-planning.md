# Bulk resume planning

The pipeline retains its existing fingerprint, status, compatibility and full
SHA-256 checks. This optimization changes how their metadata is read, not the
criteria for accepting a previous result. No database migration is needed.

- Screening joins each parent's active result into its initial library query.
- MolScrub reads latest per-parent history and existing states in batches of
  256 parents. Latest means greatest rowid, including failed/interrupted records;
  it does not mean latest successful attempt. Each SDF is still checked exactly
  where the original predicate requires it, including legacy generation hashes.
- Meeko joins conformer status into the initial read and batches active output
  artifact lookups. It announces preflight before filesystem work and checks
  cancellation while planning. The pending count is reported at planning completion.
- Docking batches prepared-artifact metadata and scored completed history per
  receptor and 256 state occurrences. Candidate order remains current first,
  then newest rowid. All original fingerprint and raw-output checks still run.
  Duplicate state occurrences use live lookups after the first occurrence so
  they can observe results just completed in the same run.

Metadata batches close their connections before hashing files, yielding work or
running tools. Actual completion/failure transactions and project locking remain
unchanged. These are bounded metadata reads rather than a persistent work queue;
scientific execution and integrity verification still occur as before.

## Meeko destinations

Every output and log destination is still checked for collisions before any
worker writes, including destinations belonging to reusable outputs. Existing
paths are resolved with `Path.resolve()` so symlinks, junctions and Windows
aliases retain their physical interpretation.

Missing outputs share resolution of missing ancestors in a preflight-only cache.
The nearest existing ancestor is resolved normally; missing suffixes are appended
without repeatedly asking the filesystem to resolve the same directories.
Windows lexical normalization is performed before physical resolution, matching
Windows `realpath` behavior for `..`; POSIX component ordering is preserved.
Other filesystem errors fall back to normal pathlib behavior. The cache is
discarded before workers write and is never used for artifact-integrity decisions.

## Validation

`tests/test_bulk_planning.py` compares decisions with reference queries and
predicates from the previous pipeline. Cases include completed/failed/stale and
interrupted history, missing and inactive state/artifact metadata, legacy
generation fingerprints, missing/corrupt files, retries and cancellation.
Tests also bound cache-read connection counts by batch count and require every
previously required SDF hash check to remain present.

Docking equivalence tests in `tests/test_idempotency.py` compare exact selected
run records from the original per-item lookup and the batch metadata reader,
including historical results, legacy compatibility and damaged pose metadata.
The existing parallel/collision, multi-receptor, invalidation and stop/resume
regressions remain applicable. Path-resolution tests compare pathlib results,
physical alias collisions and reuse of missing ancestors; real symlink coverage
is skipped on hosts that do not allow creating symlinks.

A local planning experiment used a snapshot of the 5,043-parent project, with
no scientific jobs executed and no student project writes:

| Check | Previous planning replay, median of 3 passes | Implemented planning, one pass |
|---|---:|---:|
| Screening | 31.40 s | 0.23 s |
| MolScrub | 36.31 s | 4.50 s |
| Meeko | 28.05 s | 4.61 s |

Screening reused all 5,043 results. MolScrub reused 1,425 parents and skipped
629 known failures while verifying the same 2,045 SDFs. Meeko reused 48 outputs
and selected the same 1,997 pending conformers. Previous timings are from the
investigation's simplified read-only replay, not a full application profiler;
the implemented measurements use the real stage methods for screening/Meeko
and the real batch reader plus unchanged predicates for MolScrub. Database
snapshots and suppressed scientific dispatch prevent modifying project data.
Filesystem caching and machine conditions affect these timings; they are not
end-to-end performance guarantees or test thresholds.

Persistent completion signatures, invalidation projections and indexed pending
work are a separate proposal in
[stage-signatures-and-pending-work.md](architecture/stage-signatures-and-pending-work.md).
