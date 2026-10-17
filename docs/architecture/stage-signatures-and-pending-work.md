# Persistent stage signatures and pending work

Status: **Proposed; separate from the bulk-planning implementation.**

The current optimization batches the existing database reads and shares missing
destination-directory resolution during Meeko preflight. It changes no stored
fingerprint formats, integrity policy, failure retry policy or schema. This
proposal is a subsequent architectural change and must be reviewed separately.

## Objective and boundary

Resume should discover pending scientific work from durable records rather than
reconstructing every eligibility decision. Completed results remain reusable
only for their exact effective scientific inputs and verified outputs. Existing
historical docking results must remain available after parameters are reverted.

A durable signature proves recorded provenance, **not external file integrity**.
Under the existing strict policy, each relevant artifact must still be checked
against its stored SHA-256 on resume. A clean stage signature cannot justify
skipping that verification. This design reduces database planning and repeated
scientific fingerprint computation; it does not promise work proportional only
to the pending count while also detecting arbitrary corruption of every file.

Size/mtime-only verification would be a separate, explicitly selected weaker
policy, outside this proposal's default behavior.

## Identities and signatures

Use a versioned canonical signature envelope containing the stage's scientific
contract version, normalized effective parameters, upstream artifact identities
and hashes, relevant tool identities, and the scientific entity identity.
Separate scheduling and display metadata from scientific identity.

| Stage | Scope and inputs | Exclusions / distinctions |
|---|---|---|
| Screening | Parent ID and immutable source SMILES; enabled rules/policy; RDKit and screening contract version | Notes and annotations do not invalidate chemistry. |
| MolScrub success | Parent; effective chemical settings; source SMILES; generation toolchain identity | Worker count and timeout do not change a successful chemical result. Preserve existing compatibility rules when migrating. |
| MolScrub failure eligibility | Same parent and attempted chemical settings/toolchain; timeout and failure policy | Increasing timeout may retry failures without discarding successful states. Explicit retry bypasses matching known failures. |
| Meeko | State/conformer; verified SDF hash; effective preparation settings; Meeko and relevant toolchain identity | Worker count is excluded. Adding this signature strengthens the present output-status/hash guard and requires its own compatibility policy. |
| Docking | State/conformer and receptor profile; verified ligand/receptor hashes; normalized effective backend parameters; engine/kernel hashes | Receptor name/reference-ligand annotations are excluded. CPU and GPU have distinct parameter sets. Preserve historical matching and numeric normalization. |
| Export | Source docking run/hash; requested export mode and pose limit; export toolchain; receptor/reference export inputs | Output presence alone cannot prove all requested poses and bundle files are complete. |

Define explicit contract versions rather than incorporating an entire application
commit into every signature. A documentation or UI change must not invalidate
all scientific results. Record application commit/version separately for audit.
Specify toolchain-version coverage with tests before changing today's semantics;
do not silently relabel old records as if missing provenance had been captured.

## Proposed additive records

Existing scientific tables and attempt history remain authoritative evidence.
New tables are indexed projections for planning, not replacements for history.
The following is a logical schema, not a migration ready to execute.

`work_targets`

- Stable scope: `(stage, entity_type, entity_id, receptor_profile_id)`; use an
  empty receptor ID for stages without receptors, avoiding nullable-unique-key
  behavior. Conformer identity must be explicit where multiple conformers are
  supported; the current Meeko single-output restriction remains until redesigned.
- `desired_signature`, `signature_contract_version`, `plan_generation`.
- `status`: `pending`, `running`, `completed`, `failed`, `interrupted`, or
  `blocked` (for an unavailable upstream dependency, not an ambiguous failure).
- `result_id`, `current_attempt_id`, `failure_signature`, failure details,
  timestamps and invalidation reason.
- Unique scope key; index `(stage, receptor_profile_id, status, entity_id)` for
  bounded pending reads. Failed rows are selected according to stage retry policy.

`work_results`

- Immutable result identity, scope, scientific signature and contract version.
- References to existing stage/docking result records and a complete output
  manifest of artifact IDs/hashes and expected result cardinality.
- Index `(stage, entity_type, entity_id, receptor_profile_id, signature)` so a
  reverted configuration can find an older matching result.
- Import provenance distinguishes verified legacy compatibility from native
  records. Never overwrite an old result when a new attempt fails.

`work_dependencies`

- Downstream target and upstream entity/result/artifact identity.
- Reverse indexes on upstream identity for targeted invalidation, including
  dependencies on a desired upstream result that is not available yet.
- A membership revision for changing collections, such as a parent's retained
  states or an enabled receptor set. Removing a state/profile must obsolete its
  dependent targets without deleting historical scientific evidence.

`stage_completions`

- Stage/scope, library-membership revision, settings/tool revisions,
  completion signature, pending/failed/completed counts and plan generation.
- Written only when every required target is terminal under the defined policy.
  Scientific success and "all work considered, with failures" are distinct.
- A valid completion record allows reuse of the recorded plan; strict artifact
  verification remains necessary before treating results as reusable.

`project_revisions` / input import record

- Separate library membership, source/input metadata, chemical settings,
  receptor settings and toolchain revisions. Notes-only changes leave chemical
  revisions unchanged.
- Last successfully imported CSV content hash and import contract version. Hash
  the file before deciding it is unchanged; timestamps alone are insufficient.
- Update the successful-import record only in the transaction that commits the
  validated cumulative import. Invalid CSV input must not advance it.
- Track all supported import entry points; a cached file signature must never
  hide a library change made through another import path.

## Planning, execution and invalidation

1. Acquire the existing project execution lock. Capture desired revisions and a
   new plan generation. Synchronize changed inputs and recover interrupted work.
2. Apply revision/dependency invalidations. New parents dirty their own branch;
   changed receptor/box settings dirty only relevant receptor/state targets;
   annotations or worker-count changes do not dirty scientific results.
3. Verify relevant recorded artifacts against full hashes. Missing or changed
   content invalidates only affected results/dependents. A cached verification
   is limited to the run and invalidated on application writes; do not silently
   make it persistent or trust a database signature as a file check.
4. Read desired pending/eligible-failed targets in indexed bounded pages. Reuse
   compatible verified `work_results` before scheduling a new attempt.
5. Claim each attempt with a compare-and-swap from the expected status and plan
   generation. Workers compute; the coordinating thread owns database writes.
6. Finalize outputs and result evidence, then atomically update the target and
   dependent invalidations. Publish aggregate reuse counts and pending progress.
7. Commit a stage completion projection only for the captured revision set and
   plan generation. A changed desired revision prevents stale completion.

Persistent work is a set of identities, not a last-completed-ID cursor.
Parallel jobs finish out of order, failures create holes, and one parent may
produce multiple states for multiple receptors. No prefix assumption is valid.

Stage retry behavior initially remains compatible: MolScrub skips matching known
failures unless explicitly retried or its failure signature changes; Meeko and
docking currently retry failed work. A later explicit retry UI may revise that
policy, but the queue migration must not silently change it.

## Atomicity, cancellation and recovery

Write artifact content to temporary files and finalize it before committing its
manifest and completed result in one short database transaction. Repository
methods need transaction-aware variants so result/manifest/target changes do not
commit independently. Never hold a transaction open across a subprocess run.

A crash before file finalization leaves disposable temporary files. A crash
after finalization but before the completion transaction leaves untrusted orphan
files; their presence must not mark a target complete. Interrupted attempts are
recovered as pending/eligible under their desired signature; completed results
are not rewritten. Preserve attempt lineage rather than replacing it on retry.

Keep the existing OS project lock as the single-coordinator boundary. Persistent
status is not a replacement for that lock. If cross-process worker claims are
introduced later, add explicit leases, owner identity and expiry/heartbeat
semantics; a timestamp-based lease must not steal a live worker's files.

On Stop, cease claims, drain or cancel active processes with the current safety
rules, and commit each finished result. Unstarted targets remain pending. A late
completion is accepted only for its claimed attempt/signature/generation; it may
be retained as historical evidence but cannot overwrite a newer desired target.

## Compatibility and rollout

1. Add the projection schema and write records alongside existing execution.
   Continue using today's eligibility code to make decisions.
2. Backfill in bounded chunks using existing fingerprint/compatibility checks
   and strict artifact verification. Keep historical docking ordering/current
   selection. Unknown provenance remains unproven, with an explicit compatibility
   classification; do not manufacture a complete modern signature.
3. Run the new planner in comparison mode against existing decisions, recording
   discrepancies without executing a second set of jobs. Test on copied projects.
4. Enable indexed pending planning after equivalence and migration tests pass.
   Retain a fallback to the current batch planner for unmigrated/unsupported
   records, and make backfill resumable and versioned.
5. Review Meeko preparation and export signature strengthening as explicit
   behavior changes, including how previously accepted artifacts are treated.

Historical tables are never destructively rewritten to make the new projections
look complete. Rebuilding projections must be recoverable from original evidence.

## Required validation and observability

- Compare reuse identities and selected historical run IDs with the existing
  batch planner across mixed statuses, current/historical results, legacy
  compatibility, settings reversion, inactive data and missing/corrupt files.
- Test adding one ligand, notes-only changes, timeout/worker changes, changed
  SDF coordinates, changed tools/kernels, receptor enable/disable and removal,
  multiple states and interruption at every file/transaction boundary.
- Test noncontiguous completions, repeated failures, explicit retries, stale
  generations, cancellation during planning and restart during backfill.
- Require bounded query/connection counts per planning page and zero scientific
  tool calls on a fully reusable resume. Avoid machine-specific timing assertions.
- Record CSV sync, metadata planning, artifact verification, collision checks,
  scientific execution and UI time separately; report pending/reused/known-failed
  counts without delivering a GUI/log event for every already-completed item.

Design decisions to resolve before implementation: precise toolchain
coverage per signature contract, supported conformer multiplicity, legacy
compatibility treatment for Meeko/export, and any future verification or retry
policy change. None of these changes are introduced by the present optimization.
