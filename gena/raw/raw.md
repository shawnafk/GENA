# Frozen processing entry point

For production raw conversion use `gena.raw.pipeline.process_raw` or `python -m gena.raw.frozen_raw`. It writes a new immutable package with packet inventory, L0/L1/L1A/L2 products, `Q_*` quality masks, `manifest.json`, and rejection diagnostics. The fixed level meanings and v1 thresholds are in `docs/raw_audit_20261027/processing_rules_v1.md`.

The older `decode()`/`raw_2_pqt` path remains for historical product reproduction. It uses append-style daily writes and its L1A→L1B repair branch is disabled by `conversion.l1a_2_l1b`; it is not the v1 production path.
