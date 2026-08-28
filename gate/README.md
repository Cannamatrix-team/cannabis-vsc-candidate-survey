# Final validation gate

The manuscript-active gate has two inputs:

1. Motif support passes operationally when a protein is `confirmed` or `partial`, or when it has the explicit singleton passthrough.
2. Expected-domain support passes when the configured architecture is `complete`. A protein with no historical domain record and no configured post-resolution rule passes by default.

The historical physicochemical-stability and contamination switches were off, and selection context was annotation-only. They cannot change the deposited catalog and are not included in the public gate.
