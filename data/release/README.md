# Release assets

Large files are distributed through the GitHub release rather than committed to Git history.

Release asset for `v1.0.0`:

- `cs-no-10x-filtered-matrices-v1.0.0.tar.gz`

Verify downloads with:

```bash
cd data/release
shasum -a 256 -c SHA256SUMS
```

The manifest records the size and SHA-256 digest of every source file included in the assets.
