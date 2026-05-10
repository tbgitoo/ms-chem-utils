# tests/conftest.py
import pytest

def _has_rdkit() -> bool:
    try:
        import rdkit  # noqa: F401
        return True
    except Exception:
        return False

requires_rdkit = pytest.mark.skipif(not _has_rdkit(), reason="RDKit not installed")