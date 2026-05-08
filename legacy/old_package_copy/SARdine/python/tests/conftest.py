"""Global test configuration for SARdine Python tests.

Critical: SARdine reads several environment flags at import time.

Therefore we set strict defaults as soon as pytest imports this conftest module
(during collection, before test modules import `sardine`).

We use `setdefault` so callers can override intentionally.
"""

import os


os.environ.setdefault("SARDINE_STRICT_CLAMP", "1")
os.environ.setdefault("SARDINE_SERDE_ONLY", "1")
os.environ.setdefault("SARDINE_REQUIRE_SUBSWATHS", "1")

