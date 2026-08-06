"""Read processing and per-cell base counting."""

from .pileup import ShardResult, plan_shards, scan_shard
from .processors import build_tasks, process_shards
from .readers import BAMReader

__all__ = [
    "BAMReader",
    "ShardResult",
    "build_tasks",
    "plan_shards",
    "process_shards",
    "scan_shard",
]
