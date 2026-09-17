import mmap
import pathlib
import sys
from collections import Counter

import pytest

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "tools"))

from compare_breakpoint_loading import decode_block, read_header, read_index
from ibd_to_binary import (
    BLOCK_CHECKPOINT,
    BinaryWriter,
    Indices,
    read_segments,
    sweep_and_write,
)


@pytest.mark.parametrize("checkpoint_interval", [1, 2, 3, 1000])
def test_zero_length_segments_cancel_in_checkpoints_and_deltas(tmp_path, checkpoint_interval):
    source = tmp_path / "segments.txt"
    source.write_text(
        "a 1 c 1 1 100 100 4\n"  # Zero-length segment at the first checkpoint.
        "a 1 b 1 1 200 300 4\n"
        "a 2 b 2 1 200 400 4\n"  # Preserve overlapping segment multiplicity.
        "a 1 b 1 1 250 250 4\n"  # Same pair as active nonzero-length segments.
        "a 1 c 1 1 250 250 4\n"  # Pair without any active segments.
        "a 1 b 1 1 300 400 4\n"  # A different segment starts as one ends.
        "a 1 c 1 1 350 350 4\n"
        "a 1 c 1 1 500 500 4\n"  # Preserve the final empty breakpoint.
    )
    segments, events, names = read_segments(str(source), Indices("hap-ibd"), 3.0)
    path = tmp_path / "output.ibdf"
    writer = BinaryWriter(str(path))
    max_active, _, _ = sweep_and_write(segments, events, writer, checkpoint_interval)
    writer.finalize(3, checkpoint_interval)

    pair = (names.index("a"), names.index("b"))
    expected = {
        100: Counter(),
        200: Counter({pair: 2}),
        250: Counter({pair: 2}),
        300: Counter({pair: 2}),
        350: Counter({pair: 2}),
        400: Counter(),
        500: Counter(),
    }
    assert max_active == 2
    with path.open("rb") as handle, mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as data:
        index = read_index(data, read_header(data))
        assert [pos for pos, _, _ in index] == list(expected)
        active = Counter()
        for pos, offset, _ in index:
            block_type, block = decode_block(data, offset)
            if block_type == BLOCK_CHECKPOINT:
                _cm, p1, p2 = block
                active = Counter(zip(p1, p2))
            else:
                _add_cm, add_p1, add_p2, _del_cm, del_p1, del_p2 = block
                # Match IBDMap's deletion-before-addition order and reject
                # deletions for inactive pairs, even if an addition follows.
                for pair in zip(del_p1, del_p2):
                    assert active[pair] > 0
                    active[pair] -= 1
                    if active[pair] == 0:
                        del active[pair]
                active.update(zip(add_p1, add_p2))
            assert active == expected[pos]
