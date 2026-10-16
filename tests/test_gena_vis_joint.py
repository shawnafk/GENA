import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from plots.basic import gena as vis


def test_joint_uses_bounded_integer_colorbar_ticks_for_sparse_counts():
    values = np.arange(20, dtype=float)
    frame = pd.DataFrame(
        {
            "TX_1s": values,
            "TX_2s": values + 1,
            "TX_1e": values + 2,
            "TX_2e": values + 3,
            "TY_1s": values + 4,
            "TY_2s": values + 5,
            "TY_1e": values + 6,
            "TY_2e": values + 7,
        }
    )

    figure, axes = vis.joint(frame)
    figure.canvas.draw()

    assert len(axes) == 4
    colorbar_axes = [axis for axis in figure.axes if axis not in axes]
    assert len(colorbar_axes) == 4
    for axis in colorbar_axes:
        assert axis.get_xticks().tolist() == [1.0]
        labels = [label.get_text() for label in axis.get_xticklabels()]
        assert labels == ["1"]
