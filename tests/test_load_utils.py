import numpy as np
import pandas as pd

import common
import common.frame as common_frame
from load.loaders.gena import _clip_variable_dict


def test_clip_variable_dict_is_private_to_load_and_clips_inclusively():
    raw = {
        "spectrum": {
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": np.arange(8).reshape(4, 2),
            "v": np.array([10.0, 20.0]),
            "v1": np.arange(12).reshape(4, 3),
            "units": "counts",
        },
        "missing_y": {"x": [0.0, 1.0]},
    }

    result = _clip_variable_dict(raw, [1.0, 2.0])

    assert list(result) == ["spectrum"]
    assert result["spectrum"]["x"].equals(pd.to_datetime([1.0, 2.0], unit="s"))
    np.testing.assert_array_equal(result["spectrum"]["y"], [[2, 3], [4, 5]])
    np.testing.assert_array_equal(result["spectrum"]["v"], [10.0, 20.0])
    np.testing.assert_array_equal(result["spectrum"]["v1"], [[3, 4, 5], [6, 7, 8]])
    assert result["spectrum"]["units"] == "counts"
    assert len(raw["spectrum"]["x"]) == 4
    assert "clip_variable_dict" not in common_frame.__all__
    assert not hasattr(common, "clip_variable_dict")
