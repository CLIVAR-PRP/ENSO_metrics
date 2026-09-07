import numpy as np
import numpy.ma as ma
import pytest

from lib.EnsoUvcdatToolsLib import MV2concatenate, MV2zeros, get_year_by_year
from lib.XarrayCompat import create_axis, create_rect_grid, create_variable


def _time_axis(start_year, start_month, n_months):
    # 360_day calendar keeps every month exactly 30 days long
    days = (start_month - 1) * 30 + 15 + np.arange(n_months) * 30
    return create_axis(
        days,
        id="time",
        units="days since %d-01-01" % start_year,
        attributes={"axis": "T", "calendar": "360_day"},
    )


def _monthly_series(start_year, start_month, n_months):
    time = _time_axis(start_year, start_month, n_months)
    return create_variable(ma.arange(n_months, dtype=float) + 1, axes=[time], id="sst")


def _monthly_map(start_year, start_month, n_months):
    time = _time_axis(start_year, start_month, n_months)
    lat = create_axis([-0.5, 0.5], id="lat", units="degrees_north", attributes={"axis": "Y"})
    lon = create_axis([10.0, 11.0, 12.0], id="lon", units="degrees_east", attributes={"axis": "X"})
    data = np.arange(n_months, dtype=float)[:, None, None] + 1 + np.zeros((n_months, 2, 3))
    mask = np.zeros(data.shape, dtype=bool)
    mask[:, 0, 0] = True
    return create_variable(
        ma.array(data, mask=mask),
        axes=[time, lat, lon],
        grid=create_rect_grid(lat, lon),
        id="sst",
        attributes={"units": "K"},
    )


@pytest.mark.parametrize("pad_first", [False, True])
def test_concatenate_with_axisless_padding(pad_first):
    var = _monthly_series(2000, 1, 4)
    pad = MV2zeros([1])
    seq = (pad, var) if pad_first else (var, pad)

    out = MV2concatenate(seq)

    assert out.shape == (5,)
    assert out.id == "sst"
    assert out.getAxis(0) is None


def test_concatenate_with_axisless_padding_keeps_spatial_metadata():
    var = _monthly_map(2000, 1, 4)
    pad = MV2zeros([1, 2, 3])

    out = MV2concatenate((pad, var))

    assert out.shape == (5, 2, 3)
    assert out.id == "sst"
    assert out.units == "K"
    assert out.getAxis(0) is None
    assert [ax.id for ax in out.getAxisList()[1:]] == ["lat", "lon"]
    assert out.getGrid() is not None


def test_concatenate_keeps_axis_when_all_inputs_have_one():
    a = _monthly_series(2000, 1, 2)
    b = _monthly_series(2000, 3, 2)

    out = MV2concatenate((a, b))

    assert out.shape == (4,)
    assert out.getAxis(0).id == "time"
    assert len(out.getAxis(0)) == 4


def test_year_by_year_series_ending_mid_year():
    var = _monthly_series(2000, 1, 30)

    out = get_year_by_year(var, frequency="monthly")

    assert out.shape == (3, 12)
    assert [ax.id for ax in out.getAxisList()] == ["years", "months"]
    mask = ma.getmaskarray(out._data)
    assert not mask[:2].any()
    assert not mask[2, :6].any()
    assert mask[2, 6:].all()
    assert np.allclose(out._data[0], np.arange(1, 13))
    assert np.allclose(out._data[2, :6], np.arange(25, 31))


def test_year_by_year_series_starting_mid_year():
    var = _monthly_series(2000, 7, 30)

    out = get_year_by_year(var, frequency="monthly")

    assert out.shape == (3, 12)
    mask = ma.getmaskarray(out._data)
    assert mask[0, :6].all()
    assert not mask[0, 6:].any()
    assert not mask[1:].any()
    assert np.allclose(out._data[0, 6:], np.arange(1, 7))
    assert np.allclose(out._data[2], np.arange(19, 31))


def test_year_by_year_map_ending_mid_year():
    var = _monthly_map(2000, 1, 30)

    out = get_year_by_year(var, frequency="monthly")

    assert out.shape == (3, 12, 2, 3)
    assert [ax.id for ax in out.getAxisList()] == ["years", "months", "lat", "lon"]
    assert out.getGrid() is not None
    mask = ma.getmaskarray(out._data)
    assert mask[2, 6:].all()
    assert mask[:, :, 0, 0].all()
    assert not mask[:2, :, 1, :].any()
    assert np.allclose(out._data[2, :6, 1, 1], np.arange(25, 31))


def test_year_by_year_full_years_unchanged():
    var = _monthly_series(2000, 1, 24)

    out = get_year_by_year(var, frequency="monthly")

    assert out.shape == (2, 12)
    assert not ma.getmaskarray(out._data).any()
