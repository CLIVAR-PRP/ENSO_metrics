import numpy as np
import numpy.ma as ma

from lib.EnsoUvcdatToolsLib import MV2concatenate, MV2zeros, get_year_by_year
from lib.XarrayCompat import create_axis, create_variable


def _monthly_series(start_year, start_month, n_months):
    # 360_day calendar keeps every month exactly 30 days long
    days = (start_month - 1) * 30 + 15 + np.arange(n_months) * 30
    time = create_axis(
        days,
        id="time",
        units="days since %d-01-01" % start_year,
        attributes={"axis": "T", "calendar": "360_day"},
    )
    return create_variable(ma.arange(n_months, dtype=float) + 1, axes=[time], id="sst")


def test_concatenate_with_axisless_padding_after():
    var = _monthly_series(2000, 1, 4)
    pad = MV2zeros([1])

    out = MV2concatenate((var, pad))

    assert out.shape == (5,)
    assert out.id == "sst"
    assert out.getAxis(0) is None


def test_concatenate_with_axisless_padding_before():
    var = _monthly_series(2000, 1, 4)
    pad = MV2zeros([1])

    out = MV2concatenate((pad, var))

    assert out.shape == (5,)
    assert out.id == "sst"
    assert out.getAxis(0) is None


def test_concatenate_keeps_axis_when_all_inputs_have_one():
    a = _monthly_series(2000, 1, 2)
    b = _monthly_series(2000, 3, 2)

    out = MV2concatenate((a, b))

    assert out.shape == (4,)
    assert out.getAxis(0).id == "time"
    assert len(out.getAxis(0)) == 4


def test_year_by_year_series_ending_mid_year():
    # Jan 2000 to Jun 2002 covers 30 months over 3 calendar years
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
    # Jul 2000 to Dec 2002 covers 30 months over 3 calendar years
    var = _monthly_series(2000, 7, 30)

    out = get_year_by_year(var, frequency="monthly")

    assert out.shape == (3, 12)
    mask = ma.getmaskarray(out._data)
    assert mask[0, :6].all()
    assert not mask[0, 6:].any()
    assert not mask[1:].any()
    assert np.allclose(out._data[0, 6:], np.arange(1, 7))
    assert np.allclose(out._data[2], np.arange(19, 31))


def test_year_by_year_full_years_unchanged():
    var = _monthly_series(2000, 1, 24)

    out = get_year_by_year(var, frequency="monthly")

    assert out.shape == (2, 12)
    assert not ma.getmaskarray(out._data).any()
