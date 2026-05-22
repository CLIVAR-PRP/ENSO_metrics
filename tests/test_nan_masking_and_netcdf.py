import numpy as np
import numpy.ma as ma
import xarray as xr

from lib.EnsoUvcdatToolsLib import (
    AverageTemporal,
    AverageZonal,
    Correlation,
    GENUTILlinearregression,
    SaveNetcdf,
)
from lib.XarrayCompat import create_axis, create_rect_grid, create_variable


def _axis(values, axis, name, units):
    return create_axis(
        values,
        id=name,
        units=units,
        attributes={"axis": axis},
    )


def test_zonal_average_masks_when_half_or_more_inputs_are_nan():
    lat = _axis([-0.5, 0.5], "Y", "lat", "degrees_north")
    lon = _axis([10.0, 11.0], "X", "lon", "degrees_east")
    data = ma.array(
        [[1.0, 2.0], [3.0, 4.0]],
        mask=[[True, False], [False, False]],
    )
    var = create_variable(data, axes=[lat, lon], grid=create_rect_grid(lat, lon), id="sst")

    averaged, err = AverageZonal(var)

    assert err is None
    assert ma.getmaskarray(averaged._data).tolist() == [True, False]
    assert np.isclose(float(averaged._data[1]), 3.5)


def test_temporal_average_masks_when_half_or_more_inputs_are_nan():
    time = _axis([0, 1, 2, 3], "T", "time", "days since 2000-01-01")
    var = create_variable(
        ma.array([1.0, 2.0, 3.0, 4.0], mask=[True, False, True, False]),
        axes=[time],
        id="nino",
    )

    averaged, err = AverageTemporal(var)

    assert err is None
    assert ma.is_masked(averaged._data)


def test_same_series_correlation_and_regression_are_exact_one():
    time = _axis([0, 1, 2, 3], "T", "time", "days since 2000-01-01")
    var = create_variable([1.0, 2.0, 4.0, 8.0], axes=[time], id="nino")

    corr = Correlation(var, var, axis=0)
    slope, stderr = GENUTILlinearregression(var, x=var, error=1, nointercept=1)

    assert float(corr) == 1.0
    assert float(slope[0, 0]) == 1.0
    assert float(stderr[0, 0]) == 0.0


def test_save_netcdf_preserves_variables_with_different_same_named_dims(tmp_path):
    time_a = _axis([0, 1], "T", "time", "days since 2000-01-01")
    lat_a = _axis([-0.5, 0.5], "Y", "lat", "degrees_north")
    lon_a = _axis([10.0, 11.0], "X", "lon", "degrees_east")
    var_a = create_variable(
        np.ones((2, 2, 2)),
        axes=[time_a, lat_a, lon_a],
        grid=create_rect_grid(lat_a, lon_a),
        id="a",
    )

    time_b = _axis([10, 11, 12], "T", "time", "days since 2000-01-01")
    lat_b = _axis([1.5], "Y", "lat", "degrees_north")
    lon_b = _axis([20.0, 21.0, 22.0], "X", "lon", "degrees_east")
    var_b = create_variable(
        np.ones((3, 1, 3)) * 2.0,
        axes=[time_b, lat_b, lon_b],
        grid=create_rect_grid(lat_b, lon_b),
        id="b",
    )

    path = tmp_path / "multi.nc"
    SaveNetcdf(str(path), var1=var_a, var1_name="a", var2=var_b, var2_name="b")

    with xr.open_dataset(path, decode_times=False) as ds:
        assert ds["a"].shape == (2, 2, 2)
        assert ds["b"].shape == (3, 1, 3)
        assert ds["a"].dims == ("time", "lat", "lon")
        assert ds["b"].dims == ("time_b", "lat_b", "lon_b")


def test_append_netcdf_preserves_variables_with_different_same_named_dims(tmp_path):
    time_a = _axis([0, 1], "T", "time", "days since 2000-01-01")
    lat_a = _axis([-0.5, 0.5], "Y", "lat", "degrees_north")
    lon_a = _axis([10.0, 11.0], "X", "lon", "degrees_east")
    var_a = create_variable(
        np.ones((2, 2, 2)),
        axes=[time_a, lat_a, lon_a],
        grid=create_rect_grid(lat_a, lon_a),
        id="a",
    )

    time_b = _axis([10, 11, 12], "T", "time", "days since 2000-01-01")
    lat_b = _axis([1.5], "Y", "lat", "degrees_north")
    lon_b = _axis([20.0, 21.0, 22.0], "X", "lon", "degrees_east")
    var_b = create_variable(
        np.ones((3, 1, 3)) * 2.0,
        axes=[time_b, lat_b, lon_b],
        grid=create_rect_grid(lat_b, lon_b),
        id="b",
    )

    path = tmp_path / "append.nc"
    SaveNetcdf(str(path), var1=var_a, var1_name="a")
    SaveNetcdf(str(path), var1=var_b, var1_name="b")

    with xr.open_dataset(path, decode_times=False) as ds:
        assert ds["a"].shape == (2, 2, 2)
        assert ds["b"].shape == (3, 1, 3)
        assert ds["a"].dims == ("time", "lat", "lon")
        assert ds["b"].dims == ("time_b", "lat_b", "lon_b")
