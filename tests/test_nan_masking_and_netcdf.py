import numpy as np
import numpy.ma as ma
import pytest
import xarray as xr

from lib.EnsoUvcdatToolsLib import (
    AverageTemporal,
    AverageZonal,
    CheckUnits,
    Correlation,
    GENUTILlinearregression,
    LinearRegressionTsAgainstMap,
    LinearRegressionTsAgainstTs,
    SaveNetcdf,
    open_file,
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


def test_map_regression_uses_only_jointly_valid_pairs():
    time = _axis([0, 1, 2], "T", "time", "days since 2000-01-01")
    lon = _axis([10.0, 11.0], "X", "lon", "degrees_east")
    x = create_variable([1.0, 2.0, 3.0], axes=[time], id="nino")
    y = create_variable(
        ma.array(
            [[2.0, 3.0], [4.0, 6.0], [6.0, 9.0]],
            mask=[[False, False], [False, False], [True, False]],
        ),
        axes=[time, lon],
        id="sst",
    )

    slope, stderr = LinearRegressionTsAgainstMap(y, x)

    assert np.allclose(slope._data, [2.0, 3.0])
    assert np.allclose(stderr._data, [0.0, 0.0])


def test_ts_regression_aligns_when_index_starts_before_field():
    months = [np.datetime64(f"{year}-{month:02d}-15") for year in range(2002, 2012) for month in range(1, 13)]
    time = _axis(months, "T", "time", "")
    years = _axis([np.datetime64(f"{year}-12-15") for year in range(2000, 2012)], "T", "time", "")
    x_all = np.linspace(-2.0, 2.0, 12)
    y = create_variable(np.repeat(3.0 * x_all[2:], 12), axes=[time], id="sst")
    x = create_variable(x_all, axes=[years], id="enso")

    slope = LinearRegressionTsAgainstTs(y, x, 2, return_stderr=False, frequency="monthly")

    assert slope.shape == (24,)
    assert ma.count(slope._data) == 24
    assert not np.allclose(slope._data, 0.0)


def test_ts_regression_broadcasts_index_against_hovmoeller_field():
    months = [np.datetime64(f"{year}-{month:02d}-15") for year in range(2002, 2012) for month in range(1, 13)]
    time = _axis(months, "T", "time", "")
    lon = _axis([150.0, 170.0], "X", "lon", "degrees_east")
    years = _axis([np.datetime64(f"{year}-12-15") for year in range(2000, 2012)], "T", "time", "")
    x_all = np.linspace(-2.0, 2.0, 12)
    base = np.repeat(x_all[2:], 12)
    y = create_variable(
        np.column_stack([2.0 * base, -1.5 * base]),
        axes=[time, lon],
        id="sst",
    )
    x = create_variable(x_all, axes=[years], id="enso")

    slope = LinearRegressionTsAgainstTs(y, x, 2, return_stderr=False, frequency="monthly")

    assert slope.shape == (24, 2)
    assert ma.count(slope._data) == slope.size
    assert not np.allclose(slope._data, 0.0)


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


def test_open_file_accepts_region_subset_on_curvilinear_lat_lon(tmp_path):
    time = np.array([0.0, 30.0])
    j = np.arange(4)
    i = np.arange(5)
    lat2d = np.array(
        [
            [-8.0, -8.0, -8.0, -8.0, -8.0],
            [-3.0, -3.0, -3.0, -3.0, -3.0],
            [3.0, 3.0, 3.0, 3.0, 3.0],
            [8.0, 8.0, 8.0, 8.0, 8.0],
        ]
    )
    lon2d = np.array(
        [
            [200.0, 220.0, 240.0, 260.0, 280.0],
            [200.0, 220.0, 240.0, 260.0, 280.0],
            [200.0, 220.0, 240.0, 260.0, 280.0],
            [200.0, 220.0, 240.0, 260.0, 280.0],
        ]
    )
    data = np.arange(2 * 4 * 5, dtype=float).reshape(2, 4, 5)
    ds = xr.Dataset(
        {
            "zos": (
                ("time", "j", "i"),
                data,
                {"units": "m", "standard_name": "sea_surface_height_above_geoid"},
            )
        },
        coords={
            "time": (
                "time",
                time,
                {"units": "days since 2000-01-01", "calendar": "standard"},
            ),
            "j": j,
            "i": i,
            "lat": (("j", "i"), lat2d, {"units": "degrees_north"}),
            "lon": (("j", "i"), lon2d, {"units": "degrees_east"}),
        },
    )
    path = tmp_path / "curvilinear_zos.nc"
    ds.to_netcdf(path)

    var = open_file(str(path))("zos", latitude=(-5.0, 5.0), longitude=(210.0, 270.0))

    assert var.shape == (2, 10, 60)
    assert [ax.axis for ax in var.getAxisList()] == ["T", "Y", "X"]
    assert np.allclose(var.getLatitude()[:], np.arange(-4.5, 5.0, 1.0))
    assert np.allclose(var.getLongitude()[:], np.arange(210.5, 270.0, 1.0))
    assert np.isfinite(var.filled(np.nan)).any()


def test_open_file_regrids_mpas_like_unstructured_cells(tmp_path):
    time = np.array([0.0, 30.0])
    lat_deg = np.repeat(np.array([-2.0, 0.0, 2.0]), 3)
    lon_deg = np.tile(np.array([220.0, 222.0, 224.0]), 3)
    data = np.stack(
        [
            lat_deg + lon_deg / 100.0,
            lat_deg + lon_deg / 100.0 + 1.0,
        ]
    )
    ds = xr.Dataset(
        {
            "zos": (
                ("time", "nCells"),
                data,
                {"units": "m", "standard_name": "sea_surface_height_above_geoid"},
            ),
            "latCell": (
                "nCells",
                np.deg2rad(lat_deg),
                {"units": "radians", "standard_name": "latitude"},
            ),
            "lonCell": (
                "nCells",
                np.deg2rad(lon_deg),
                {"units": "radians", "standard_name": "longitude"},
            ),
        },
        coords={
            "time": (
                "time",
                time,
                {"units": "days since 2000-01-01", "calendar": "standard"},
            ),
            "nCells": np.arange(lat_deg.size),
        },
    )
    path = tmp_path / "mpas_like_zos.nc"
    ds.to_netcdf(path)

    var = open_file(str(path))("zos", latitude=(-2.0, 2.0), longitude=(220.0, 224.0))

    assert var.shape == (2, 4, 4)
    assert [ax.axis for ax in var.getAxisList()] == ["T", "Y", "X"]
    assert np.allclose(var.getLatitude()[:], [-1.5, -0.5, 0.5, 1.5])
    assert np.allclose(var.getLongitude()[:], [220.5, 221.5, 222.5, 223.5])
    assert np.isfinite(var.filled(np.nan)).any()


def test_check_units_infers_zos_meters_when_units_missing():
    time = _axis([0], "T", "time", "days since 2000-01-01")
    lat = _axis([0.5], "Y", "lat", "degrees_north")
    lon = _axis([220.5], "X", "lon", "degrees_east")
    var = create_variable(
        np.ones((1, 1, 1)),
        axes=[time, lat, lon],
        grid=create_rect_grid(lat, lon),
        id="zos",
        attributes={"standard_name": "sea_surface_height_above_geoid"},
    )

    with pytest.warns(UserWarning, match="Missing units.*zos.*inferring meters"):
        out, units, keyerror = CheckUnits(
            var,
            "sea surface height",
            "zos",
            "",
            return_tab_only=False,
        )

    assert out is var
    assert units == "m"
    assert keyerror is None


def test_check_units_infers_ssh_map_meters_when_units_missing():
    time = _axis([0], "T", "time", "days since 2000-01-01")
    lat = _axis([0.5], "Y", "lat", "degrees_north")
    lon = _axis([220.5], "X", "lon", "degrees_east")
    var = create_variable(
        np.ones((1, 1, 1)),
        axes=[time, lat, lon],
        grid=create_rect_grid(lat, lon),
        id="ssh_map__ACCESS1-0_r1i1p1",
        attributes={"standard_name": "sea_surface_height_above_geoid"},
    )

    with pytest.warns(UserWarning, match="Missing units.*ssh_map.*inferring meters"):
        _, units, keyerror = CheckUnits(
            var,
            "sea surface height",
            "ssh_map__ACCESS1-0_r1i1p1",
            "",
            return_tab_only=False,
        )

    assert units == "m"
    assert keyerror is None
