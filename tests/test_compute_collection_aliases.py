from lib.EnsoComputeMetricsLib import _dataset_variable_entry
from lib.EnsoUvcdatToolsLib import MyDerive
from lib.XarrayCompat import create_axis, create_variable
import lib.EnsoComputeMetricsLib as compute_lib


def test_ssh_dataset_entry_can_use_zos_alias():
    zos_entry = {
        "path + filename": "/tmp/zos.nc",
        "varname": "zos",
    }
    dataset = {"zos": zos_entry}

    assert _dataset_variable_entry(dataset, "ssh") is zos_entry


def test_canonical_dataset_entry_takes_precedence_over_alias():
    ssh_entry = {
        "path + filename": "/tmp/ssh.nc",
        "varname": "ssh",
    }
    zos_entry = {
        "path + filename": "/tmp/zos.nc",
        "varname": "zos",
    }
    dataset = {"ssh": ssh_entry, "zos": zos_entry}

    assert _dataset_variable_entry(dataset, "ssh") is ssh_entry


def test_sst_dataset_entry_can_use_ts_alias():
    ts_entry = {
        "path + filename": "/tmp/ts.nc",
        "varname": "ts",
    }
    dataset = {"ts": ts_entry}

    assert _dataset_variable_entry(dataset, "sst") is ts_entry


def test_thf_dataset_entry_can_use_precomputed_netflux_alias():
    netflux_entry = {
        "path + filename": "/tmp/netflux.nc",
        "varname": "netflux",
    }
    dataset = {"netflux": netflux_entry}

    assert _dataset_variable_entry(dataset, "thf") is netflux_entry


def test_myderive_accepts_precomputed_thf_alias():
    time = create_axis([0, 1], id="time", units="days since 2000-01-01", axis_type="T")
    netflux = create_variable([1.0, 2.0], axes=[time], id="netflux")

    out, err = MyDerive("CMIP", "thf", {"netflux": netflux})

    assert err is None
    assert out is netflux


def test_compute_collection_passes_zos_alias_downstream(monkeypatch):
    calls = []

    def fake_collection(_):
        return {
            "long_name": "test",
            "description": "test",
            "metrics_list": {
                "EnsoFbTauxSsh": {
                    "variables": ["taux", "ssh"],
                    "regions": {"taux": "nino4", "ssh": "nino3"},
                }
            },
        }

    def fake_compute_metric(*args, **kwargs):
        calls.append((args, kwargs))
        return (
            {"metric": {}, "diagnostic": {"MODEL": {"value": 1.0, "value_error": None}}},
            {"metric": {"name": "EnsoFbTauxSsh"}, "diagnostic": {"MODEL": {"name": "MODEL"}}},
            {},
            {},
        )

    monkeypatch.setattr(compute_lib, "defCollection", fake_collection)
    monkeypatch.setattr(compute_lib, "ComputeMetric", fake_compute_metric)

    dict_datasets = {
        "model": {
            "MODEL": {
                "taux": {
                    "path + filename": "/tmp/tauu.nc",
                    "varname": "tauu",
                    "path + filename_area": None,
                    "areaname": None,
                    "path + filename_landmask": None,
                    "landmaskname": None,
                },
                "zos": {
                    "path + filename": "/tmp/zos.nc",
                    "varname": "zos",
                    "path + filename_area": None,
                    "areaname": None,
                    "path + filename_landmask": None,
                    "landmaskname": None,
                },
            }
        },
        "observations": {},
    }

    values, _ = compute_lib.ComputeCollection("test", dict_datasets, "MODEL")

    assert "EnsoFbTauxSsh" in values["value"]
    assert calls
    args, kwargs = calls[0]
    assert args[1] == "EnsoFbTauxSsh"
    assert args[3] == "/tmp/tauu.nc"
    assert args[4] == "tauu"
    assert kwargs["modelFile2"] == "/tmp/zos.nc"
    assert kwargs["modelVarName2"] == "zos"


def test_compute_collection_passes_ts_alias_downstream_for_model_and_obs(monkeypatch):
    calls = []

    def fake_collection(_):
        return {
            "long_name": "test",
            "description": "test",
            "metrics_list": {
                "BiasSstLonRmse": {
                    "variables": ["sst"],
                    "regions": {"sst": "equatorial_pacific"},
                }
            },
        }

    def fake_compute_metric(*args, **kwargs):
        calls.append((args, kwargs))
        return (
            {"metric": {"OBS": {"value": 1.0, "value_error": None}}, "diagnostic": {"MODEL": {}, "OBS": {}}},
            {
                "metric": {"name": "BiasSstLonRmse", "units": "degC"},
                "diagnostic": {
                    "MODEL": {"name": "MODEL", "nyears": 1, "time_period": "model"},
                    "OBS": {"name": "OBS", "nyears": 1, "time_period": "obs"},
                },
            },
            {},
            {},
        )

    monkeypatch.setattr(compute_lib, "defCollection", fake_collection)
    monkeypatch.setattr(compute_lib, "ComputeMetric", fake_compute_metric)

    dict_datasets = {
        "model": {
            "MODEL": {
                "ts": {
                    "path + filename": "/tmp/model_ts.nc",
                    "varname": "ts",
                    "path + filename_area": None,
                    "areaname": None,
                    "path + filename_landmask": None,
                    "landmaskname": None,
                }
            }
        },
        "observations": {
            "OBS": {
                "tos": {
                    "path + filename": "/tmp/obs_tos.nc",
                    "varname": "tos",
                    "path + filename_area": None,
                    "areaname": None,
                    "path + filename_landmask": None,
                    "landmaskname": None,
                }
            }
        },
    }

    compute_lib.ComputeCollection("test", dict_datasets, "MODEL")

    assert calls
    args, _ = calls[0]
    assert args[3] == "/tmp/model_ts.nc"
    assert args[4] == "ts"
    assert args[5] == ["OBS"]
    assert args[6] == ["/tmp/obs_tos.nc"]
    assert args[7] == ["tos"]
