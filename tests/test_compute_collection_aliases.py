from lib.EnsoComputeMetricsLib import _dataset_variable_entry
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
