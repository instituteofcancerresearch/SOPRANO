import pathlib
import tempfile

import pytest

from SOPRANO.sh_utils import setup_conda


@pytest.mark.dependency(name="create_rc")
def test__create_new_condarc():
    with tempfile.TemporaryDirectory() as tmp_dir:
        condarc_path = pathlib.Path(tmp_dir).joinpath(".condarc").as_posix()
        mock_env_location = pathlib.Path.cwd().as_posix()

        expected_lines = [
            "envs_dirs:\n",
            "  - ~/.conda/envs\n",
            f"  - {mock_env_location}\n",
        ]

        setup_conda._create_new_condarc(mock_env_location, condarc_path)
        with open(condarc_path, "r") as f:
            written_lines = f.readlines()

        for expected, written in zip(expected_lines, written_lines):
            assert expected == written, (expected, written)


@pytest.mark.dependency(name="update_rc")
def test__update_condarc():
    mock_env_location = pathlib.Path.cwd().as_posix()

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Test 1: Existing envs block but no soprano env loc

        condarc_path = pathlib.Path(tmp_dir).joinpath(".condarc").as_posix()

        indentation = "  "

        partial_lines = ["envs_dirs:\n", f"{indentation}- ~/.conda/envs\n"]

        with open(condarc_path, "w") as f:
            f.writelines(partial_lines)

        expected_lines = partial_lines + [
            f"{indentation}- {mock_env_location}\n"
        ]

        setup_conda._update_condarc(mock_env_location, condarc_path)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Test 2: Existing rc but no envs block

        condarc_path = pathlib.Path(tmp_dir).joinpath(".condarc").as_posix()
        pkgs_lines = ["pkgs_dirs:\n", "  - /opt/anaconda/pkgs\n"]

        with open(condarc_path, "w") as f:
            f.writelines(pkgs_lines)

        setup_conda._update_condarc(mock_env_location, condarc_path)

        # pkgs lines should be prepended by envs lines
        expected_lines = expected_lines + pkgs_lines

        with open(condarc_path, "r") as f:
            written_lines = f.readlines()

        for expected, written in zip(expected_lines, written_lines):
            assert expected == written, (expected, written)


@pytest.mark.dependency(depends=["create_rc", "update_rc"])
def test_prepare_condarc():
    pass  # assert False


@pytest.mark.dependency(name="has_conda")
def test__has_conda():
    assert setup_conda._has_conda() or setup_conda._has_mamba()


@pytest.mark.dependency(depends=["has_conda"])
def test__build_conda_env():
    pass  # assert False
