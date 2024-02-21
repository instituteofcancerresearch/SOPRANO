from SOPRANO.core import tidy


def _quick_build(tmp_path):
    files_to_ignore = {"foo.results.tsv", "bar.log"}
    files_to_find = {"spam.txt", "eggs.cat"}

    for filename in files_to_ignore.union(files_to_find):
        path = tmp_path / filename
        path.touch()

    class Parameters:
        # tmp override...
        cache_dir = tmp_path

    return files_to_ignore, files_to_find, Parameters()


def test_find_files_to_tar_gz(tmp_path):
    files_to_ignore, files_to_find, tmp_pars = _quick_build(tmp_path)
    detected = tidy.find_files_to_tar_gz(tmp_pars)

    for filename in files_to_find:
        assert (tmp_path / filename).as_posix() in detected

    for filename in detected:
        assert filename not in files_to_ignore


def test_tar_and_compress(tmp_path):
    files_to_ignore, files_to_find, tmp_pars = _quick_build(tmp_path)

    tidy.tar_and_compress(tmp_pars)

    assert tmp_pars.cache_dir.joinpath("intermediate.data.tar.gz").exists()

    for file in files_to_ignore:
        assert tmp_pars.cache_dir.joinpath(file).exists()

    for file in files_to_find:
        assert not tmp_pars.cache_dir.joinpath(file).exists()
