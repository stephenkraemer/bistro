import mqc.utils
import pytest
import gzip


def test_fails_when_file_not_found():
    with pytest.raises(FileNotFoundError):
        mqc.utils.open_gzip_or_plain_file('/path/to/chrom.fa.gz')


def test_fails_when_file_empty(tmpdir):
    tmp_path = tmpdir.join('chrom1.fa')
    tmp_path.write('')
    error_message = f'File {tmp_path} is empty'
    with pytest.raises(IOError) as e:
        mqc.utils.open_gzip_or_plain_file(str(tmp_path))
    assert str(e.value) == error_message


def test_fails_with_no_access_permissions(tmpdir):
    tmp_path = tmpdir.join('chrom1.fa')
    content = 'ACGC\nACGC\n'
    tmp_path.write(content)
    tmp_path.chmod(0o000)
    error_message = f"Can't open {tmp_path} as plain text file"
    with pytest.raises(OSError) as e:
        mqc.utils.open_gzip_or_plain_file(str(tmp_path))
    assert str(e.value) == error_message


def test_recognizes_plain(tmpdir):
    tmp_fa_path = tmpdir.join('chrom1.fa')
    content = 'ACGC\nACGC\n'
    tmp_fa_path.write(content)
    read_content = (mqc.utils.open_gzip_or_plain_file(str(tmp_fa_path))
                    .read())
    assert content == read_content


def test_recognizes_gzip_reads_text(tmpdir):
    tmp_path = tmpdir.join('chrom1.fa')
    content = 'ACGC\nACGC\n'
    with gzip.open(str(tmp_path), 'wt') as fobj:
        fobj.write(content)
    read_content = (mqc.utils.open_gzip_or_plain_file(str(tmp_path))
                    .read())
    assert content == read_content
