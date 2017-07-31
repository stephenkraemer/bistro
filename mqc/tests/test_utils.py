from mqc.utils import open_gzip_or_plain_file
import pytest
import gzip

class TestOpenGzipOrPlainFile:
    def test_fails_when_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            open_gzip_or_plain_file('/path/to/chrom.fa.gz')


    def test_fails_when_file_empty(self, tmpdir):
        tmp_path = tmpdir.join('chrom1.fa')
        tmp_path.write('')
        error_message = f'File {tmp_path} is empty'
        with pytest.raises(IOError) as e:
            open_gzip_or_plain_file(str(tmp_path))
        assert str(e.value) == error_message


    def test_fails_with_no_access_permissions(self, tmpdir):
        tmp_path = tmpdir.join('chrom1.fa')
        content = 'ACGC\nACGC\n'
        tmp_path.write(content)
        tmp_path.chmod(0o000)
        error_message = f"Can't open {tmp_path}"
        with pytest.raises(OSError) as e:
            open_gzip_or_plain_file(str(tmp_path))
        assert str(e.value) == error_message


    def test_recognizes_plain_text_files(self, tmpdir):
        tmp_fa_path = tmpdir.join('chrom1.fa')
        content = 'ACGC\nACGC\n'
        tmp_fa_path.write(content)
        read_content = (open_gzip_or_plain_file(str(tmp_fa_path))
                        .read())
        assert content == read_content


    def test_recognizes_gzip_reads_as_unicode_str(self, tmpdir):
        tmp_path = tmpdir.join('chrom1.fa')
        content = 'ACGC\nACGC\n'
        with gzip.open(str(tmp_path), 'wt') as fobj:
            fobj.write(content)
        read_content = (open_gzip_or_plain_file(str(tmp_path))
                        .read())
        assert content == read_content
