import urllib.request
import os.path

cache_path = ".porch"

class FileCache:
    """
    Object handling file names
    """
    def __init__(self, file_name):
        self.file_name = file_name

    def get_file_name(self):
        return self.file_name


class UrlFileCache(FileCache):
    """
    Object handling files with urls wich are downloaded and cashed locally
    """
    def __init__(self, file_name, url):
        self.file_name = file_name
        self.url = url

    def track_dl(self):
        response = requests.get(self.url, stream=True)
        with open(self.file_name, "wb") as handle:
            for data in response.iter_content():
                handle.write(data)

    def get_file_name(self):
        if not os.path.isfile(self.file_name):
            self.track_dl()
        return self.file_name

class TsvFileTracker(FileCache):
    """
    Object handling tab separated files
    """
    def __init__(self, file_name, filler):
        self.file_name = file_name
        self.filler = filler

    def get_file_name(self):
        return self.file_name

    def read_file(self):
        if not os.path.isfile(self.get_file_name()):
            df = self.filler()
            self.write_file(df)
            return df
        return pd.read_csv(self.get_file_name(), sep="\t", index_col=0)

    def write_file(self,df):
        df.to_csv(self.file_name, sep="\t")

def download_file(path, url):
    "This function downloads a file, path, from an url, if the file is not already cached"
    if not os.path.isfile(path):
        stream = urllib.request.urlopen(url)
        with open(path,'wb') as output:
            output.write(stream.read())
    return path
