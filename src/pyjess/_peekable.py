import typing


class PeekableFile(typing.TextIO):
    """A buffered file that supports peeking.
    """

    def __init__(self, file):
        self.file = file
        self._buffer = ""

    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if not line:
            raise StopIteration
        return line

    def read(self, n=-1):
        l = len(self._buffer)
        if l == 0:
            output = self.file.read(n)
        elif n is None or n == -1:
            output = self._buffer + self.file.read()
            self._buffer = ""
        elif n < l:
            output = self._buffer[:n]
            self._buffer = self._buffer[n:]
        else:
            output = self._buffer + self.file.read(n - l)
            self._buffer = ""
        return output

    def readline(self):
        i = self._buffer.find("\n")
        if i == -1:
            line = self._buffer + self.file.readline()
            self._buffer = ""
        else:
            line = self._buffer[:i+1]
            self._buffer = self._buffer[i+1:]
        return line

    def peek(self, n):
        l = len(self._buffer)
        if l < n:
            self._buffer += self.file.read(n - l)
        return self._buffer[:n]

    def close(self):
        self.file.close()