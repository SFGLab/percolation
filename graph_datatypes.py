import numpy as np
import time

from numba import int32, jitclass


@jitclass([
    ('_parent', int32[:]),
    ('_size', int32[:]),
    ('n_clusters', int32),
    ('_track_cc', int32),
    ('cc1_root', int32),
    ('cc2_root', int32)
])
class DisjointSetPlus(object):
    def __init__(self, size: int, track_cc: int = 1):
        self._parent = np.empty(size, dtype=np.int32)
        self._size = np.empty(size, dtype=np.int32)
        self.n_clusters = -1
        self._track_cc = -1
        self.cc1_root = -1
        self.cc2_root = -1
        self.reset(track_cc)

    def reset(self, track_cc):
        self._parent.fill(0)
        self._size.fill(1)
        self._track_cc = track_cc
        if track_cc >= 1:
            self.cc1_root = 0
        if track_cc >= 2:  # TODO: > 2?
            self.cc2_root = 1
        self.n_clusters = len(self._parent)

    @property
    def track_cc(self):
        return self._track_cc

    @property
    def cc1_size(self) -> int:
        return self._size[self.cc1_root]

    @property
    def cc2_size(self) -> int:
        return 0 if self.n_clusters == 1 else self._size[self.cc2_root]  # TODO: fix?

    def find(self, x: int) -> int:
        parent = self._parent
        parent_x = parent[x]
        while parent_x != x:
            grandpa_x = parent[parent_x]
            parent[x] = grandpa_x
            x, parent_x = parent_x, grandpa_x
        return x

    def _update_cc1_and_cc2(self, root_x: int, root_y: int, new_size: int) -> None:  # size of root_x >= root_y
        _size = self._size
        if new_size >= _size[self.cc1_root]:  # cc1 is replaced
            if root_x == self.cc1_root:  # cc1 is one of the merged clusters
                if root_y == self.cc2_root:  # cc1 and cc2 merged, need to find replacement for cc2
                    _size[root_x] = -1  # temp. remove cc1 from search
                    self.cc2_root = _size.argmax()
                    _size[root_x] = new_size
                # else: cc1 grows, but cc2 remains
            elif new_size > _size[self.cc1_root]:  # cc1 is replaced, and it becomes new cc2
                self.cc2_root = self.cc1_root
                self.cc1_root = root_x
            else:  # new_size == _size[self.cc1_root]:
                # for stability in this case let's keep the old cc1, and the new cluster becomes cc2
                self.cc2_root = root_x
        elif new_size > _size[self.cc2_root]:  # cc2 is replaced
            self.cc2_root = root_x

    def union(self, x: int, y: int) -> bool:
        root_x = self.find(x)
        root_y = self.find(y)

        if root_x == root_y:
            return False

        _size = self._size
        if _size[root_x] < _size[root_y]:
            root_x, root_y = root_y, root_x
        new_size = _size[root_x] + _size[root_y]

        if self._track_cc == 2:
            self._update_cc1_and_cc2(root_x, root_y, new_size)
        elif self._track_cc == 1 and new_size > self._size[self.cc1_root]:  # cc1 is replaced
            self.cc1_root = root_x

        self._parent[root_y] = root_x
        _size[root_x] = new_size
        self.n_clusters -= 1
        return True

    def get_cluster(self, x):
        root = self.find(x)
        size = self._size[root]
        return root, size

    def get_size(self, x):
        root = self.find(x)
        return self._size[root]


class TimingMessage(object):
    def __init__(self, name: str, print_fun=print, print_on_start: str = False):
        self.name = name
        self._t0 = None
        self._print_fun = print_fun
        self._print_on_start = print_on_start
        self._already_printed = False

    def print_raw(self, verb, extra, timed=True):
        verb = verb if not timed else f'{verb} ({self.elapsed():.2f}s)'
        sep = ': ' if verb and extra else ''
        self._print_fun(f'[{self.name}] {verb}{sep}{extra}')

    def started(self, msg=''):
        self.print_raw('Started', msg, timed=False)

    def finished(self, msg=''):
        self.print_raw('Finished', msg)
        self._already_printed = True

    def exception(self, msg):
        self.print_raw('Interrupted', msg)

    def message(self, msg, verb='Running', timed=True):
        self.print_raw(verb, msg, timed=timed)

    def print(self, msg, verb='', timed=False):
        self.print_raw(verb, msg, timed=timed)

    def elapsed(self):
        return time.time() - self._t0

    def __enter__(self):
        if self._print_on_start:
            self.started()
        self._t0 = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            if not self._already_printed:
                self.finished()
        else:
            self.exception(str(exc_val))
