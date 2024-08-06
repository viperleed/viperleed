"""Cases for test_itertool_utils.py."""


__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'


from itertools import chain


class BaseSequence:
    """Base class for the sequences below."""

    def __init__(self, seq):
        self.seq = seq
        self.ind = 0


class BaseSequenceWithIter(BaseSequence):
    """An almost-sequence with only __iter__."""

    def __iter__(self):
        return self


class BaseSequenceWithNext(BaseSequence):
    """An almost-sequence with only __next__."""

    def __next__(self):
        if self.ind >= len(self.seq):
            raise StopIteration
        item = self.seq[self.ind]
        self.ind += 1
        return item


class CasesSequence:
    """Sequence of various kind."""

    def case_generator(self):
        """A regular generator from a sequence."""
        def _make(seq):
            for item in seq:
                yield item
        return _make

    def case_get_item(self):
        """A sequence returning items via __getitem__."""
        class _GetItem(BaseSequence):
            def __getitem__(self_, index):
                return self_.seq[index]
            def __len__(self_):
                return len(self_.seq)
        return _GetItem

    def case_iterator(self):
        """A sequence that uses the iterator protocol."""
        class _Iterator(BaseSequenceWithNext, BaseSequenceWithIter):
            pass
        return _Iterator

    def case_generator_sequence(self):
        """A sequence using the iterator protocol defined with a generator."""
        generate = self.case_generator()
        class _Generator(BaseSequence):
            def __iter__(self_):
                yield from generate(self_.seq)
        return _Generator

    def case_stop_immediately(self):
        """A sequence whose __next__ stops right away."""
        class _StopImmediately(BaseSequenceWithIter):
            def __next__(self_):
                raise StopIteration
        return _StopImmediately

    def case_nested(self):
        """A combined sequence of multiple cases."""
        inner_to_outer = (self.case_get_item(),
                          self.case_iterator(),
                          self.case_generator())
        def _make(seq):
            result = seq
            for seq_type in inner_to_outer:
                result = seq_type(result)
            return chain(item for item in result)
        return _make


class CasesRaises:
    """Cases that have some implementation problems."""

    def case_misses_iterator_methods(self):
        """A class that has __next__, but no __getitem__ nor __iter__."""
        class _NoIteratorMethods(BaseSequenceWithNext):
            exc = TypeError
        return _NoIteratorMethods

    def case_misses_next(self):
        """A class with __iter__, but no __next__."""
        class _MissesNext(BaseSequenceWithIter):
            exc = TypeError
        return _MissesNext

    def case_except_immediately(self):
        """A class that raises ZeroDivisionError in __next__."""
        class _ExceptImmediately(BaseSequenceWithIter):
            exc = ZeroDivisionError
            def __next__(self):
                raise ZeroDivisionError
        return _ExceptImmediately

    def case_except_at_second(self):
        """A class that raises ZeroDivisionError at the second iteration."""
        class _ExceptAtSecond(BaseSequenceWithNext, BaseSequenceWithIter):
            exc = ZeroDivisionError
            def __next__(self):
                if self.ind == 2:
                    raise ZeroDivisionError
                return super().__next__()
        return _ExceptAtSecond
