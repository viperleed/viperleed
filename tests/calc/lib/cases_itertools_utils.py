"""Cases for test_itertool_utils.py."""


__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'


from itertools import chain


class BaseSequence:  # pylint: disable=too-few-public-methods  # Mix-in
    """Base class for the sequences below."""

    def __init__(self, seq):
        """Initialize this instance from some kind of sequence."""
        self.seq = seq
        self.ind = 0

# pylint: disable-next=too-few-public-methods                  # Mix-in
class BaseSequenceWithIter(BaseSequence):
    """An almost-sequence with only __iter__."""

    # pylint: disable-next=non-iterator-returned               # Mix-in
    def __iter__(self):
        return self


# pylint: disable-next=too-few-public-methods                  # Mix-in
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

    @staticmethod
    def case_generator():
        """A regular generator from a sequence."""
        def _make(seq):
            yield from seq
        return _make

    @staticmethod
    def case_get_item():
        """A sequence returning items via __getitem__."""
        class _GetItem(BaseSequence):
            def __getitem__(self, index):
                return self.seq[index]
            def __len__(self):
                return len(self.seq)
        return _GetItem

    @staticmethod
    def case_iterator():
        """A sequence that uses the iterator protocol."""
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _Iterator(BaseSequenceWithNext, BaseSequenceWithIter):
            pass
        return _Iterator

    @classmethod
    def case_generator_sequence(cls):
        """A sequence using the iterator protocol defined with a generator."""
        generate = cls.case_generator()
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _Generator(BaseSequence):
            def __iter__(self):
                yield from generate(self.seq)
        return _Generator

    @staticmethod
    def case_stop_immediately():
        """A sequence whose __next__ stops right away."""
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _StopImmediately(BaseSequenceWithIter):
            def __next__(self):
                raise StopIteration
        return _StopImmediately

    @classmethod
    def case_nested(cls):
        """A combined sequence of multiple cases."""
        inner_to_outer = (cls.case_get_item(),
                          cls.case_iterator(),
                          cls.case_generator())
        def _make(seq):
            result = seq
            for seq_type in inner_to_outer:
                result = seq_type(result)
            return chain(item for item in result)
        return _make


class CasesRaises:
    """Cases that have some implementation problems."""

    @staticmethod
    def case_misses_iterator_methods():
        """A class that has __next__, but no __getitem__ nor __iter__."""
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _NoIteratorMethods(BaseSequenceWithNext):
            exc = TypeError
        return _NoIteratorMethods

    @staticmethod
    def case_misses_next():
        """A class with __iter__, but no __next__."""
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _MissesNext(BaseSequenceWithIter):
            exc = TypeError
        return _MissesNext

    @staticmethod
    def case_except_immediately():
        """A class that raises ZeroDivisionError in __next__."""
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _ExceptImmediately(BaseSequenceWithIter):
            exc = ZeroDivisionError
            def __next__(self):
                raise ZeroDivisionError
        return _ExceptImmediately

    @staticmethod
    def case_except_at_second():
        """A class that raises ZeroDivisionError at the second iteration."""
        # pylint: disable-next=too-few-public-methods          # Mix-in
        class _ExceptAtSecond(BaseSequenceWithNext, BaseSequenceWithIter):
            exc = ZeroDivisionError
            def __next__(self):
                # pylint: disable-next=magic-value-comparison
                if self.ind == 2:
                    raise ZeroDivisionError
                return super().__next__()
        return _ExceptAtSecond
